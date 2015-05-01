/*++

Module Name:

    bigalloc.cpp

Abstract:

    Allocator that uses big pages where appropriate and possible.

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

Revision History:

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "exit.h"
#include "Error.h"

bool BigAllocUseHugePages = false;


#ifdef PROFILE_BIGALLOC

struct ProfileEntry
{
    ProfileEntry() : caller(NULL), total(0), count(0) {}
    const char*   caller;
    size_t  total;
    size_t  count;
};

static const int MaxCallers = 1000;
static int NCallers = 0;
static int LastCaller = 0;
static ProfileEntry AllocProfile[1000];

static ProfileEntry ProfileTotal;
static ProfileEntry LastPrintProfile;

void *BigAllocInternal(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated,
        bool        reserveOnly = FALSE,
        size_t      *pageSize = NULL);

void RecordAllocProfile(size_t bytes, const char* caller)
{
    if (caller) {
        if (LastCaller >= NCallers || strcmp(AllocProfile[LastCaller].caller, caller)) {
            LastCaller = NCallers;
            for (int i = 0; i < NCallers; i++) {
                if (0 == strcmp(AllocProfile[i].caller, caller)) {
                    LastCaller = i;
                    break;
                }
            }
            if (LastCaller == NCallers && NCallers < MaxCallers) {
                NCallers++;
                char* buffer = (char*) malloc(strlen(caller) + 1);
                strcpy(buffer, caller);
                AllocProfile[LastCaller].caller = buffer;
                AllocProfile[LastCaller].total = AllocProfile[LastCaller].count = 0;
            }
        }
        if (LastCaller < MaxCallers) {
            AllocProfile[LastCaller].count++;
            AllocProfile[LastCaller].total += bytes;
        }
    }
    ProfileTotal.count++;
    ProfileTotal.total += bytes;
    if (ProfileTotal.count - LastPrintProfile.count >= 1000 || ProfileTotal.total - LastPrintProfile.total >= ((size_t)1 << 30)) {
        fprintf(stderr, "BigAllocProfile %lld allocs, %lld total; caller %s alloc %lld\n", ProfileTotal.count, ProfileTotal.total, caller ? caller : "?", bytes);
        LastPrintProfile = ProfileTotal;
    }
}

void *BigAllocProfile(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated,
        const char  *caller)
{
    RecordAllocProfile(sizeToAllocate, caller);
    return BigAllocInternal(sizeToAllocate, sizeAllocated);
}

#endif

#ifdef _MSC_VER

//
// Assert an NT privilege for this thread.
//
BOOL
AssertPrivilege(
    IN  LPCSTR PrivilegeName
    )

{
    BOOL                b;
    HANDLE              hThread;
    HANDLE              hProcess;
    TOKEN_PRIVILEGES    tokenPrivileges, oldTokenPrivileges;
    DWORD               oldPrivilegesLength;


    b = OpenThreadToken(GetCurrentThread(), TOKEN_ADJUST_PRIVILEGES |
                        TOKEN_QUERY, TRUE, &hThread);
    if (!b) {
        if (GetLastError() != ERROR_NO_TOKEN) {
            return b;
        }

        b = OpenProcessToken(GetCurrentProcess(), TOKEN_DUPLICATE, &hProcess);
        if (!b) {
            return b;
        }

        b = DuplicateTokenEx(hProcess, TOKEN_ADJUST_PRIVILEGES | TOKEN_QUERY |
                             TOKEN_IMPERSONATE, NULL, SecurityImpersonation,
                             TokenImpersonation, &hThread);
        if (!b) {
            CloseHandle(hProcess);
            return b;
        }

        b = SetThreadToken(NULL, hThread);
        if (!b) {
            CloseHandle(hProcess);
            CloseHandle(hThread);
            return b;
        }

        CloseHandle(hProcess);
    }

    ZeroMemory(&tokenPrivileges, sizeof(tokenPrivileges));

    b = LookupPrivilegeValue(NULL,PrivilegeName,&tokenPrivileges.Privileges[0].Luid);
    if (!b) {
        return b;
    }

    tokenPrivileges.PrivilegeCount = 1;
    tokenPrivileges.Privileges[0].Attributes = SE_PRIVILEGE_ENABLED;

    b = AdjustTokenPrivileges(hThread, FALSE, &tokenPrivileges,
                              sizeof(tokenPrivileges), &oldTokenPrivileges,
                              &oldPrivilegesLength);

    CloseHandle(hThread);

    return b;
}


void *BigAllocInternal(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated,
        bool        reserveOnly,
        size_t      *pageSize)
/*++

Routine Description:

    Allocate memory, using large pages if both appropriate and possible, and always using
    VirtualAlloc (meaning that this will always use at least one VM page, so you shouldn't
    use it for small stuff, only gigantic data structures for which you want to reduce TLB
    misses and cache misses on the page table).  Use malloc or new for ordinary allocations.

Arguments:

    sizeToAllocate      - The amount of memory that is needed
    sizeAllocated       - Optional parameter that if provided returns the amount of memory actually allocated, which
                          will always be >= sizeToAllocate (unless the allocation fails).
    reserveOnly         - If TRUE, will only reserve address space, must call BigCommit to commit memory
    pageSize            - Optional parameter that if provided returns the page size (not large page size)

Return Value:

    pointer to the memory allocated, or NULL if the allocation failed.

--*/
{
    if (sizeToAllocate == 0) {
       sizeToAllocate = 1;
    }

    static bool warningPrinted = false;

    void *allocatedMemory;

    SYSTEM_INFO systemInfo[1];
    GetSystemInfo(systemInfo);

    size_t virtualAllocSize = ((sizeToAllocate + systemInfo->dwPageSize - 1) / systemInfo->dwPageSize) * systemInfo->dwPageSize;
    if (pageSize != NULL) {
        *pageSize = systemInfo->dwPageSize;
    }

    //
    // Try to do the VirtualAlloc using large pages if the size we're getting is at last one large page.
    // Callers should have asserted the SeLockMemoryPrivilege if they want large pages.
    //

    size_t largePageSize = GetLargePageMinimum();
    DWORD commitFlag = reserveOnly ? 0 : MEM_COMMIT;
    if (0 != largePageSize && virtualAllocSize >= largePageSize) {
        //
        // Start by asserting the SeLockMemoryPrivilege, which is necessary for large page allocations.  It's overkill to
        // do this every time, it only has to happen once/thread.  However, a BigAllocation is a big deal and shouldn't be
        // happening very much, so we just don't worry about the extra cost.
        //
        BOOL assertPrivilegeWorked = AssertPrivilege("SeLockMemoryPrivilege");
        DWORD assertPrivilegeError = GetLastError();

        size_t largePageSizeToAllocate = ((virtualAllocSize + largePageSize - 1) / largePageSize) * largePageSize;

#if     _DEBUG
        largePageSizeToAllocate += largePageSize;   // For the guard page.
#endif  // DEBUG

        allocatedMemory = (BYTE *)VirtualAlloc(0,largePageSizeToAllocate,commitFlag|MEM_RESERVE|((BigAllocUseHugePages && !reserveOnly) ? MEM_LARGE_PAGES : 0),PAGE_READWRITE);

        if (NULL != allocatedMemory) {
#if     _DEBUG
            DWORD oldProtect;
            if (!VirtualProtect((char *)allocatedMemory + virtualAllocSize, systemInfo->dwPageSize, PAGE_NOACCESS, &oldProtect)) {
                static bool printedVirtualProtectedWarning = false;
                if (! printedVirtualProtectedWarning) {
                    //WriteErrorMessage("VirtualProtect for guard page failed, %d\n", GetLastError());
                    printedVirtualProtectedWarning = true;
                }
            }
            largePageSizeToAllocate -= largePageSize;   // Back out the guard page
#endif  // DEBUG
            if (NULL != sizeAllocated) {
                *sizeAllocated = largePageSizeToAllocate;
            }
            return allocatedMemory;
        } else if (!warningPrinted) {
            //
            // The first time we fail, print out a warning and then fall back to VirtualAlloc.  We want be able to use
            // the fallback because the caller might not be able to assert the appropriate privilege and we'd still like
            // to run.  The check for printing only once isn't thread safe, so you might get more than one printed
            // if multiple threads fail at the same time.
            //
            warningPrinted = true;
            WriteErrorMessage("BigAlloc: WARNING: Unable to allocate large page memory, %d.  Falling back to VirtualAlloc.  Performance may be adversely affected.  Size = %lld\n", GetLastError(), largePageSizeToAllocate);
            if (!assertPrivilegeWorked || GetLastError() == 1314) { // TODO: Look up the error code name for 1314.
                WriteErrorMessage("BigAlloc: Unable to assert the SeLockMemoryPrivilege (%d), which is probably why it failed.\n"
                                  "Try secpol.msc, then SecuritySettings, Local Policies, User Rights Assignment.\n"
                                  "Then double click 'Lock Pages in Memory,' add the current user directly or by being\n"
                                  "In a group and then reboot (you MUST reboot) for it to work.\n", GetLastError());
            }
        }
    }
    
    allocatedMemory = (BYTE *)VirtualAlloc(0,virtualAllocSize,commitFlag|MEM_RESERVE,PAGE_READWRITE);

    if (NULL != allocatedMemory && NULL != sizeAllocated) {
        *sizeAllocated = virtualAllocSize;
    }

    if (NULL == allocatedMemory) {
        WriteErrorMessage("BigAlloc of size %lld failed.\n", sizeToAllocate);
#ifdef PROFILE_BIGALLOC
        PrintBigAllocProfile();
#endif
        soft_exit(1);
    }

    return allocatedMemory;

}

#ifndef PROFILE_BIGALLOC
void *BigAlloc(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated)
{
    return BigAllocInternal(sizeToAllocate, sizeAllocated, FALSE, NULL);
}
#endif

void BigDealloc(void *memory)
/*++

Routine Description:

    Free memory allocated by BigAlloc.

Arguments:

    memory  - address of the memory to free.


--*/
{
    if (NULL == memory) return;
    VirtualFree(memory,0,MEM_RELEASE);
}

#ifdef PROFILE_BIGALLOC
void *BigReserveProfile(
    size_t      sizeToReserve,
    size_t      *sizeReserved,
    size_t      *pageSize,
    const char* caller)
{
    char buffer[1000];
    strncpy(buffer, caller, sizeof(buffer));
    strncat(buffer, "(RESERVE)", sizeof(buffer));
    RecordAllocProfile(sizeToReserve, buffer);
    return BigAllocInternal(sizeToReserve, sizeReserved, TRUE, pageSize);
}

bool BigCommitProfile(
    void        *memoryToCommit,
    size_t      sizeToCommit,
    const char* caller)
{
    char buffer[1000];
    strncpy(buffer, caller, sizeof(buffer));
    strncat(buffer, "(COMMIT)", sizeof(buffer));
    RecordAllocProfile(sizeToCommit, buffer);
    void* allocatedMemory = VirtualAlloc(memoryToCommit, sizeToCommit, MEM_COMMIT, PAGE_READWRITE);
    if (allocatedMemory == NULL) {
        WriteErrorMessage("BigCommit VirtualAlloc failed with error 0x%x\n", GetLastError());
    }
    return allocatedMemory != NULL;
}
#else
void *BigReserve(
        size_t      sizeToReserve,
        size_t      *sizeReserved,
        size_t      *pageSize)
{
    return BigAllocInternal(sizeToReserve, sizeReserved, TRUE, pageSize);
}

bool BigCommit(
    void        *memoryToCommit,
    size_t      sizeToCommit)
{
    void* allocatedMemory = VirtualAlloc(memoryToCommit, sizeToCommit, MEM_COMMIT, PAGE_READWRITE);
    if (allocatedMemory == NULL) {
        WriteErrorMessage("BigCommit VirtualAlloc failed with error 0x%x\n", GetLastError());
    }
    return allocatedMemory != NULL;
}
#endif


#else /* no _MSC_VER */

#ifdef PROFILE_BIGALLOC
void *BigAllocInternal(
#else
void *BigAlloc(
#endif
        size_t      sizeToAllocate,
        size_t      *sizeAllocated)
{
    // Make space to include the allocated size at the start of our region; this is necessary
    // so that we can BigDealloc the memory later.
    sizeToAllocate += sizeof(size_t);

    const size_t ALIGN_SIZE = 4096;
    if (sizeToAllocate % ALIGN_SIZE != 0) {
        sizeToAllocate += ALIGN_SIZE - (sizeToAllocate % ALIGN_SIZE);
    }
    if (sizeAllocated != NULL) {
      *sizeAllocated = sizeToAllocate - sizeof(size_t);
    }

    int flags = MAP_PRIVATE|MAP_ANONYMOUS;
#ifdef USE_HUGETLB
    flags |= MAP_HUGETLB;
#endif
    char *mem = (char *) mmap(NULL, sizeToAllocate, PROT_READ|PROT_WRITE, flags, -1, 0);
    if (mem == MAP_FAILED) {
        perror("mmap");
        soft_exit(1);
    }

#if (defined(MADV_HUGEPAGE) && !defined(USE_HUGETLB))
    // Tell Linux to use huge pages for this range
    if (BigAllocUseHugePages) {
        if (madvise(mem, sizeToAllocate, MADV_HUGEPAGE) == -1) {
            WriteErrorMessage("WARNING: failed to enable huge pages -- your kernel may not support it\n"); 
        }
    }
#endif

    // Remember the size allocated in the first sizeof(size_t) bytes
    *((size_t *) mem) = sizeToAllocate;
    return (void *) (mem + sizeof(size_t));
}


void BigDealloc(void *memory)
{
    if (NULL == memory) return;
    // Figure out the size we had allocated
    char *startAddress = ((char *) memory) - sizeof(size_t);
    size_t sizeAllocated = *((size_t *) startAddress);
    if (munmap(startAddress, sizeAllocated) != 0) {
        perror("munmap");
        soft_exit(1);
    }
}

void *BigReserve(
        size_t      sizeToReserve,
        size_t      *sizeReserved,
        size_t      *pageSize)
{
    // TODO: use actual reserve/commit API; this is a temporary hack
    if (pageSize != NULL) {
        *pageSize = 4096;
    }
#ifdef PROFILE_BIGALLOC
    return BigAllocInternal(sizeToReserve, sizeReserved);
#else
    return BigAlloc(sizeToReserve, sizeReserved);
#endif
}

bool BigCommit(
        void        *memoryToCommit,
        size_t      sizeToCommit)
{
    // TODO: use actual reserve/commit API; this is a temporary hack
    return true;
}



#endif /* _MSC_VER */

BigAllocator::BigAllocator(size_t i_maxMemory, size_t i_allocationGranularity) : maxMemory(i_maxMemory), allocationGranularity(i_allocationGranularity)
{
#if     _DEBUG
    maxMemory += maxCanaries * sizeof(unsigned);
#endif  // DEBUG
    basePointer = (char *)BigAlloc(__max(maxMemory, 2 * 1024 * 1024)); // The 2MB minimum is to assure this lands in a big page
    allocPointer = basePointer;

#if     _DEBUG
    //
    // Stick a canary at the beginning of the array so that we can detect underflows for whatever's allocated first.
    //
    canaries[0] = (unsigned *) allocPointer;
    *canaries[0] = canaryValue;
    nCanaries = 1;
    allocPointer += sizeof(unsigned);
#endif  // _DEBUG
}

BigAllocator::~BigAllocator()
{
    BigDealloc(basePointer);
}

void *
BigAllocator::allocate(size_t amountToAllocate)
{
    //
    // Round up to the allocation granularity.
    //
    if ((size_t)allocPointer % allocationGranularity != 0) {
        allocPointer = (char *)((size_t)allocPointer + allocationGranularity - (size_t)allocPointer % allocationGranularity);
        _ASSERT((size_t)allocPointer % allocationGranularity == 0);
    }

    if (allocPointer + amountToAllocate > basePointer + maxMemory) {
        WriteErrorMessage("BigAllocator: allocating too much memory, %lld > %lld\n", allocPointer + amountToAllocate  - basePointer , maxMemory);
        soft_exit(1);
    }
 
    void *retVal = allocPointer;
    allocPointer += amountToAllocate;

#if     _DEBUG
    if (nCanaries < maxCanaries) {
        _ASSERT(allocPointer + sizeof(unsigned) <= basePointer + maxMemory);
        canaries[nCanaries] = (unsigned *)allocPointer;
        *canaries[nCanaries] = canaryValue;
        nCanaries++;
        allocPointer += sizeof(unsigned);
    }
#endif  // DEBUG
    return retVal;
}
 
#if     _DEBUG
    void
BigAllocator::checkCanaries()
{
    bool allOK = true;
    for (unsigned i = 0; i < nCanaries; i++) {
        if (*canaries[i] != canaryValue) {
            WriteErrorMessage("Memory corruption detected: canary at 0x%llx has value 0x%llx\n",canaries[i], *canaries[i]);
            allOK = false;
        }
    }
    _ASSERT(allOK);
}
#endif  // DEBUG


    void *
CountingBigAllocator::allocate(size_t sizeToAllocate)
{            
    size += sizeToAllocate + allocationGranularity - 1; // Add in the max roundoff

    Allocation *allocation = new Allocation;
    allocation->next = allocations;
    allocation->ptr = malloc(sizeToAllocate);
    allocations = allocation;
    return allocation->ptr;
}

CountingBigAllocator::~CountingBigAllocator()
{
    while (NULL != allocations) {
        Allocation *allocation = allocations;
        allocations = allocation->next;
        free(allocation->ptr);
        delete allocation;
    }
}

void PrintBigAllocProfile()
{
#ifdef PROFILE_BIGALLOC
    WriteStatusMessage("BigAlloc usage\n");
    for (int i = 0; i < NCallers; i++) {
        WriteStatusMessage("%7.1f Mb %7lld %s\n", 
            AllocProfile[i].total * 1e-6, AllocProfile[i].count, AllocProfile[i].caller);
    }
#endif
}

void* zalloc(void* opaque, unsigned items, unsigned size)
{
    size_t bytes = items * (size_t) size;
    void* result = ((ThreadHeap*) opaque)->alloc(bytes);
    static int printed = 0;
    if ((! result) && printed++ < 10) {
        WriteErrorMessage("warning: zalloc using malloc for %lld bytes\n", bytes);
    }
    return result ? result : malloc(bytes);
}

void zfree(void* opaque, void* p)
{
    if (! ((ThreadHeap*) opaque)->free(p)) {
        free(p);
    }
}
