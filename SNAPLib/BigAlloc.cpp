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
        bool        reserveOnly = FALSE,
        size_t      *pageSize = NULL)
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
        allocatedMemory = (BYTE *)VirtualAlloc(0,largePageSizeToAllocate,commitFlag|MEM_RESERVE|MEM_LARGE_PAGES,PAGE_READWRITE);

        if (NULL != allocatedMemory) {
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
            warningPrinted= true;
            fprintf(stderr,"BigAlloc: WARNING: Unable to allocate large page memory, %d.  Falling back to VirtualAlloc.  Performance may be adversely affected.\n",GetLastError());
            if (!assertPrivilegeWorked || GetLastError() == 1314) { // TODO: Look up the error code name for 1314.
                fprintf(stderr,"BigAlloc: Unable to assert the SeLockMemoryPrivilege (%d), which is probably why it failed.\n",assertPrivilegeError);
                fprintf(stderr,"Try secpol.msc, then SecuritySettings, Local Policies, User Rights Assignment.\n");
                fprintf(stderr,"Then double click 'Lock Pages in Memory,' add the current user directly or by being\n");
                fprintf(stderr,"In a group and then reboot (you MUST reboot) for it to work.\n");
            }
        }
    }
    
    allocatedMemory = (BYTE *)VirtualAlloc(0,virtualAllocSize,commitFlag|MEM_RESERVE,PAGE_READWRITE);

    if (NULL != allocatedMemory && NULL != sizeAllocated) {
        *sizeAllocated = virtualAllocSize;
    }

    if (NULL == allocatedMemory) {
        fprintf(stderr,"BigAlloc of size %lld failed.\n",sizeToAllocate);
        exit(1);
    }

    return allocatedMemory;

}

#ifdef PROFILE_BIGALLOC

struct ProfileEntry
{
    char*   caller;
    size_t  total;
    size_t  count;
};

static const int MaxCallers = 1000;
static int NCallers = 0;
static int LastCaller = 0;
static ProfileEntry AllocProfile[1000];

void *BigAllocProfile(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated,
        char        *caller)
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
                AllocProfile[LastCaller].caller = caller;
                AllocProfile[LastCaller].total = AllocProfile[LastCaller].count = 0;
            }
        }
        if (LastCaller < MaxCallers) {
            AllocProfile[LastCaller].count++;
            AllocProfile[LastCaller].total += sizeToAllocate;
        }
    }
    return BigAllocInternal(sizeToAllocate, sizeAllocated);
}

void PrintAllocProfile()
{
    printf("BigAlloc usage\n");
    for (int i = 0; i < NCallers; i++) {
        printf("%7.1f Mb %7lld %s\n", 
            AllocProfile[i].total * 1e-6, AllocProfile[i].count, AllocProfile[i].caller);
    }
}

#else
void *BigAlloc(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated)
{
    return BigAllocInternal(sizeToAllocate, sizeAllocated);
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
    VirtualFree(memory,0,MEM_RELEASE);
}

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
        printf("BigCommit VirtualAlloc failed with error 0x%x\n", GetLastError());
    }
    return allocatedMemory != NULL;
}

#else /* no _MSC_VER */

void *BigAlloc(
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
        exit(1);
    }

#if (defined(MADV_HUGEPAGE) && !defined(USE_HUGETLB))
    // Tell Linux to use huge pages for this range
    if (madvise(mem, sizeToAllocate, MADV_HUGEPAGE) == -1) {
        fprintf(stderr, "WARNING: failed to enable huge pages -- your kernel may not support it\n"); 
    }
#endif

    // Remember the size allocated in the first sizeof(size_t) bytes
    *((size_t *) mem) = sizeToAllocate;
    return (void *) (mem + sizeof(size_t));
}


void BigDealloc(void *memory)
{
    // Figure out the size we had allocated
    char *startAddress = ((char *) memory) - sizeof(size_t);
    size_t sizeAllocated = *((size_t *) startAddress);
    if (munmap(startAddress, sizeAllocated) != 0) {
        perror("munmap");
        exit(1);
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
    return BigAlloc(sizeToReserve, sizeReserved);
}

bool BigCommit(
        void        *memoryToCommit,
        size_t      sizeToCommit)
{
    // TODO: use actual reserve/commit API; this is a temporary hack
    return true;
}



#endif /* _MSC_VER */
