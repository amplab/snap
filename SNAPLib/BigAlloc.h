/*++

Module Name:

    bigalloc.h

Abstract:

    Headers for an allocator that uses big pages where appropriate and possible.

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

Revision History:

--*/

#pragma once

inline unsigned RoundUpToPageSize(unsigned size)
{
    const unsigned pageSize = 4096;
    return ((size + pageSize - 1) / pageSize) * pageSize;
}

//#define PROFILE_BIGALLOC

#ifdef PROFILE_BIGALLOC

#define BigAlloc(s) BigAllocProfile((s), NULL, __FUNCTION__)
#define BigAlloc2(s,p) BigAllocProfile((s), (p), __FUNCTION__)
#define BigReserve(s) BigReserveProfile((s), NULL, NULL, __FUNCTION__)
#define BigCommit(p, s) BigCommitProfile((p), (s), __FUNCTION__)

void *BigAllocProfile(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated = NULL,
        const char* caller = NULL);

void *BigReserveProfile(
    size_t      sizeToReserve,
    size_t      *sizeReserved = NULL,
    size_t      *pageSize = NULL,
    const char* caller = NULL);

bool BigCommitProfile(
    void        *memoryToCommit,
    size_t      sizeToCommit,
    const char* caller = NULL);

#else

void *BigAlloc(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated = NULL);
#define BigAlloc2(s,p) BigAlloc((s), (p))

void *BigReserve(
    size_t      sizeToReserve,
    size_t      *sizeReserved = NULL,
    size_t      *pageSize = NULL);

bool BigCommit(
    void        *memoryToCommit,
    size_t      sizeToCommit);

#endif

void PrintBigAllocProfile();

void BigDealloc(void *memory);

//
// This class is used to allocate a group of objects all onto a single set of big pages.  It requires knowing
// the amount of memory to be allocated when it's created.  It does not support deleting memory other than
// all at once.
//
class BigAllocator {
public:
    BigAllocator(size_t i_maxMemory, size_t i_allocationGranularity = 8);
    ~BigAllocator();

    virtual void *allocate(size_t amountToAllocate);

#if     _DEBUG
    void checkCanaries();
#else  // DEBUG
    void checkCanaries() {}
#endif  // DEBUG
private:

    char    *basePointer;
    char    *allocPointer;
    size_t  maxMemory;
    size_t  allocationGranularity;

#if     _DEBUG
    //
    // Stick a canary between each allocation and 
    unsigned    nCanaries;
    static const unsigned maxCanaries = 100;
    static const unsigned canaryValue = 0xca4a71e5;
    unsigned    *canaries[maxCanaries];
#endif  // DEBUG
};

//
// An allocator that doesn't actually allocate, it just counts how much it would allocate.  The idea is that
// you can write allocations in a fairly normal looking way, call them with this to see how much would be
// allocated, then create a real BigAllocator with that amount of memory.  That way, you don't need to
// keep in sync the actual allocation and the code that knows how much memory will be needed.
//
class CountingBigAllocator : public BigAllocator
{
public:
    CountingBigAllocator(size_t i_allocationGranularity = 8) :size(0), allocations(NULL), BigAllocator(0), allocationGranularity(i_allocationGranularity) {}
    ~CountingBigAllocator();

    virtual void *allocate(size_t amountToAllocate);
    virtual void assertAllMemoryUsed() {}
    size_t getMemoryUsed() {return size;}

private:
    size_t  size;
    size_t  allocationGranularity;

    struct Allocation {
        void *ptr;
        Allocation *next;
    } *allocations;
};

extern bool BigAllocUseHugePages;


// trivial per-thread heap for use in zalloc
struct ThreadHeap
{
    char* start;
    char* end;
    char* next;
    ThreadHeap(size_t bytes)
    {
        next = start = (char*) BigAlloc(bytes);
        end = start + bytes;
    }
    void* alloc(size_t bytes)
    {
        if (next + bytes <= end) {
            void* result = next;
            next += bytes;
            return result;
        }
        return NULL;
    }
    bool free(void* p)
    {
        return (char*)p >= start && (char*) p <= end;
    }
    void reset()
    {
        next = start;
    }
    ~ThreadHeap()
    {
        BigDealloc(start);
    }
};

void* zalloc(void* opaque, unsigned items, unsigned size);

void zfree(void* opaque, void* p);
