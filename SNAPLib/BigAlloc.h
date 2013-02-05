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

//#define PROFILE_BIGALLOC

#ifdef PROFILE_BIGALLOC
#define BigAlloc(s) BigAllocProfile((s), NULL, __FUNCTION__)

void *BigAllocProfile(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated = NULL,
        char*       caller = NULL);

void PrintAllocProfile();

#else
void *BigAlloc(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated = NULL);
#endif

void BigDealloc(void *memory);

void *BigReserve(
        size_t      sizeToReserve,
        size_t      *sizeReserved = NULL,
        size_t      *pageSize = NULL);

bool BigCommit(
        void        *memoryToCommit,
        size_t      sizeToCommit);

//
// This class is used to allocate a group of objects all onto a single set of big pages.  It requires knowing
// the amount of memory to be allocated when it's created.  It does not support deleting memory other than
// all at once.
//
class BigAllocator {
public:
    BigAllocator(size_t i_maxMemory);
    ~BigAllocator();

    void *allocate(size_t amountToAllocate);
    void assertAllMemoryUsed();

#if     _DEBUG
    void checkCanaries();
#else  // DEBUG
    void checkCanaries() {}
#endif  // DEBUG
private:

    char    *basePointer;
    char    *allocPointer;
    size_t  maxMemory;

#if     _DEBUG
    //
    // Stick a canary between each allocation and 
    unsigned    nCanaries;
    static const unsigned maxCanaries = 100;
    static const unsigned canaryValue = 0xca4a71e5;
    unsigned    *canaries[maxCanaries];
#endif  // DEBUG
};


