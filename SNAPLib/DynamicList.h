/*++

Module Name:

    DynamicList.h

Abstract:

    Dynamically growing list that can be efficiently appended by multiple threads.

Authors:

    Ravi Pandya, June 2012

Environment:

    User mode service.
    
    DynamicList has specific thread-safe APIs for threads to allocate blocks of space.
    DynamicListWriter provides a per-thread object that manages allocating buffers as needed.
    The elements probably will not be contigious since each thread will only fill
    part of its last block.

Revision History:

--*/

#pragma once
#include "stdafx.h"
#include "Compat.h"
#include "BigAlloc.h"
#include "DynamicMemory.h"

class DynamicList
{
public:
    // reserve and grow
    DynamicList(size_t reserveCount, unsigned growCount, unsigned i_blockSize, unsigned i_elementSize);

    // allocate fixed size, cannot grow - for pre-built data
    DynamicList(void** io_start, size_t count, unsigned i_elementSize);

    DynamicList(DynamicList& b);

    DynamicList& operator=(DynamicList& b);

    ~DynamicList();

    // base address for all elements
    inline void* getBase() { return memory.getBase(); }

    inline unsigned getElementSize() { return elementSize; }

    inline unsigned getBlockSize() { return blockSize; }

    inline size_t getAllocated() { return allocated; }

    inline DynamicMemory* getMemory() { return &memory; }

    inline void clear()
    { memset(getBase(), 0, allocated * elementSize); allocated = 0; }

private:

    friend class DynamicListWriter;

    unsigned allocate(void** o_elements, _uint32 minElements = 0);

    ExclusiveLock   lock;           // controls access to memory
    DynamicMemory   memory;         // allocates large blocks of memory
    const unsigned  elementSize;    // size of each list element
    const unsigned  blockSize;      // number of elements to allocate for each thread at a time

    size_t          allocated;      // number of elements allocated to all threads so far
};

class DynamicListWriter
{
public:
    DynamicListWriter(DynamicList* i_list);

    void* next(_uint32 nelements = 1);

    unsigned getUnused(void** o_elements);

private:
    DynamicList*    list;           // shared list
    void*           elements;       // next available element
    unsigned        count;          // count of remaining elements
};
