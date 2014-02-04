/*++

Module Name:

    DynamicList.cpp

Abstract:

    Dynamically growing list that can be efficiently appended by multiple threads.

Authors:

    Ravi Pandya, June 2012

Environment:

Revision History:

--*/

#include "stdafx.h"
#include "Compat.h"
#include "exit.h"
#include "BigAlloc.h"
#include "DynamicMemory.h"
#include "DynamicList.h"

using std::max;

DynamicList::DynamicList(
    size_t reserveCount,
    unsigned growCount,
    unsigned i_blockSize,
    unsigned i_elementSize)
    : memory(max(reserveCount, (size_t) growCount) * i_elementSize, growCount * (size_t) i_elementSize),
    blockSize(i_blockSize),
    elementSize(i_elementSize),
    allocated(0)
{
    InitializeExclusiveLock(&lock);
}

DynamicList::DynamicList(
    void** io_start,
    size_t count,
    unsigned i_elementSize)
    : memory(*io_start, count * elementSize),
    elementSize(i_elementSize),
    blockSize(0),
    allocated(count)
{
    *io_start = ((char*) *io_start) + count * elementSize;
}

DynamicList::DynamicList(DynamicList& b)
    : memory(b.memory),
    elementSize(b.elementSize),
    blockSize(b.blockSize),
    allocated(b.allocated)
{
    InitializeExclusiveLock(&lock);
}

DynamicList& DynamicList::operator=(DynamicList& b)
{
    memory = b.memory;
    _ASSERT(elementSize == b.elementSize && blockSize == b.blockSize);
    allocated = b.allocated;
    return *this;
}

DynamicList::~DynamicList()
{
    DestroyExclusiveLock(&lock);
}

    unsigned
DynamicList::allocate(
    void** o_elements,
    _uint32 minElements)
{
    AcquireExclusiveLock(&lock);
    
    // try growing array, ensure it succeeds
    _uint32 block = max(blockSize, minElements);
    for (int i = 0; i < 2; i++) {
        if (elementSize * (size_t) (allocated + block) > memory.getCommitted()) {
            if (i == 0) {
                memory.grow();
            } else {
                fprintf(stderr, "DynamicList: unable to allocate more elements\n");
				soft_exit(1);
            }
        }
    }
    
    // return pointer to new array segment
    *o_elements = ((char*) memory.getBase()) + elementSize * allocated;

    // figure out increase
    size_t available = memory.getCommitted() / elementSize;
    _ASSERT(available >= allocated);
    unsigned increase = (unsigned) min((size_t) blockSize, available - allocated);
    allocated += increase;
    
    ReleaseExclusiveLock(&lock);
    return increase;
}

DynamicListWriter::DynamicListWriter(
    DynamicList* i_list)
    : list(i_list)
{
    count = list->allocate(&elements);
}

    void*
DynamicListWriter::next(
    _uint32 nelements)
{
    if (count < nelements) {
        count = list->allocate(&elements, nelements);
        if (count == 0) {
            printf("cannot allocate more elements\n");
            soft_exit(1);
        }
    }
    void* result = elements;
    elements = ((char*) elements) + list->getElementSize() * nelements;
    count -= nelements;
    return result;
}
    unsigned
DynamicListWriter::getUnused(
    void** o_elements)
{
    *o_elements = elements;
    return count;
}
