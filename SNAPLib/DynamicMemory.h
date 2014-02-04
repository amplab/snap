/*++

Module Name:

    DynamicMemory.h

Abstract:

    Dynamically growing large memory allocator. Reserves a maximum-size block of address space,
    and then commits chunks of it as needed.

Authors:

    Ravi Pandya, May 2012

Environment:

    User mode service.
    
    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

--*/

#pragma once
#include "stdafx.h"
#include "Compat.h"
#include "BigAlloc.h"

using std::min;

struct DynamicMemory
{
public:

    // reserve maximum space, set size to grow at each call
    DynamicMemory(size_t i_reserve, size_t i_growSize);

    // allocate fixed region of memory, cannot grow
    DynamicMemory(void* i_base, size_t i_size);

    DynamicMemory(DynamicMemory& b);

    DynamicMemory& operator=(DynamicMemory& b);

    ~DynamicMemory();

    // base address of memory block
    inline void* getBase() { return base; }

    // number of bytes committed so far
    inline size_t getCommitted() { return committed; }

    // grow by growSize, return true if succeeded
    // NOT threadsafe, must be called in an exclusive lock
    bool grow();

private:

    // base address
    void* base;

    // max bytes reserved
    size_t reserved;

    // bytes to commit per call to grow(); 0 if pre-allocated memory
    size_t growSize;

    // bytes committed so far
    size_t committed;
};
