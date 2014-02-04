/*++

Module Name:

    DynamicMemory.cpp

Abstract:

    Dynamically growable allocator for a large memory block

Authors:

    Ravi Pandya, May 2012

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

--*/

#include "stdafx.h"
#include "BaseAligner.h"
#include "Compat.h"
#include "BigAlloc.h"
#include "DynamicMemory.h"

using std::min;

DynamicMemory::DynamicMemory(
    size_t i_reserve,
    size_t i_growSize)
{
    _ASSERT(i_reserve >= i_growSize && i_growSize > 0);
    size_t pageSize;
    base = BigReserve(i_reserve, &reserved, &pageSize);
    growSize = pageSize * ((i_growSize + pageSize - 1) / pageSize);
    committed = 0;
    grow();
}

DynamicMemory::DynamicMemory(
    void* i_base,
    size_t i_size)
{
    base = i_base;
    reserved = committed = i_size;
    growSize = 0;
}

DynamicMemory::DynamicMemory(DynamicMemory& b)
{
    if (growSize > 0) {
        BigDealloc(base);
    }
    base = b.base;
    reserved = b.reserved;
    committed = b.committed;
    growSize = b.growSize;
    b.growSize = 0; // don't free memory when b is deleted
}

DynamicMemory& DynamicMemory::operator=(DynamicMemory& b)
{
    if (growSize > 0) {
        BigDealloc(base);
    }
    base = b.base;
    reserved = b.reserved;
    committed = b.committed;
    growSize = b.growSize;
    b.growSize = 0; // don't free memory when b is deleted
    return *this;
}

DynamicMemory::~DynamicMemory()
{
    if (growSize > 0) {
        BigDealloc(base);
    }
}

    bool
DynamicMemory::grow()
{
    size_t growth =  min(growSize, reserved - committed);
    if (growth == 0) {
        return false;
    }
    bool ok = BigCommit(((char*)base) + committed, growth);
    if (ok) {
        committed += growth;
    }
    return ok;
}

