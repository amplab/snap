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

void *BigAlloc(
        size_t      sizeToAllocate,
        size_t      *sizeAllocated = NULL);

void BigDealloc(void *memory);

void *BigReserve(
        size_t      sizeToReserve,
        size_t      *sizeReserved = NULL,
        size_t      *pageSize = NULL);

bool BigCommit(
        void        *memoryToCommit,
        size_t      sizeToCommit);


