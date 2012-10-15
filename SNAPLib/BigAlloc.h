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


