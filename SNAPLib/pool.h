/*++

Module Name:

	pool.h

Abstract:

	Fixed size memory allocator headers.

Author:

	Bill Bolosky		[bolosky]		1993

Revision History:

    7/30/2015 - taken from SRS/Farsite (bolosky)

--*/
#pragma once

struct PoolEntry;
struct PoolBlob;

class Pool {
public:
                Pool(unsigned objectSize);

				~Pool();

    void		preAllocate(
			    unsigned		 					n);

    void		*allocate(void);

    void	 	free(
				    void							*object);

private:

    PoolEntry		*getEntry(void);

    void		 releaseEntry(
			    PoolEntry		*entry);

    void		 allocateMoreObjects(void);

    unsigned		 objectSize;
    struct PoolEntry	*entries;		// PoolEntries with vaid data attached to them
    struct PoolEntry	*freeEntries;		// PoolEntries without valid data attached to them
    struct PoolBlob	*entriesBlobHead;	// The head of the blob list for PoolEntries
    unsigned		 entriesPerBlob;	// How many entries in an entry blob
    struct PoolBlob	*objectsBlobHead;	// The head of the blob list for the allocated objects
    unsigned		 objectsPerBlob;	// How many objects in an object blob

    unsigned		 numFree;		// Current size of free list

	//
	// Is this pool one in which we overlap the PoolEntries (ie., free list) and the pool objects.
	// Doing this means that we're not using extra memory for pool entries, but it's not always possible.
	// If we're using a provided allocator, then we can't mess with the memory that we get back, and if
	// we've got objects that are too small to hold PoolEntries, then it doesn't work, either.
	//
	bool			 overlapEntriesAndObjects;

};

