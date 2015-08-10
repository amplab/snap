/*++

Module Name:

	pool.cpp

Abstract:

	Fixed size memory allocator.

Author:

	Bill Bolosky		[bolosky]		1993

Revision History:

    7/30/2015 - adapted from SRS/Farsite (bolosky)

--*/

#include "stdafx.h"
#include "pool.h"
#include "BigAlloc.h"


struct PoolEntry {
    void		*object;
    struct PoolEntry	*next;
};

struct PoolBlob {
    struct PoolBlob	*next;
    void		*data;
};


// This version of the pool constructor uses the old kind of allocator function that only returns one object.  Object blobs have one object,
// and the blob size for entry blobs is smaller so that we have finer grain memory allocation.
Pool::Pool(
    unsigned		 objectSize)
{
    _ASSERT(objectSize > 0);

    this->objectSize = objectSize;
    entries = NULL;
    freeEntries = NULL;
    entriesBlobHead = NULL;
    objectsBlobHead = NULL;

    unsigned blobSize = 1024 - 50;	// Our default allocation size; we leave the 50 byte headroom for the underlying allocator

    entriesPerBlob = blobSize / sizeof(PoolEntry);
    _ASSERT(entriesPerBlob > 0);

    objectsPerBlob = 1;

    numFree = 0;

	overlapEntriesAndObjects = objectSize >= sizeof(PoolEntry);
}

Pool::~Pool(void)
{
    // Just delete the blob lists.  All objects that have been allocated from this pool will be destroyed.
    
    while (entriesBlobHead) {
        _ASSERT(!overlapEntriesAndObjects);

		PoolBlob *blob = entriesBlobHead;
        _ASSERT(blob->data);
        BigDealloc(blob->data);
		entriesBlobHead = blob->next;
        delete blob;
    }

    while (objectsBlobHead) {
		PoolBlob *blob = objectsBlobHead;
        _ASSERT(blob->data);

        delete [] blob->data;

		objectsBlobHead = blob->next;
        delete blob;
    }

}


    void
Pool::allocateMoreObjects(void)
{
    _ASSERT(objectsPerBlob);

    PoolBlob *blob = new PoolBlob;

    blob->data = BigAlloc(objectSize * objectsPerBlob);

    if (!blob->data) {
		delete blob;
		return;
    }

    blob->next = objectsBlobHead;
    objectsBlobHead = blob;

    // Now put them on the free list.

	if (overlapEntriesAndObjects) {
		//
		// Turn the new objects into entries that point at themselves.
		//
		for (unsigned i = 0; i < objectsPerBlob; i++) {
			PoolEntry *entry = (PoolEntry *)((char *)blob->data + i * objectSize);
			entry->object = entry;
			entry->next = entries;
			entries = entry;
		}

	} else {
	    for (unsigned i = 0; i < objectsPerBlob; i++) {
    	    PoolEntry *entry = getEntry();
			if (!entry) {
			    return;		// This is kinda bogus, because it might leave some allocated objects unreachable.
			}
			entry->object = (void *)(((char *)blob->data) + i * objectSize);
			entry->next = entries;
			entries = entry;
	    }
	}

	numFree += objectsPerBlob;
}


// Allocate entries until the free list is of size n (or until an allocation fails).
    void
Pool::preAllocate(
    unsigned		 n)
{
    _ASSERT(n);

    while (numFree < n) {
		unsigned oldNumFree = numFree;
		allocateMoreObjects();
		if (oldNumFree == numFree) {
		    // We can't allocate more; punt

		    return;
		}
    }
}

    PoolEntry *
Pool::getEntry(void)
{
	_ASSERT(!overlapEntriesAndObjects);	// This doesn't make sense in this case

    PoolEntry *entry = NULL;
    if (freeEntries) {
	    entry = freeEntries;
	    freeEntries = entry->next;
	    _ASSERT(entry->object == NULL);
    } else {
	    // Allocate a new entry blob and fill it in.
	    PoolBlob *blob = new PoolBlob;

        if (blob) {
            PoolEntry *blobEntries = new PoolEntry[entriesPerBlob];
	        if (blobEntries) {
		        blob->data = (void *)blobEntries;
		        // Release all of the newly allocated entries except the first one, which we'll return.
		        for (unsigned i = 1; i < entriesPerBlob; i++) {
		            releaseEntry(&blobEntries[i]);
		        }
		        entry = &blobEntries[0];

		        // Stick the new blob on the entries blob list.
		        blob->next = entriesBlobHead;
		        entriesBlobHead = blob;
	        } else {
		        // Give up; we couldn't get memory
		        delete blob;
	        }
	    }
    }
    return(entry);
}

    void
Pool::releaseEntry(
    PoolEntry 		*entry)
{
	_ASSERT(!overlapEntriesAndObjects);	// This doesn't make sense in this case.

    _ASSERT(entry);
    entry->object = NULL;
    entry->next = freeEntries;
    freeEntries = entry;
}

    void *
Pool::allocate(void)
{
    _ASSERT((numFree == 0) == (entries == NULL));

    if (!entries) {
	    allocateMoreObjects();
    }

	struct PoolEntry *thisEntry = entries;
	entries = entries->next;
	void *object = thisEntry->object;

    _ASSERT(object);

	if (overlapEntriesAndObjects) {
        _ASSERT(thisEntry->object == (PVOID)thisEntry);
	} else {
		releaseEntry(thisEntry);
	}

    _ASSERT(0 != numFree);
	numFree--;

	return(object);
}

    void
Pool::free(
    void		*object)
{
    _ASSERT(object);

    // No way to assert that this is the right kind (size) of object...

    struct PoolEntry *entry;
	if (overlapEntriesAndObjects) {
		entry = (PoolEntry *)object;
	} else {
	    // Get a PoolEntry.
	    entry = getEntry();
	    if (!entry) {
			// We couldn't get an entry, so we can't add this object to the free list.  Leak it.
			return;
	    }
	}

	numFree++;

    entry->object = object;
    entry->next = entries;
    entries = entry;
}
