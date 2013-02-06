/*++

Module Name:

    MultiInputReadSupplier.cp

Abstract:

    A read supplier that combines other read suppliers.  It's used when there are muliple input files to process.

Authors:

    Bill Bolosky, November, 2012

Environment:

    User mode service.

Revision History:


--*/

#pragma once
#include "stdafx.h"
#include "Read.h"
#include "Compat.h"
#include "MultiInputReadSupplier.h"


MultiInputReadSupplier::MultiInputReadSupplier(int i_nReadSuppliers, ReadSupplier **i_readSuppliers)
{
    readSuppliers = i_readSuppliers;    // We get to own the array
    nRemainingReadSuppliers = nReadSuppliers = i_nReadSuppliers;
    nextReadSupplier = 0;
    activeReadSuppliers = new ActiveRead[nReadSuppliers];
    for (int i = 0; i < nReadSuppliers; i++) {
        activeReadSuppliers[i].index = i;
        activeReadSuppliers[i].lastBatch = DataBatch();
        activeReadSuppliers[i].firstReadInNextBatch = NULL;
    }
}

MultiInputReadSupplier::~MultiInputReadSupplier()
{
    for (int i = 0; i < nReadSuppliers; i++) {
        delete readSuppliers[i];
        readSuppliers[i] = NULL;
    }
    delete [] readSuppliers;
    delete [] activeReadSuppliers;
}

    Read *
MultiInputReadSupplier::getNextRead()
{
    if (0 == nRemainingReadSuppliers) {
        return NULL;
    }

    _ASSERT(nextReadSupplier < nRemainingReadSuppliers);

    Read *read = readSuppliers[activeReadSuppliers[nextReadSupplier].index]->getNextRead();
    bool nextBatch = read != NULL && read->getBatch() != activeReadSuppliers[nextReadSupplier].lastBatch;
    if (nextBatch) {
        activeReadSuppliers[nextReadSupplier].firstReadInNextBatch = read;
        activeReadSuppliers[nextReadSupplier].lastBatch = read->getBatch();
    }
    if (read == NULL || nextBatch) {
        // end of batch from current supplier, round-robin through suppliers
        while (true) {
            if (nextBatch) {
                nextReadSupplier = (nextReadSupplier + 1) % nRemainingReadSuppliers;
                nextBatch = false;
            }
            if (activeReadSuppliers[nextReadSupplier].firstReadInNextBatch != NULL) {
                read = activeReadSuppliers[nextReadSupplier].firstReadInNextBatch;
                activeReadSuppliers[nextReadSupplier].firstReadInNextBatch = NULL;
                break;
            }
            read = readSuppliers[activeReadSuppliers[nextReadSupplier].index]->getNextRead();
            if (read != NULL) {
                break;
            }
            nRemainingReadSuppliers--;
            if (0 == nRemainingReadSuppliers) {
                return NULL;
            }
            //
            // This supplier is done.  Update our array to pull the
            // last live read supplier into the slot that we just vacated (this will result
            // in violating a strict round robin, but we don't promise any such thing
            // anyway). Can't delete because it might be retaining read data in use downstream.
            //
            activeReadSuppliers[nextReadSupplier] = activeReadSuppliers[nRemainingReadSuppliers];
        }
    }
    // update batch info
    activeReadSuppliers[nextReadSupplier].lastBatch = read->getBatch();
    read->setBatch(DataBatch(1 + activeReadSuppliers[nextReadSupplier].index, read->getBatch().batchID));

    return read;
}

    void
MultiInputReadSupplier::releaseBefore(
    DataBatch batch)
{
    int index = batch.fileID - 1;
    _ASSERT(index >= 0 && index < nReadSuppliers);
    readSuppliers[index]->releaseBefore(DataBatch(batch.batchID));
}

MultiInputPairedReadSupplier::MultiInputPairedReadSupplier(int i_nReadSuppliers, PairedReadSupplier **i_pairedReadSuppliers)
{
    pairedReadSuppliers = i_pairedReadSuppliers;    // We get to own the array
    nRemainingReadSuppliers = nReadSuppliers = i_nReadSuppliers;
    nextReadSupplier = 0;
    activeReadSuppliers = new ActiveRead[nReadSuppliers];
    for (int i = 0; i < nReadSuppliers; i++) {
        activeReadSuppliers[i].index = i;
        activeReadSuppliers[i].lastBatch[0] = activeReadSuppliers[i].lastBatch[1] = DataBatch();
        activeReadSuppliers[i].firstReadInNextBatch[0] = activeReadSuppliers[i].firstReadInNextBatch[1] = NULL;
    }
}
 
MultiInputPairedReadSupplier::~MultiInputPairedReadSupplier()
{
    for (int i = 0; i < nReadSuppliers; i++) {
        delete pairedReadSuppliers[i];
        pairedReadSuppliers[i] = NULL;
    }

    delete [] pairedReadSuppliers;
    delete [] activeReadSuppliers;
}

    bool 
MultiInputPairedReadSupplier::getNextReadPair(Read **read0, Read **read1)
{
    if (0 == nRemainingReadSuppliers) {
        return NULL;
    }

    _ASSERT(nextReadSupplier < nRemainingReadSuppliers);

    ActiveRead* active = &activeReadSuppliers[nextReadSupplier];
    bool hasReads = pairedReadSuppliers[active->index]->getNextReadPair(read0, read1);
    bool nextBatch = hasReads && ((*read0)->getBatch() != active->lastBatch[0] || (*read1)->getBatch() != active->lastBatch[1]);
    if (nextBatch) {
        active->firstReadInNextBatch[0] = *read0;
        active->firstReadInNextBatch[1] = *read1;
    }
    if (nextBatch || ! hasReads) {
        while (true) {
            // end of batch from current supplier, round-robin through suppliers
            if (nextBatch) {
                nextReadSupplier = (nextReadSupplier + 1) % nRemainingReadSuppliers;
                nextBatch = false;
            }
            active = &activeReadSuppliers[nextReadSupplier];
            if (active->firstReadInNextBatch[0] != NULL) {
                _ASSERT(active->firstReadInNextBatch[1] != NULL);
                *read0 = active->firstReadInNextBatch[0];
                *read1 = active->firstReadInNextBatch[1];
                active->firstReadInNextBatch[0] = active->firstReadInNextBatch[1] = NULL;
                break;
            }
            if (pairedReadSuppliers[active->index]->getNextReadPair(read0, read1)) {
                break;
            }
            nRemainingReadSuppliers--;
            if (0 == nRemainingReadSuppliers) {
                return false;
            }
            //
            // This supplier is done.  Update our array to pull the
            // last live read supplier into the slot that we just vacated (this will result
            // in violating a strict round robin, but we don't promise any such thing
            // anyway). Can't delete because it might be retaining read data in use downstream.
            //
            activeReadSuppliers[nextReadSupplier] = activeReadSuppliers[nRemainingReadSuppliers];
        }
    }
    active->lastBatch[0] = (*read0)->getBatch();
    (*read0)->setBatch(DataBatch(active->index + 1, active->lastBatch[0].batchID));
    active->lastBatch[1] = (*read1)->getBatch();
    (*read1)->setBatch(DataBatch(active->index + 1, active->lastBatch[1].batchID));

    return true;
}

    void
MultiInputPairedReadSupplier::releaseBefore(
    DataBatch batch)
{
    int index = batch.fileID - 1;
    _ASSERT(index >= 0 && index < nReadSuppliers);
    pairedReadSuppliers[index]->releaseBefore(DataBatch(batch.batchID));
}


MultiInputReadSupplierGenerator::MultiInputReadSupplierGenerator(int i_nReadSuppliers, ReadSupplierGenerator **i_readSupplierGenerators)
{
    nReadSuppliers = i_nReadSuppliers;
    readSupplierGenerators = i_readSupplierGenerators;  // We take ownership of the array
}

MultiInputReadSupplierGenerator::~MultiInputReadSupplierGenerator()
{
    for (int i = 0; i < nReadSuppliers; i++) {
        delete readSupplierGenerators[i];
    }
    delete [] readSupplierGenerators;
    readSupplierGenerators = NULL;
}

    ReadSupplier *
MultiInputReadSupplierGenerator::generateNewReadSupplier()
{
    ReadSupplier **readSuppliers = new ReadSupplier *[nReadSuppliers];
    for (int i = 0; i < nReadSuppliers; ) {
        readSuppliers[i] = readSupplierGenerators[i]->generateNewReadSupplier();
		if (NULL == readSuppliers[i]) {
			nReadSuppliers--;
		} else {
			i++;
		}
    }

	if (0 == nReadSuppliers) {
		delete [] readSuppliers;
		return NULL;
	}

    return new MultiInputReadSupplier(nReadSuppliers,readSuppliers);    // The Supplier owns the array and suppliers we created
}



MultiInputPairedReadSupplierGenerator::MultiInputPairedReadSupplierGenerator(int i_nReadSuppliers, PairedReadSupplierGenerator **i_readSupplierGenerators)
{
    nReadSuppliers = i_nReadSuppliers;
    readSupplierGenerators = i_readSupplierGenerators;  // We own the array and the generators.

}

MultiInputPairedReadSupplierGenerator::~MultiInputPairedReadSupplierGenerator()
{
    for (int i = 0; i < nReadSuppliers; i++) {
        delete readSupplierGenerators[i];
        readSupplierGenerators[i] = NULL;
    }
    delete [] readSupplierGenerators;
    readSupplierGenerators = NULL;
}

    PairedReadSupplier *
MultiInputPairedReadSupplierGenerator::generateNewPairedReadSupplier()
{
    PairedReadSupplier **readSuppliers = new PairedReadSupplier *[nReadSuppliers];
    for (int i = 0; i < nReadSuppliers; ) {
        readSuppliers[i] = readSupplierGenerators[i]->generateNewPairedReadSupplier();
		if (NULL == readSuppliers[i]) {
			nReadSuppliers--;
		} else {
			i++;
		}
    }
	if (0 == nReadSuppliers) {
		delete [] readSuppliers;
		return NULL;
	}
    return new MultiInputPairedReadSupplier(nReadSuppliers,readSuppliers);
}
