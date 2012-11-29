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


MultiInputReadSupplier::MultiInputReadSupplier(int nReadSuppliers, ReadSupplier **i_readSuppliers)
{
    readSuppliers = i_readSuppliers;    // We get to own the array
    nRemainingReadSuppliers = nReadSuppliers;
    nextReadSupplier = 0;
}
MultiInputReadSupplier::~MultiInputReadSupplier()
{
    for (int i = 0; i < nRemainingReadSuppliers; i++) {
        delete readSuppliers[i];
        readSuppliers[i] = NULL;
    }
    delete [] readSuppliers;
}

    Read *
MultiInputReadSupplier::getNextRead()
{
    if (0 == nRemainingReadSuppliers) {
        return NULL;
    }

    _ASSERT(nextReadSupplier < nRemainingReadSuppliers);

    Read *read;
    while (NULL == (read = readSuppliers[nextReadSupplier]->getNextRead())) {
        nRemainingReadSuppliers--;
        if (0 == nRemainingReadSuppliers) {
            return NULL;
        }
        //
        // This supplier is done.  Delete it, and then update our array to pull the
        // last live read supplier into the slot that we just vacated (this will result
        // in violating a strict round robin, but we don't promise any such thing
        // anyway).
        //
        delete readSuppliers[nextReadSupplier];
        readSuppliers[nextReadSupplier] = readSuppliers[nRemainingReadSuppliers];
        readSuppliers[nRemainingReadSuppliers] = NULL;
        nextReadSupplier = 0;
    }

    nextReadSupplier = (nextReadSupplier + 1) % nRemainingReadSuppliers;
    return read;
}



MultiInputPairedReadSupplier::MultiInputPairedReadSupplier(int nReadSuppliers, PairedReadSupplier **i_pairedReadSuppliers)
{
    pairedReadSuppliers = i_pairedReadSuppliers;    // We get to own the array
    nRemainingReadSuppliers = nReadSuppliers;
    nextReadSupplier = 0;
}
 
MultiInputPairedReadSupplier::~MultiInputPairedReadSupplier()
{
    for (int i = 0; i < nRemainingReadSuppliers; i++) {
        delete pairedReadSuppliers[i];
        pairedReadSuppliers[i] = NULL;
    }

    delete [] pairedReadSuppliers;
}

    bool 
MultiInputPairedReadSupplier::getNextReadPair(Read **read0, Read **read1)
{
    if (0 == nRemainingReadSuppliers) {
        return NULL;
    }

    _ASSERT(nextReadSupplier < nRemainingReadSuppliers);

    while (!pairedReadSuppliers[nextReadSupplier]->getNextReadPair(read0, read1)) {
        nRemainingReadSuppliers--;
        if (0 == nRemainingReadSuppliers) {
            return false;
        }
        //
        // This supplier is done.  Delete it, and then update our array to pull the
        // last live read supplier into the slot that we just vacated (this will result
        // in violating a strict round robin, but we don't promise any such thing
        // anyway).
        //
        delete pairedReadSuppliers[nextReadSupplier];
        pairedReadSuppliers[nextReadSupplier] = pairedReadSuppliers[nRemainingReadSuppliers];
        pairedReadSuppliers[nRemainingReadSuppliers] = NULL;
        nextReadSupplier = 0;
    }

    nextReadSupplier = (nextReadSupplier + 1) % nRemainingReadSuppliers;
    return true;
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
    for (int i = 0; i < nReadSuppliers; i++) {
        readSuppliers[i] = readSupplierGenerators[i]->generateNewReadSupplier();
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
    for (int i = 0; i < nReadSuppliers; i++) {
        readSuppliers[i] = readSupplierGenerators[i]->generateNewPairedReadSupplier();
    }
    return new MultiInputPairedReadSupplier(nReadSuppliers,readSuppliers);
}
