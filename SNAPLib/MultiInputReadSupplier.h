/*++

Module Name:

    MultiInputReadSupplier.h

Abstract:

    Headers for a read supplier that combines other read suppliers.  It's used when there are muliple input files to process.

Authors:

    Bill Bolosky, November, 2012

Environment:

    User mode service.

Revision History:


--*/

#pragma once
#include "Read.h"
#include "Compat.h"

class MultiInputReadSupplier: public ReadSupplier {
public:
    MultiInputReadSupplier(int nReadSuppliers, ReadSupplier **i_readSuppliers);
    virtual ~MultiInputReadSupplier();

    virtual Read *getNextRead();


private:
    int                 nRemainingReadSuppliers;
    int                 nextReadSupplier;
    ReadSupplier        **readSuppliers;
};

class MultiInputPairedReadSupplier: public PairedReadSupplier {
public:
    MultiInputPairedReadSupplier(int nReadSuppliers, PairedReadSupplier **i_pairedReadSuppliers);
    virtual ~MultiInputPairedReadSupplier();

    virtual bool getNextReadPair(Read **read0, Read **read1);

private:

    int                 nRemainingReadSuppliers;
    int                 nextReadSupplier;
    PairedReadSupplier  **pairedReadSuppliers;
};

class MultiInputReadSupplierGenerator: public ReadSupplierGenerator
{
public:
    MultiInputReadSupplierGenerator(int i_nReadSuppliers, ReadSupplierGenerator **i_readSupplierGenerators);
    virtual ~MultiInputReadSupplierGenerator();

    virtual ReadSupplier *generateNewReadSupplier();

private:

    int nReadSuppliers;
    ReadSupplierGenerator **readSupplierGenerators;
};

class MultiInputPairedReadSupplierGenerator: public PairedReadSupplierGenerator
{
public:
    MultiInputPairedReadSupplierGenerator(int i_nReadSuppliers, PairedReadSupplierGenerator **i_readSupplierGenerators);
    virtual ~MultiInputPairedReadSupplierGenerator();

    virtual PairedReadSupplier *generateNewPairedReadSupplier();

private:

    int nReadSuppliers;
    PairedReadSupplierGenerator **readSupplierGenerators;
};
