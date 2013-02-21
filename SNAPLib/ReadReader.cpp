/*++

Module Name:

    ReadReader.cpp

Abstract:

    Concrete file reading classes

Environment:

    User mode service.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "Read.h"
#include "Aligner.h"
#include "FileFormat.h"

class SimpleReadReader : public ReadReader
{
public:
    SimpleReadReader(const FileFormat* i_format, DataReader* i_data, ReadClippingType i_clipping) 
        : format(i_format), data(i_data), clipping(i_clipping)
    {}

    virtual ~SimpleReadReader()
    {
        delete data;
    }

    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

    virtual bool getNextRead(Read *readToUpdate);
    
    virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, unsigned *genomeLocation, bool *isRC, unsigned *mapQ,
                    unsigned *flag, const char **cigar)
    {
        // return getNextRead(read,alignmentResult,genomeLocation,isRC,mapQ,flag,false,cigar);
    }

    void releaseBefore(DataBatch batch)
    { data->releaseBefore(batch); }


private:
    const FileFormat* format;
    DataReader* data;
    _int64 headerSize;
    const ReadClippingType clipping;
    const Genome* genome;
};
