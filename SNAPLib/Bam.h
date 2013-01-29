/*++

Module Name:

    Bam.h

Abstract:

    Binary Alignment Map (BAM) file writer.

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

--*/

#pragma once

#include "Compat.h"
#include "LandauVishkin.h"
#include "PairedEndAligner.h"
#include "VariableSizeVector.h"
#include "BufferedAsync.h"
#include "Read.h"

// BAM format layout
// SAM Format Specification v1.4-r985

// header information in each BGZF compression block
#pragma pack(push, 1)
struct BAMHeaderRefSeq;
struct BAMHeader
{
    static const _uint32 BAM_MAGIC = 0x014d4142; // 'BAM\1'
    _uint32     magic;
    
    _int32      l_text;
    
    char*       text() // not necessarily null terminated
    { return sizeof(magic) + sizeof(l_text) + (char*) this; }
    
    _int32      n_ref()
    { return * (_int32*) (l_text + text()); }

    BAMHeaderRefSeq* firstRefSeq()
    { return (BAMHeaderRefSeq*) (size(l_text) + (char*) this); }

    BAMHeader()
        : magic(BAM_MAGIC), l_text(0)
    {}

    // bytes required for a given text length, not including reference sequence blocks
    static size_t size(_int32 ltext)
    { return sizeof(_uint32) + 2 * sizeof(_int32) + ltext; }
};

// header information for each reference sequence record
struct BAMHeaderRefSeq
{
    _int32      l_name;
    
    char*       name()
    { return 4 + (char*) this; }

    _int32      l_ref()
    { return * (_int32*) (l_name + name()); }

    BAMHeaderRefSeq*    next()
    { return (BAMHeaderRefSeq*) (size(l_name) + (char*) this); }

    static size_t size(_int32 l_name)
    { return sizeof(l_name) + sizeof(_int32) + l_name; }
};

// information for each alignment record
struct BAMAlignAux;
struct BAMAlignment
{
    _int32      block_size;
    _int32      refID;
    _int32      pos;
    _uint8      l_read_name;
    _uint8      MAPQ;
    _uint16     bin;
    _uint16     n_cigar_op;
    _uint16     FLAG;
    _int32      l_seq;
    _int32      next_pos;
    _int32      tlen;
    
    char*       read_name()
    { return sizeof(tlen) + (char*) &this->tlen; }
    
    _uint32*    cigar()
    { return (_uint32*) (l_read_name + read_name()); }
    
    _uint8*     seq()
    { return (_uint8*) ((l_seq + 1) / 2 + cigar()); }

    char*       qual()
    { return (char*) (l_seq + seq()); }

    BAMAlignAux*    firstAux()
    { return (BAMAlignAux*) (l_seq + qual()); }
};

#define INT8_VAL_TYPE       'c'
#define UINT8_VAL_TYPE      'C'
#define INT16_VAL_TYPE      's'
#define UINT16_VAL_TYPE     'S'
#define INT32_VAL_TYPE      'i'
#define UINT32_VAL_TYPE     'I'
#define FLOAT_VAL_TYPE      'f'
#define CHAR_VAL_TYPE       'A'
#define STRING_VAL_TYPE     'Z'
#define HEX_VAL_TYPE        'H'
#define ARRAY_VAL_TYPE      'B'

// header for each auxiliary data field
struct BAMAlignAux
{
    char        tag[2];
    char        val_type;

    // accessors for single-valued fields

    void*       value()
    { return 3 + (char*) this; }
    
    _int8       int8Value()
    { _ASSERT(val_type == INT8_VAL_TYPE); return * (_int8*) value(); }

    _uint8      uint8Value()
    { _ASSERT(val_type == UINT8_VAL_TYPE); return * (_uint8*) value(); }
    
    _int16      int16Value()
    { _ASSERT(val_type == INT16_VAL_TYPE); return * (_int16*) value(); }

    _uint16     uint16Value()
    { _ASSERT(val_type == UINT16_VAL_TYPE); return * (_uint16*) value(); }
    
    _int32      int32Value()
    { _ASSERT(val_type == INT32_VAL_TYPE); return * (_int32*) value(); }

    _uint32     uint32Value()
    { _ASSERT(val_type == UINT32_VAL_TYPE); return * (_uint32*) value(); }
    
    float       floatValue()
    { _ASSERT(val_type == FLOAT_VAL_TYPE); return * (float*) value(); }
    
    // accessors for array fields

    _int32      count()
    { _ASSERT(val_type == ARRAY_VAL_TYPE); return * (_uint32*) (1 + (char*) value()); }

    char        arrayValType()
    { _ASSERT(val_type == ARRAY_VAL_TYPE); return * (char*) value(); }

    void*       data()
    { _ASSERT(val_type == ARRAY_VAL_TYPE); return 5 + (char*) value(); }
    
    _int8*      int8Array()
    { _ASSERT(arrayValType() == INT8_VAL_TYPE); return (_int8*) data(); }

    _uint8*     uint8Array()
    { _ASSERT(arrayValType() == UINT8_VAL_TYPE); return (_uint8*) data(); }
    
    _int16*     int16Array()
    { _ASSERT(arrayValType() == INT16_VAL_TYPE); return (_int16*) data(); }

    _uint16*    uint16Array()
    { _ASSERT(arrayValType() == UINT16_VAL_TYPE); return (_uint16*) data(); }
    
    _int32*     int32Array()
    { _ASSERT(arrayValType() == INT32_VAL_TYPE); return (_int32*) data(); }

    _uint32*    uint32Array()
    { _ASSERT(arrayValType() == UINT32_VAL_TYPE); return (_uint32*) data(); }
    
    float*      floatArray()
    { _ASSERT(arrayValType() == FLOAT_VAL_TYPE); return (float*) data(); }
    
    // compute overall size

    size_t      size()
    { return val_type != ARRAY_VAL_TYPE ? size(val_type) : arraySize(arrayValType(), count()); }
    
    static size_t size(char val_type)
    { return 3 + valueSize(val_type); }

    static size_t arraySize(char array_val_type, _uint32 count)
    { return 8 + valueSize(array_val_type) * (size_t) count; }

    static size_t valueSize(char val_type)
    {
        switch (val_type) {

        case INT8_VAL_TYPE:
        case UINT8_VAL_TYPE:
        case CHAR_VAL_TYPE:
            return 1;

        case INT16_VAL_TYPE:
        case UINT16_VAL_TYPE:
            return 2;

        case INT32_VAL_TYPE:
        case UINT32_VAL_TYPE:
        case FLOAT_VAL_TYPE:
            return 4;

        case ARRAY_VAL_TYPE:
        case STRING_VAL_TYPE:
        case HEX_VAL_TYPE:
        default:
            // todo: remove? log?
            _ASSERT(false);
            return 1;
        }
    }
};
#pragma pack(pop)


class BAMReader : public PairedReadReader, public ReadReader {
public:
        virtual ~BAMReader();

        virtual bool getNextRead(Read *readToUpdate);
    
        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, unsigned *genomeLocation, bool *isRC, unsigned *mapQ,
                        unsigned *flag, const char **cigar)
        {
            return getNextRead(read,alignmentResult,genomeLocation,isRC,mapQ,flag,false,cigar);
        }

        //
        // In getNextReadPair mapQ points to an array of two unsigneds.
        //
        virtual bool getNextReadPair(Read *read1, Read *read2, PairedAlignmentResult *alignmentResult, unsigned *mapQ, const char **cigar);

        //
        // The PairedReadReader version of getNextReadPair, which throws away the alignment, mapQ and cigar values.
        //
            bool
        getNextReadPair(Read *read1, Read *read2) {
            PairedAlignmentResult pairedAlignmentResult;
            unsigned mapQ[2];

            return getNextReadPair(read1,read2,&pairedAlignmentResult,mapQ,NULL);
        }

        virtual void readDoneWithBuffer(unsigned *referenceCount);

        static BAMReader* create(const char *fileName, const Genome *genome, _int64 startingOffset, _int64 amountOfFileToProcess, 
                                 ReadClippingType clipping = ClipBack);

        virtual ReadReader *getReaderToInitializeRead(int whichHalfOfPair) {
            _ASSERT(0 == whichHalfOfPair || 1 == whichHalfOfPair);
            return this;
        }

        static ReadSupplierGenerator *createReadSupplierGenerator(const char *fileName, int numThreads, const Genome *genome, ReadClippingType clipping = ClipBack);
        static PairedReadSupplierGenerator *createPairedReadSupplierGenerator(const char *fileName, int numThreads, const Genome *genome, ReadClippingType clipping = ClipBack);
protected:

        // result and fieldLengths must be of size nSAMFields
        static bool parseHeader(const char *fileName, char *firstLine, char *endOfBuffer, const Genome *genome, size_t *headerSize);
        
        static bool parseLine(char *line, char *endOfBuffer, char *result[], size_t *lineLength, size_t fieldLengths[]);

        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, 
                        unsigned *genomeLocation, bool *isRC, unsigned *mapQ, unsigned *flag, bool ignoreEndOfRange, const char **cigar);

        static void getReadFromLine(const Genome *genome, char *line, char *endOfBuffer, Read *read, AlignmentResult *alignmentResult,
                        unsigned *genomeLocation, bool *isRC, unsigned *mapQ, 
                        size_t *lineLength, unsigned *flag, const char **cigar, ReadClippingType clipping);
};
