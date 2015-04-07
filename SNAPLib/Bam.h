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
#include "SAM.h"
#include "Read.h"
#include "DataReader.h"

// for debugging file I/O, validate BAM records on input & output
//#define VALIDATE_BAM

// BAM format layout
// SAM Format Specification v1.4-r985

// header information in each BGZF compression block

#define BAM_BLOCK 65536

#pragma pack(push, 1)
struct BAMHeaderRefSeq;
struct BAMHeader
{
    static const _uint32 BAM_MAGIC = 0x014d4142; // 'BAM\1'
    _uint32     magic;
    
    _int32      l_text;
    
    char*       text() // not necessarily null terminated
    { return sizeof(magic) + sizeof(l_text) + (char*) this; }
    
    _int32&     n_ref()
    { return * (_int32*) (l_text + text()); }

    BAMHeaderRefSeq* firstRefSeq()
    { return (BAMHeaderRefSeq*) (size(l_text) + (char*) this); }

    BAMHeader()
        : magic(BAM_MAGIC), l_text(0)
    {}

    size_t size()
    { return size(l_text); }

    // bytes required for a given text length, not including reference sequence blocks
    static size_t size(_int32 ltext)
    { return sizeof(BAMHeader) + sizeof(_int32)/*n_ref*/ + ltext; }
};

// header information for each reference sequence record
struct BAMHeaderRefSeq
{
    _int32      l_name;
    
    char*       name()
    { return 4 + (char*) this; }

    _int32&     l_ref()
    { return * (_int32*) (l_name + name()); }

    BAMHeaderRefSeq*    next()
    { return (BAMHeaderRefSeq*) (size(l_name) + (char*) this); }

    static size_t size(_int32 l_name)
    { return sizeof(l_name) + sizeof(_int32)/*l_ref*/ + l_name; }
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
    _int32      next_refID;
    _int32      next_pos;
    _int32      tlen;
    
    char*       read_name()
    { return sizeof(tlen) + (char*) &this->tlen; }
    
    _uint32*    cigar()
    { return (_uint32*) (l_read_name + read_name()); }
    
    _uint8*     seq()
    { return (_uint8*) (n_cigar_op + cigar()); }

    char*       qual()
    { return (char*) ((l_seq + 1) / 2 + seq()); }

    BAMAlignAux*    firstAux()
    { return (BAMAlignAux*) (l_seq + qual()); }

    unsigned    auxLen()
    { return (unsigned) (size() - ((char*) firstAux() - (char*) this)); }

    BAMAlignAux*    endAux()
    { return (BAMAlignAux*) (auxLen() + (char*) firstAux()); }

    size_t size()
    { return block_size + sizeof(block_size); }

    static size_t size(unsigned l_read_name, unsigned n_cigar_op, unsigned l_seq, unsigned l_aux)
    { return sizeof(BAMAlignment) + l_read_name + n_cigar_op * sizeof(_uint32) + (l_seq + 1) / 2 + l_seq + l_aux; }

    // conversions

    static const char* CodeToSeq;
    static const char* CodeToSeqRC;
    static _uint16 CodeToSeqPair[256];
    static _uint16 CodeToSeqPairRC[256];
    static _uint8 SeqToCode[256];
    static const char* CodeToCigar;
    static _uint8 CigarToCode[256];
    static _uint8 CigarCodeToRefBase[9];
    static int GetCigarOpCode(_uint32 op) { return op & 0xf; }
    static int GetCigarOpCount(_uint32 op) { return op >> 4; }
    
    static void decodeSeq(char* o_sequence, const _uint8* nibbles, int bases);
    static void decodeQual(char* o_qual, char* quality, int bases);
    static void decodeSeqRC(char* o_sequence, const _uint8* nibbles, int bases);
    static void decodeQualRC(char* o_qual, char* quality, int bases);
    static bool decodeCigar(char* o_cigar, int cigarSize, _uint32* cigar, int ops);
    static void getClippingFromCigar(_uint32 *cigar, int ops, unsigned *o_frontClipping, unsigned *o_backClipping, unsigned *o_frontHardClipping, unsigned *o_backHardClipping);

    static void encodeSeq(_uint8* nibbles, char* ascii, int length);

    int l_ref(); // length of reference aligned to read

    class _init { public: _init(); };
    static _init _init_;

    // binning

    static const _uint32 BAM_EXTRA_BIN = 37450; // extra bin for metadata

    /* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
    static int reg2bin(int beg, int end);
    /* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
    static const int MAX_BIN = (((1<<18)-1)/7);
    static int reg2bins(int beg, int end, _uint16* list/*[MAX_BIN]*/);

    // absoluate genome locations

    GenomeLocation getLocation(const Genome* genome) const
    { return genome == NULL || pos < 0 || refID < 0 || refID >= genome->getNumContigs() || (FLAG & SAM_UNMAPPED)
        ? UINT32_MAX : (genome->getContigs()[refID].beginningLocation + pos); }

    GenomeLocation getNextLocation(const Genome* genome) const
    { return next_pos < 0 || next_refID < 0 || (FLAG & SAM_NEXT_UNMAPPED) ? UINT32_MAX : (genome->getContigs()[next_refID].beginningLocation + next_pos); }

#ifdef VALIDATE_BAM
    void validate();
#else
    inline void validate() {} // inline noop
#endif
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
    {
        return val_type == STRING_VAL_TYPE ? strlen((const char*) value()) + 4
            : val_type == ARRAY_VAL_TYPE ? size(arrayValType(), count())
            : size(val_type);
    }
    
    static size_t size(char val_type)
    { return 3 + valueSize(val_type); }

    static size_t size(char array_val_type, _uint32 count)
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

    BAMAlignAux* next()
    { return (BAMAlignAux*) (size() + (char*) this); }
};

struct BgzfExtra
{
    _uint8 SI1;
    _uint8 SI2;
    _uint16 SLEN;
    
    void* data()
    { return this + 1; }
    
    BgzfExtra* nextExtra()
    { return (BgzfExtra*) (SLEN + (char*) data()); }
};

struct BgzfHeader
{
    _uint8 ID1;		// 0x1f
    _uint8 ID2;		// 0x8b (magic numbers)
    _uint8 CM;		// Compression method (8 == deflate)
    _uint8 FLG;		// flags
    _uint32 MTIME;	// Mod time
    _uint8 XFL;		// extra flags
    _uint8 OS;		// operating system 
    _uint16 XLEN;	// extra length

    BgzfExtra* firstExtra()
    { return (BgzfExtra*) (sizeof(_uint16) + (char*) &this->XLEN); }

    _uint16 BSIZE()
    {
        for (BgzfExtra* x = firstExtra(); (char*) x < XLEN + (char*) firstExtra(); x = x->nextExtra()) {
            if (x->SI1 == 66 && x->SI2 == 67) {
                _ASSERT(x->SLEN == 2);
                return * (_uint16*) x->data();
            }
        }
        _ASSERT(false);
        return 0;
    }

    _uint32 ISIZE()
    {
        _uint32 result = * (_uint32*) (BSIZE() - 3 + (char*) this);
        _ASSERT(result <= BAM_BLOCK);
        return result;
    }

    bool validate(size_t compressed, size_t uncompressed);

    static bool validate(char* buffer, size_t bytes);
};


#pragma pack(pop)


class BAMReader : public PairedReadReader, public ReadReader {
public:

        BAMReader(const ReaderContext& i_context);

        virtual ~BAMReader();

        void init(const char *fileName, int bufferCount, _int64 startingOffset, _int64 amountOfFileToProcess);

        virtual bool getNextRead(Read *readToUpdate)
        {
            return getNextRead(readToUpdate, NULL, NULL, NULL, NULL, NULL, false, NULL);
        }
    
        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, GenomeLocation *genomeLocation, bool *isRC, unsigned *mapQ,
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
        bool getNextReadPair(Read *read1, Read *read2)
        {
            return getNextReadPair(read1,read2,NULL,NULL,NULL);
        }
        
        void holdBatch(DataBatch batch)
        { data->holdBatch(batch); }

        bool releaseBatch(DataBatch batch)
        { return data->releaseBatch(batch); }

        virtual ReaderContext* getContext()
        { return ((ReadReader*)this)->getContext(); }

        static BAMReader* create(const char *fileName, int bufferCount,
            _int64 startingOffset, _int64 amountOfFileToProcess, 
            const ReaderContext& context);
        
        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);
        
        static ReadSupplierGenerator *createReadSupplierGenerator(const char *fileName, int numThreads, const ReaderContext& context);
        
        static PairedReadSupplierGenerator *createPairedReadSupplierGenerator(const char *fileName, int numThreads, bool quicklyDropUnmatchedReads, 
            const ReaderContext& context, int matchBufferSize = 5000);

        static const int MAX_SEQ_LENGTH;
        static const int MAX_RECORD_LENGTH;

protected:

        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, 
                        GenomeLocation *genomeLocation, bool *isRC, unsigned *mapQ, unsigned *flag, bool ignoreEndOfRange, const char **cigar);

        void getReadFromLine(const Genome *genome, char *line, char *endOfBuffer, Read *read, AlignmentResult *alignmentResult,
                        GenomeLocation *genomeLocation, bool *isRC, unsigned *mapQ, 
                        size_t *lineLength, unsigned *flag, const char **cigar, ReadClippingType clipping);

private:
        void readHeader(const char* fileName);


        char* getExtra(_int64 bytes);

        DataReader*         data;
        //unsigned            n_ref; // number of reference sequences
        //unsigned*           refOffset; // array mapping ref sequence ID to contig location
        _int64              extraOffset; // offset into extra data
};
