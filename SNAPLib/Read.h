/*++


Module Name:

    Read.h

Abstract:

    Headers for the Read class for the SNAP sequencer

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once

#include "Compat.h"
#include "Tables.h"
#include "DataReader.h"
#include "DataWriter.h"
#include "directions.h"

class FileFormat;

class Genome;
class GTFReader;

struct PairedAlignmentResult;

enum AlignmentResult {NotFound, SingleHit, MultipleHits, UnknownAlignment}; // BB: Changed Unknown to UnknownAlignment because of a conflict w/Windows headers

bool isAValidAlignmentResult(AlignmentResult result);
// constant for small/medium/large reads
#define MAX_READ_LENGTH 500
//#define MAX_READ_LENGTH 1000
///#define MAX_READ_LENGTH 20000

//
// Here's a brief description of the classes for input in SNAP:
// Read:
//      A Read is some data that's come from a NGS machine.  It includes some bases and associated quality score, as well as an ID.
//      Reads may be clipped (because the sequencing machine was unsure of some bases).  They may be switched between forward and
//      reverse complement sense.  They may or may not own their own memory for the various fields.
//
// ReadReader:
//      A ReadReader understands how to generate reads from some input source (i.e., a FASTQ, SAM, BAM or CRAM file, for instance).
//      It owns the storage for the read's information (i.e., the base string), but does not own the Read object itself.  It is responsible
//      for assuring that the memory for the read data is valid for the lifetime of the ReadReader (which, in practice, means it needs
//      to use mapped files).  ReadReaders may assume that they will only be called from one thread.
//
// PairedReadReader:
//      Similar to a ReadReader, except that it gets mate pairs of Reads.
//
// ReadSupplier:
//      A class that supplies reads to a consumer.  It looks similar to a ReadReader, except that it own the storage for the
//      ReadObject.  The idea here is to allow the supplier to manage the memory that the Read object lives in so that a supplier
//      can be implemented by a parallel queue with batches of reads in it.  Supplier may, of course, also be implemented in
//      different ways, such as range splitters.  Like ReadReaders, ReadSuppliers will be called from only one thread.  In
//      practice, ReadSuppliers will have underlying ReadReaders (which might be behind a shared queue, for example).
//
// ReadSupplierGenerator:
//      A class that creates a ReadSupplier.  This has to be thread safe.  The usual pattern is that the initialization code will
//      create a read supplier generator, which will then be called on each of the threads to create a supplier, which will supply
//      the reads to be aligned.
//
// PairedReadSupplierGenerator:
//      The paired version of a ReadSupplier.
//

const int MaxReadLength = MAX_READ_LENGTH;

class Read;

enum ReadClippingType {NoClipping, ClipFront, ClipBack, ClipFrontAndBack};

struct ReaderContext
{
    const Genome*       genome;
    const Genome*       transcriptome;
    GTFReader*          gtf;
    const char*         defaultReadGroup;
    ReadClippingType    clipping;
    bool                paired;
    bool                ignoreSecondaryAlignments;   // Should we just ignore reads with the Secondary Alignment bit set?
    const char*         header; // allocated buffer for header
    size_t              headerLength; // length of string
    size_t              headerBytes; // bytes used for header in file
    bool                headerMatchesIndex; // header refseq matches current index
};

class ReadReader {
public:
    ReadReader(const ReaderContext& i_context) : context(i_context) {}

    virtual ~ReadReader() {}
        
    // reading

    virtual bool getNextRead(Read *readToUpdate) = 0;
    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess) = 0;

    virtual void releaseBatch(DataBatch batch) = 0;

protected:
    ReaderContext context;
};

class PairedReadReader {
public:
    virtual ~PairedReadReader() {}

    // reading

    virtual bool getNextReadPair(Read *read1, Read *read2) = 0;
    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess) = 0;

    virtual void releaseBatch(DataBatch batch) = 0;

    // wrap a single read source with a matcher that buffers reads until their mate is found
    static PairedReadReader* PairMatcher(ReadReader* single, bool autoRelease);
};

class ReadSupplier {
public:
    virtual Read *getNextRead() = 0;    // This read is valid until you call getNextRead, then it's done.  Don't worry about deallocating it.
    virtual ~ReadSupplier() {}

    virtual void releaseBatch(DataBatch batch) = 0;
};

class PairedReadSupplier {
public:
    // These read are valid until you call getNextRead, then they're done.  Don't worry about deallocating them.
    virtual bool getNextReadPair(Read **read0, Read **read1) = 0;
    virtual ~PairedReadSupplier() {}

    virtual void releaseBatch(DataBatch batch) = 0;
};

class ReadSupplierGenerator {
public:
    virtual ReadSupplier *generateNewReadSupplier() = 0;
    virtual ~ReadSupplierGenerator() {}
};

class PairedReadSupplierGenerator {
public:
    virtual PairedReadSupplier *generateNewPairedReadSupplier() = 0;
    virtual ~PairedReadSupplierGenerator() {}
};

class ReadWriter {
public:

    virtual ~ReadWriter() {}

    // write out header
    virtual bool writeHeader(const ReaderContext& context, bool sorted, int argc, const char **argv, const char *version, const char *rgLine) = 0;

    // write a single read, return true if successful
    virtual bool writeRead(Read *read, AlignmentResult result, int mapQuality, unsigned genomeLocation, Direction direction, bool isTranscriptome, unsigned tlocation) = 0;

    // write a pair of reads, return true if successful
    virtual bool writePair(Read *read0, Read *read1, PairedAlignmentResult *result) = 0;

    // close out this thread
    virtual void close() = 0;
};

class DataWriterSupplier;

class ReadWriterSupplier
{
public:
    virtual ReadWriter* getWriter() = 0;

    virtual void close() = 0;

    static ReadWriterSupplier* create(const FileFormat* format, DataWriterSupplier* dataSupplier,
        const Genome* genome, const Genome* transcriptome, const GTFReader* gtf);
};

#define READ_GROUP_FROM_AUX     ((const char*) -1)
    
class Read {
public:
        Read() :    
            id(NULL), data(NULL), quality(NULL), 
            localUnclippedDataBuffer(NULL), localUnclippedQualityBuffer(NULL), localBufferSize(0),
            originalUnclippedDataBuffer(NULL), originalUnclippedQualityBuffer(NULL), clippingState(NoClipping),
            upperCaseDataBuffer(NULL), upperCaseDataBufferLength(0), auxiliaryData(NULL), auxiliaryDataLength(0),
            readGroup(NULL), originalAlignedLocation(-1), originalMAPQ(-1), originalSAMFlags(0),
            originalFrontClipping(0), originalBackClipping(0), originalFrontHardClipping(0), originalBackHardClipping(0)
        {}

        Read(const Read& other) : 
            id(other.id), data(other.data), quality(other.quality), 
            idLength(other.idLength), frontClippedLength(other.frontClippedLength), dataLength(other.dataLength),
            unclippedData(other.unclippedData), unclippedLength(other.unclippedLength), unclippedQuality(other.unclippedQuality),
            localUnclippedDataBuffer(other.localUnclippedDataBuffer),
            localUnclippedQualityBuffer(other.localUnclippedQualityBuffer),
            localBufferSize(other.localBufferSize),
            originalUnclippedDataBuffer(other.originalUnclippedDataBuffer),
            originalUnclippedQualityBuffer(other.originalUnclippedQualityBuffer),
            clippingState(other.clippingState),
            batch(other.batch), 
            upperCaseDataBuffer(other.upperCaseDataBuffer), upperCaseDataBufferLength(other.upperCaseDataBufferLength),
            auxiliaryData(other.auxiliaryData), auxiliaryDataLength(other.auxiliaryDataLength), readGroup(other.readGroup),
            originalAlignedLocation(other.originalAlignedLocation), originalMAPQ(other.originalMAPQ), 
            originalSAMFlags(other.originalSAMFlags), originalFrontClipping(other.originalFrontClipping), 
            originalBackClipping(other.originalBackClipping), originalFrontHardClipping(other.originalFrontHardClipping), 
            originalBackHardClipping(other.originalBackHardClipping)
        {
            Read* o = (Read*) &other; // hack!
            o->localUnclippedDataBuffer = NULL;
            o->localUnclippedQualityBuffer = NULL;
            o->localBufferSize = 0;
            o->originalUnclippedDataBuffer = NULL;
            o->originalUnclippedQualityBuffer = NULL;
            o->upperCaseDataBuffer = NULL;
            o->upperCaseDataBufferLength = 0;
        }

        ~Read()
        {
            delete [] localUnclippedDataBuffer;
            delete [] localUnclippedQualityBuffer;
            BigDealloc(upperCaseDataBuffer);
        }

        void dispose()
        {
            delete [] localUnclippedDataBuffer;
            localUnclippedDataBuffer = NULL;
            delete [] localUnclippedQualityBuffer;
            localUnclippedQualityBuffer = NULL;
            BigDealloc(upperCaseDataBuffer);
            upperCaseDataBuffer = NULL;
        }

        void operator=(Read& other)
        {
            id = other.id;
            data = other.data;
            quality = other.quality;
            idLength = other.idLength;
            frontClippedLength = other.frontClippedLength;
            dataLength = other.dataLength;
            unclippedData = other.unclippedData;
            unclippedLength = other.unclippedLength;
            unclippedQuality = other.unclippedQuality;
            delete [] localUnclippedDataBuffer;
            localUnclippedDataBuffer = other.localUnclippedDataBuffer;
            delete [] localUnclippedQualityBuffer;
            localUnclippedQualityBuffer = other.localUnclippedQualityBuffer;
            localBufferSize = other.localBufferSize;
            other.localUnclippedDataBuffer = NULL;
            other.localUnclippedQualityBuffer = NULL;
            other.localBufferSize = 0;
            originalUnclippedDataBuffer = other.originalUnclippedDataBuffer;
            originalUnclippedQualityBuffer = other.originalUnclippedQualityBuffer;
            other.originalUnclippedDataBuffer = NULL;
            other.originalUnclippedQualityBuffer = NULL;
            clippingState = other.clippingState;
            batch = other.batch;
            BigDealloc(upperCaseDataBuffer);
            upperCaseDataBuffer = other.upperCaseDataBuffer;
            upperCaseDataBufferLength = other.upperCaseDataBufferLength;
            other.upperCaseDataBuffer = NULL;
            other.upperCaseDataBufferLength = 0;
            readGroup = other.readGroup;
            auxiliaryData = other.auxiliaryData;
            auxiliaryDataLength = other.auxiliaryDataLength;
            originalAlignedLocation = other.originalAlignedLocation;
            originalMAPQ = other.originalMAPQ;
            originalSAMFlags = other.originalSAMFlags;
            originalFrontClipping = other.originalFrontClipping;
            originalBackClipping = other.originalBackClipping;            
            originalFrontHardClipping = other.originalFrontHardClipping;
            originalBackHardClipping = other.originalBackHardClipping;
        }

        //
        // Initialize the Read.  Reads do NOT take ownership of the memory to which they
        // point, and it's the caller's responsibility to make sure that it continues to
        // exist as long as the Read does.  This is so that the caller can read a bunch of
        // read data into a buffer, and then carve Reads out of it without doing further
        // memory allocations, which would slow down the sequencing.
        //

        void init(
                const char *i_id, 
                unsigned i_idLength,
                const char *i_data, 
                const char *i_quality, 
                unsigned i_dataLength)
        {
            init(i_id, i_idLength, i_data, i_quality, i_dataLength, -1, -1, 0, 0, 0, 0, 0);
        }

        void init(
                const char *i_id, 
                unsigned i_idLength,
                const char *i_data, 
                const char *i_quality, 
                unsigned i_dataLength,
                unsigned i_originalAlignedLocation,
                unsigned i_originalMAPQ,
                unsigned i_originalSAMFlags,
                unsigned i_originalFrontClipping,
                unsigned i_originalBackClipping,
                unsigned i_originalFrontHardClipping,
                unsigned i_originalBackHardClipping)
        {
            id = i_id;
            idLength = i_idLength;
            data = unclippedData = i_data;
            quality = unclippedQuality = i_quality;
            dataLength = i_dataLength;
            unclippedLength = dataLength;
            frontClippedLength = 0;
            originalUnclippedDataBuffer = NULL;
            originalUnclippedQualityBuffer = NULL;
            clippingState = NoClipping;
            originalAlignedLocation = i_originalAlignedLocation;
            originalMAPQ = i_originalMAPQ;
            originalSAMFlags = i_originalSAMFlags;
            originalFrontClipping = i_originalFrontClipping;
            originalBackClipping = i_originalBackClipping;
            originalFrontHardClipping = i_originalFrontHardClipping;
            originalBackHardClipping = i_originalBackHardClipping;


            //
            // Check for lower case letters in the data, and convert to upper case if there are any.
            //
            unsigned anyLowerCase = 0;
            for (unsigned i = 0; i < dataLength; i++) {
                anyLowerCase |= IS_LOWER_CASE[data[i]];
            }

            if (anyLowerCase) {
                if (upperCaseDataBufferLength < dataLength) {
                    BigDealloc(upperCaseDataBuffer);
                    upperCaseDataBuffer = (char *)BigAlloc(dataLength);
                    upperCaseDataBufferLength = dataLength;
                }

                for (unsigned i = 0; i < dataLength; i++) {
                    upperCaseDataBuffer[i] = TO_UPPER_CASE[data[i]];
                }

                unclippedData = data = upperCaseDataBuffer;
            }
        }

        // For efficiency, this class holds id, data and quality pointers that are
        // *NOT* guaranteed to be to null-terminated strings; use the the length fields
        // to figure out how far to read into these strings.
        inline const char *getId() const {return id;}
        inline unsigned getIdLength() const {return idLength;}
        inline const char *getData() const {return data;}
        inline const char *getUnclippedData() const {return unclippedData;}
        inline const char *getQuality() const {return quality;}
        inline const char *getUnclippedQuality() const {return unclippedQuality;}
        inline unsigned getDataLength() const {return dataLength;}
        inline unsigned getUnclippedLength() const {return unclippedLength;}
        inline unsigned getFrontClippedLength() const {return (unsigned)(data - unclippedData);}    // number of bases clipped from the front of the read
        inline void setUnclippedLength(unsigned length) {unclippedLength = length;}
		inline ReadClippingType getClippingState() const {return clippingState;}
        inline DataBatch getBatch() { return batch; }
        inline void setBatch(DataBatch b) { batch = b; }
        inline const char* getReadGroup() const { return readGroup; }
        inline void setReadGroup(const char* rg) { readGroup = rg; }
        inline unsigned getOriginalAlignedLocation() {return originalAlignedLocation;}
        inline unsigned getOriginalMAPQ() {return originalMAPQ;}
        inline unsigned getOriginalSAMFlags() {return originalSAMFlags;}
        inline unsigned getOriginalFrontClipping() {return originalFrontClipping;}
        inline unsigned getOriginalBackClipping() {return originalBackClipping;}
        inline unsigned getOriginalFrontHardClipping() {return originalFrontHardClipping;}
        inline unsigned getOriginalBackHardClipping() {return originalBackHardClipping;}
        inline char* getAuxiliaryData(unsigned* o_length, bool * o_isSAM) const
        {
            *o_length = auxiliaryDataLength;
            *o_isSAM = auxiliaryData && auxiliaryDataLength >= 5 && auxiliaryData[2] == ':';
            return auxiliaryData;
        }
        inline void setAuxiliaryData(char* data, unsigned len)
        { auxiliaryData = data; auxiliaryDataLength = len; }

        void clip(ReadClippingType clipping, bool maintainOriginalClipping = false) {
            if (clipping == clippingState) {
                //
                // Already in the right state.
                //
                return;
            }

            //
            // Revert to unclipped, then clip to the correct state.
            //

            dataLength = unclippedLength;
            frontClippedLength = 0;
            data = unclippedData;
            quality = unclippedQuality;
            
            //
            // First clip from the back.
            //
            if (ClipBack == clipping || ClipFrontAndBack == clipping) {
                unsigned backClipping = 0;
                while (dataLength > 0 && quality[dataLength - 1] == '#') {
                    dataLength--;
                    backClipping++;
                }

                if (maintainOriginalClipping && backClipping < originalBackClipping) {
                    dataLength -= (originalBackClipping - backClipping);
                }
            }

            //
            // Then clip from the beginning.
            //
            if (ClipFront == clipping || ClipFrontAndBack == clipping) {
                frontClippedLength = 0;
                while (frontClippedLength < dataLength && quality[frontClippedLength] == '#') {
                    frontClippedLength++;
                }

                if (maintainOriginalClipping) {
                    frontClippedLength = max(frontClippedLength, originalFrontClipping);
                }
            }

            _ASSERT(frontClippedLength <= dataLength);

            dataLength -= frontClippedLength;
            data += frontClippedLength;
            quality += frontClippedLength;
 
            clippingState = clipping;
        };
        
        unsigned countOfTrailing2sInQuality() const {   // 2 here is represented in Phred+33, or ascii '#'
            unsigned count = 0;
            while (count < dataLength && quality[dataLength - 1 - count] == '#') {
                count++;
            }
            return count;
        }

        unsigned countOfNs() const {
            unsigned count = 0;
            for (unsigned i = 0; i < dataLength; i++) {
                count += IS_N[data[i]];
            }
            return count;
        }
        
        bool qualityFilter(float min_percent, unsigned min_qual, unsigned offset=33) const {        
            unsigned count = 0;
            for (unsigned i = 0; i < dataLength; i++) {
                if (quality[i]-offset >= min_qual) {
                    count++;
                }
            }
            if ((float(count)/float(dataLength))*100.f >= min_percent) {
                return true;
            }
            return false; 
        }

        void computeReverseCompliment(char *outputBuffer) { // Caller guarantees that outputBuffer is at least getDataLength() bytes
            for (unsigned i = 0; i < dataLength; i++) {
                outputBuffer[i] = COMPLEMENT[data[dataLength - i - 1]];
            }
        }

        void becomeRC()
        {
            if (NULL != originalUnclippedDataBuffer) {
                //
                // We've already RCed ourself.  Switch back.
                //
                unclippedData = originalUnclippedDataBuffer;
                unclippedQuality = originalUnclippedQualityBuffer;

                //
                // The clipping reverses as we go to/from RC.
                //
                frontClippedLength = unclippedLength - dataLength - frontClippedLength;
                data = unclippedData + frontClippedLength;
                quality = unclippedQuality + frontClippedLength;

                originalUnclippedDataBuffer = NULL;
                originalUnclippedQualityBuffer = NULL;

                return;
            }

            //
            // If we don't have local space to store the RC (and the quality RC) then allocate
            // it.  This preserves the buffer space (but not content) across calls to init().
            //
            if (unclippedLength > localBufferSize) {
                delete [] localUnclippedDataBuffer;
                delete [] localUnclippedQualityBuffer;

                localUnclippedDataBuffer = new char[unclippedLength];
                localUnclippedQualityBuffer =new char[unclippedLength];
                localBufferSize = unclippedLength;
            }

            //
            // We can't just call computeReverseCompliment() because we need to reverse the
            // unclipped portion of the buffer.
            //
            for (unsigned i = 0; i < unclippedLength; i++) {
                localUnclippedDataBuffer[i] = COMPLEMENT[unclippedData[unclippedLength - i - 1]];
                localUnclippedQualityBuffer[unclippedLength-i-1] = unclippedQuality[i];
            }

            originalUnclippedDataBuffer = unclippedData;
            originalUnclippedQualityBuffer = unclippedQuality;

            unclippedData = localUnclippedDataBuffer;
            unclippedQuality = localUnclippedQualityBuffer;

            //
            // The clipping reverses as we go to/from RC.
            //
            frontClippedLength = unclippedLength - dataLength - frontClippedLength;

            unsigned temp = originalFrontClipping;
            originalFrontClipping = originalBackClipping;
            originalBackClipping = temp;
            temp = originalFrontHardClipping;
            originalFrontHardClipping = originalBackHardClipping;
            originalBackHardClipping = temp;
            
            data = localUnclippedDataBuffer + frontClippedLength;
            quality = localUnclippedQualityBuffer + frontClippedLength;
        }


		static void checkIdMatch(Read* read0, Read* read1);

        static void computeClippingFromCigar(const char *cigarBuffer, unsigned *originalFrontClipping, unsigned *originalBackClipping, unsigned *originalFrontHardClipping, unsigned *originalBackHardClipping)
        {
            size_t cigarSize;
            const size_t cigarLimit = 1000;
            for (cigarSize = 0; cigarSize < cigarLimit && cigarBuffer[cigarSize] != '\0' && cigarBuffer[cigarSize] != '\t'; cigarSize++) {
                 // This loop body intentionally left blank.
            }

            if (cigarSize == cigarLimit) {
                fprintf(stderr, "Absurdly long cigar string.\n");
                soft_exit(1);
            }

            size_t frontHardClippingChars, backHardClippingChars, frontClippingChars, backClippingChars;

            //
            // Pull off the hard clipping first.
            //
            ExtractClipping(cigarBuffer, cigarSize, originalFrontHardClipping, originalBackHardClipping, 'H', &frontHardClippingChars, &backHardClippingChars);
            _ASSERT(frontHardClippingChars + backHardClippingChars <= cigarSize);

            //
            // Now look at what's left of the cigar string to see if there's soft clipping.
            //
            ExtractClipping(cigarBuffer + frontHardClippingChars, cigarSize - frontHardClippingChars - backHardClippingChars, originalFrontClipping, originalBackClipping,
                'S', &frontClippingChars, &backClippingChars);

        }

private:

        const char *id;
        const char *data;
        const char *unclippedData;
        const char *unclippedQuality;
        const char *quality;
        const char *readGroup;
        unsigned idLength;
        unsigned dataLength;
        unsigned unclippedLength;
        unsigned frontClippedLength;
        ReadClippingType clippingState;

        //
        // Alignment data that was in the read when it was read from a file.  While this should probably also be the place to put
        // information that'll be used by the read writer, for now it's not.  Hence, they're all called "original."
        //
        unsigned originalAlignedLocation;
        unsigned originalMAPQ;
        unsigned originalSAMFlags;
        unsigned originalFrontClipping;
        unsigned originalBackClipping;
        unsigned originalFrontHardClipping;
        unsigned originalBackHardClipping;

        //
        // Data stored in the read that's used if we RC ourself.  These are always unclipped.
        //
        char *localUnclippedDataBuffer;
        char *localUnclippedQualityBuffer;
        unsigned localBufferSize;

        //
        // If we've switched to RC, then the original data and quality are pointed to by these.
        // This allows us to switch back by just flipping pointers.
        //
        const char *originalUnclippedDataBuffer;
        const char *originalUnclippedQualityBuffer;

        // batch for managing lifetime during input
        DataBatch batch;

        //
        // If the read comes in with lower case letters, we need to convert it to upper case.  This buffer is used to store that
        // data.
        //
        char *upperCaseDataBuffer;
        unsigned upperCaseDataBufferLength;

        // auxiliary data in BAM or SAM format (can tell by looking at 3rd byte), if available
        char* auxiliaryData;
        unsigned auxiliaryDataLength;

        //
        // Pull the clipping info from the front and back of a cigar string.  
        static void ExtractClipping(const char *cigarBuffer, size_t cigarSize, unsigned *frontClipping, unsigned *backClipping, char clippingChar, size_t *frontClippingChars, size_t *backClippingChars)
        {

            *frontClipping = 0;
            const size_t bufferSize = 20;
            char buffer[bufferSize+1];  // +1 for trailing null
            unsigned i;
            for (i = 0; i < bufferSize && i < cigarSize && cigarBuffer[i] >= '0' && cigarBuffer[i] <= '9'; i++) {
                buffer[i] = cigarBuffer[i];
            }
            if (cigarBuffer[i] == clippingChar) {
                buffer[i] = '\0';
                *frontClipping = atoi(buffer);
                *frontClippingChars = i + 1;
            } else {
                *frontClippingChars = 0;
            }

            *backClipping = 0;
            *backClippingChars = 0;
            //
            // Find the end of the cigar string by looking for either the end of the string or a tab. Just start where we
            // were.
            //
            for (;i < cigarSize && cigarBuffer[i] != '\t' && cigarBuffer[i] != '\0'; i++) {
                // This loop body intentionally left blank.
            }
            if (i > 1 && cigarBuffer[i-1] == clippingChar) {
                for (i = i - 2; i >=0 && cigarBuffer[i] >= '0' && cigarBuffer[i] <= '9'; i--) {
                    // This loop body intentionally left blank.
                }
                //
                // If we've gotten back to the beginning of the string, then the whole thing is one big soft clip.  We arbitrarily
                // select that to be front clipping, and so leave the back clipping alone.
                if (i > 0) {
                    unsigned stringStart = i + 1;
                    for (i = stringStart; cigarBuffer[i] >= '0' && cigarBuffer[i] <= '9'; i++) {
                        buffer[i - stringStart] = cigarBuffer[i];
                    }
                    buffer[i - stringStart] = '\0';
                    *backClipping = atoi(buffer);
                    *backClippingChars = i - stringStart + 1;
                }
            }
        }
};

//
// Reads that copy the memory for their strings.  They're less efficient than the base
// Read class, but you can keep them around without holding references to the IO buffers
// and eventually stopping the IO.
//
class ReadWithOwnMemory : public Read {
public:
    ReadWithOwnMemory() : Read(), dataBuffer(NULL), idBuffer(NULL), qualityBuffer(NULL), auxBuffer(NULL) {}

    ReadWithOwnMemory(const Read &baseRead) {
        set(baseRead);
    }
    
    // must manually call destructor!
    void dispose() {
        delete [] dataBuffer;
        delete [] idBuffer;
        delete [] qualityBuffer;
        delete [] auxBuffer;
    }

private:

    void set(const Read &baseRead)
    {
        idBuffer = new char[baseRead.getIdLength()+1];
        dataBuffer = new char[baseRead.getUnclippedLength()+1];
        qualityBuffer = new char[baseRead.getUnclippedLength() + 1];

        memcpy(idBuffer,baseRead.getId(),baseRead.getIdLength());
        idBuffer[baseRead.getIdLength()] = '\0';    // Even though it doesn't need to be null terminated, it seems like a good idea.

        memcpy(dataBuffer,baseRead.getUnclippedData(),baseRead.getUnclippedLength());
        dataBuffer[baseRead.getUnclippedLength()] = '\0';

        memcpy(qualityBuffer,baseRead.getUnclippedQuality(),baseRead.getUnclippedLength());
        qualityBuffer[baseRead.getUnclippedLength()] = '\0';
    
        init(idBuffer,baseRead.getIdLength(),dataBuffer,qualityBuffer,baseRead.getUnclippedLength());
		clip(baseRead.getClippingState());

        setReadGroup(baseRead.getReadGroup());
        
        unsigned auxlen;
        bool auxsam;
        char* aux = baseRead.getAuxiliaryData(&auxlen, &auxsam);
        if (aux != NULL && auxlen > 0) {
            auxBuffer = new char[auxlen];
            memcpy(auxBuffer, aux, auxlen);
            setAuxiliaryData(auxBuffer, auxlen);
        } else {
            setAuxiliaryData(NULL, 0);
        }
    }
        
    char *idBuffer;
    char *dataBuffer;
    char *qualityBuffer;
    char *auxBuffer;
};
