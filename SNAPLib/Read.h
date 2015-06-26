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
#include "Error.h"
#include "Genome.h"
#include "AlignmentResult.h"

class FileFormat;

class Genome;

struct PairedAlignmentResult;



//#define LONG_READS
#ifdef LONG_READS
#define MAX_READ_LENGTH 400000
#else
#define MAX_READ_LENGTH 400
#endif

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
    const char*         defaultReadGroup;
    const char*         defaultReadGroupAux; // SAM or BAM depending on output format
    int                 defaultReadGroupAuxLen;
    ReadClippingType    clipping;
    bool                paired;
    bool                ignoreSecondaryAlignments;   // Should we just ignore reads with the Secondary Alignment bit set?
    bool                ignoreSupplementaryAlignments;  // Should we just ignore reads with the Supplementary Alignment bit set?
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

    // if you keep a read after the next call to getNextRead, you must call holdBatch
    // this increments the reference count to the batch
    virtual void holdBatch(DataBatch batch) = 0;

    // decremens hold refcount, when all holds are released the batch is no longer valid
    virtual bool releaseBatch(DataBatch batch) = 0;

    ReaderContext* getContext() { return &context; }

protected:
    ReaderContext context;
};

class PairedReadReader {
public:
    virtual ~PairedReadReader() {}

    // reading

    virtual bool getNextReadPair(Read *read1, Read *read2) = 0;
    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess) = 0;

    virtual void holdBatch(DataBatch batch) = 0;
    virtual bool releaseBatch(DataBatch batch) = 0;

    virtual ReaderContext* getContext() = 0;

    // wrap a single read source with a matcher that buffers reads until their mate is found
    static PairedReadReader* PairMatcher(ReadReader* single, bool quicklyDropUnpairedReads);
    static const int MatchBuffers = 2;
};

class ReadSupplier {
public:
    virtual Read *getNextRead() = 0;    // This read is valid until you call getNextRead, then it's done.  Don't worry about deallocating it.
    virtual ~ReadSupplier() {}

    virtual void holdBatch(DataBatch batch) = 0;
    virtual bool releaseBatch(DataBatch batch) = 0;
};

class PairedReadSupplier {
public:
    // These read are valid until you call getNextRead, then they're done.  Don't worry about deallocating them.
    virtual bool getNextReadPair(Read **read0, Read **read1) = 0;
    virtual ~PairedReadSupplier() {}

    virtual void holdBatch(DataBatch batch) = 0;
    virtual bool releaseBatch(DataBatch batch) = 0;
};

class ReadSupplierGenerator {
public:
    virtual ReadSupplier *generateNewReadSupplier() = 0;
    virtual ReaderContext* getContext() = 0;
    virtual ~ReadSupplierGenerator() {}
};

class PairedReadSupplierGenerator {
public:
    virtual PairedReadSupplier *generateNewPairedReadSupplier() = 0;
    virtual ReaderContext* getContext() = 0;
    virtual ~PairedReadSupplierGenerator() {}
};

class ReadWriter {
public:

    virtual ~ReadWriter() {}

    // write out header
	virtual bool writeHeader(const ReaderContext& context, bool sorted, int argc, const char **argv, const char *version, const char *rgLine, bool omitSQLines) = 0;

    //
    // write a batch of single reads, the first one of which is a primary alignment and the rest secondary.
    //
    virtual bool writeReads(const ReaderContext& context, Read *read, SingleAlignmentResult *results, int nResults, bool firstIsPrimary) = 0;

    //
    // Write a batch of paired alignments, including some single secondary alignments.  reads needs to be exactly two reads, singleAlignmentResult is a pointer to two arrays of
    // SingleAlignmentResults for single alignments of the respective reads, and the number of results is given by the two ints pointed to by nSingleResults.  The first paired
    // result is primary, all others are secondary.
    //
    virtual bool writePairs(const ReaderContext& context, Read **reads /* array of size 2 */, PairedAlignmentResult *result, int nResults,
        SingleAlignmentResult **singleResults /* array of size 2*/, int *nSingleResults /* array of size 2*/, bool firstIsPrimary) = 0;


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
        const Genome* genome);
};

#define READ_GROUP_FROM_AUX     ((const char*) -1)
    
class Read {
public:
        Read() :    
            id(NULL), data(NULL), quality(NULL), 
            localBufferAllocationOffset(0),
            clippingState(NoClipping), currentReadDirection(FORWARD),
            upcaseForwardRead(NULL), auxiliaryData(NULL), auxiliaryDataLength(0),
            readGroup(NULL), originalAlignedLocation(-1), originalMAPQ(-1), originalSAMFlags(0),
            originalFrontClipping(0), originalBackClipping(0), originalFrontHardClipping(0), originalBackHardClipping(0),
            originalRNEXT(NULL), originalRNEXTLength(0), originalPNEXT(0), additionalFrontClipping(0)
        {}

        Read(const Read& other) :  localBufferAllocationOffset(0)
        {
            copyFromOtherRead(other);
        }

        ~Read()
        {
        }

        void dispose()
        {
            localBufferAllocationOffset = 0;
            data = quality = unclippedData = unclippedQuality = externalData = NULL;
         }

        void operator=(const Read& other)
        {
            copyFromOtherRead(other);
        }

        void copyFromOtherRead(const Read& other)
        {
            id = other.id;
            idLength = other.idLength;
            frontClippedLength = other.frontClippedLength;
            dataLength = other.dataLength;
            externalData = other.externalData;
            externalQuality = other.externalQuality;
            currentReadDirection = other.currentReadDirection;
            localBufferAllocationOffset = 0;    // Clears out any allocations that might previously have been in the buffer
            upcaseForwardRead = rcData = rcQuality = NULL;
            unclippedLength = other.unclippedLength;

            if (other.localBufferAllocationOffset != 0) {
                //
                // Copy the other read's local buffer to us.
                //
                assureLocalBufferLargeEnough();
                _ASSERT(other.localBufferAllocationOffset <= localBufferLength);
                memcpy(localBuffer, other.localBuffer, other.localBufferAllocationOffset);
                localBufferAllocationOffset = other.localBufferAllocationOffset;

                if (NULL != other.upcaseForwardRead) {
                    //
                    // Assert that it's in the other read's local buffer.
                    //
                    _ASSERT(other.upcaseForwardRead >= other.localBuffer && other.upcaseForwardRead <= other.localBuffer + other.localBufferAllocationOffset - unclippedLength);

                    //
                    // And put ours at the same offset in our local buffer.
                    //
                    upcaseForwardRead = localBuffer + (other.upcaseForwardRead - other.localBuffer);
                }

                if (NULL != other.rcData) {
                    //
                    // Assert that it's in the other read's local buffer.
                    //
                    _ASSERT(other.rcData >= other.localBuffer && other.rcData <= other.localBuffer + other.localBufferAllocationOffset - unclippedLength);

                    //
                    // And put ours at the same offset in our local buffer.
                    //
                    rcData = localBuffer + (other.rcData - other.localBuffer);

                    //
                    // And the same for RC quality.
                    //
                    _ASSERT(other.rcQuality >= other.localBuffer && other.rcQuality <= other.localBuffer + other.localBufferAllocationOffset - unclippedLength);
                    rcQuality = localBuffer + (other.rcQuality - other.localBuffer);
                } else {
                    _ASSERT(NULL == other.rcQuality);
                }
            } else {
                _ASSERT(other.upcaseForwardRead == NULL && other.rcData == NULL && other.rcQuality == NULL);
            }

            //
            // Now set up the data, unclippedData, quality and unclippedQuality pointers.
            //
            if (NULL == other.localBuffer || other.data < other.localBuffer || other.data >= other.localBuffer + other.localBufferAllocationOffset - dataLength) {
                //
                // Not in the other read's local buffer, so it must be external.  Copy it.
                //
                data = other.data;
                _ASSERT(NULL == other.localBuffer || other.quality < other.localBuffer || other.quality >= other.localBuffer + other.localBufferAllocationOffset);
                quality = other.quality;
                _ASSERT(NULL == other.localBuffer || other.unclippedData < other.localBuffer || other.unclippedData >= other.localBuffer + other.localBufferAllocationOffset);
                unclippedData = other.unclippedData;
                _ASSERT(NULL == other.localBuffer || other.unclippedQuality < other.localBuffer || other.unclippedQuality >= other.localBuffer + other.localBufferAllocationOffset);
                unclippedQuality = other.unclippedQuality;
            } else {
                //
                // It is in the other read's local buffer.  Copy the local buffer offsets from the other read into this one.
                //
                data = localBuffer + (other.data - other.localBuffer);
                _ASSERT(other.quality >= other.localBuffer && other.quality <= other.localBuffer + other.localBufferAllocationOffset - dataLength);
                quality = localBuffer + (other.quality - other.localBuffer);
                _ASSERT(other.unclippedData >= other.localBuffer && other.unclippedData <= other.localBuffer + other.localBufferAllocationOffset - unclippedLength);
                unclippedData = localBuffer + (other.unclippedData - other.localBuffer);
                _ASSERT(other.unclippedQuality >= other.localBuffer && other.unclippedQuality <= other.localBuffer + other.localBufferAllocationOffset - unclippedLength);
                unclippedQuality = localBuffer + (other.unclippedQuality - other.localBuffer);
            }

            clippingState = other.clippingState;
            batch = other.batch;
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
            originalRNEXT = other.originalRNEXT;
            originalRNEXTLength = other.originalRNEXTLength;
            originalPNEXT = other.originalPNEXT;
            additionalFrontClipping = other.additionalFrontClipping;
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
            init(i_id, i_idLength, i_data, i_quality, i_dataLength, InvalidGenomeLocation, -1, 0, 0, 0, 0, 0, NULL, 0, 0);
        }

        void init(
                const char *        i_id, 
                unsigned            i_idLength,
                const char *        i_data, 
                const char *        i_quality, 
                unsigned            i_dataLength,
                GenomeLocation      i_originalAlignedLocation,
                unsigned            i_originalMAPQ,
                unsigned            i_originalSAMFlags,
                unsigned            i_originalFrontClipping,
                unsigned            i_originalBackClipping,
                unsigned            i_originalFrontHardClipping,
                unsigned            i_originalBackHardClipping,
                const char *        i_originalRNEXT,
                unsigned            i_originalRNEXTLength,
                unsigned            i_originalPNEXT,
                bool                allUpper = false)
        {
            id = i_id;
            idLength = i_idLength;
            data = unclippedData = externalData = i_data;
            quality = unclippedQuality = externalQuality = i_quality;
            dataLength = i_dataLength;
            unclippedLength = dataLength;
            frontClippedLength = 0;
            clippingState = NoClipping;
            additionalFrontClipping = 0;
            originalAlignedLocation = i_originalAlignedLocation;
            originalMAPQ = i_originalMAPQ;
            originalSAMFlags = i_originalSAMFlags;
            originalFrontClipping = i_originalFrontClipping;
            originalBackClipping = i_originalBackClipping;
            originalFrontHardClipping = i_originalFrontHardClipping;
            originalBackHardClipping = i_originalBackHardClipping;
            originalRNEXT = i_originalRNEXT;
            originalRNEXTLength = i_originalRNEXTLength;
            originalPNEXT = i_originalPNEXT;
            currentReadDirection = FORWARD;

            localBufferAllocationOffset = 0;    // Clears out any allocations that might previously have been in the buffer
            upcaseForwardRead = rcData = rcQuality = NULL;

            //
            // Check for lower case letters in the data, and convert to upper case if there are any.  Also convert
            // '.' to N.
            //
            if (! allUpper) {
                unsigned anyLowerCase = 0;
                for (unsigned i = 0; i < dataLength; i++) {
                    anyLowerCase |= IS_LOWER_CASE_OR_DOT[data[i]];
                }

                if (anyLowerCase) {
                    assureLocalBufferLargeEnough();
                    upcaseForwardRead = localBuffer;
                    localBufferAllocationOffset += unclippedLength;
                    for (unsigned i = 0; i < dataLength; i++) {
                        upcaseForwardRead[i] = TO_UPPER_CASE_DOT_TO_N[data[i]];
                    }

                    unclippedData = data = upcaseForwardRead;
                }
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
		inline unsigned getBackClippedLength() const {return unclippedLength - dataLength - getFrontClippedLength();}
        inline void setUnclippedLength(unsigned length) {unclippedLength = length;}
		inline ReadClippingType getClippingState() const {return clippingState;}
        inline DataBatch getBatch() { return batch; }
        inline void setBatch(DataBatch b) { batch = b; }
        inline const char* getReadGroup() const { return readGroup; }
        inline void setReadGroup(const char* rg) { readGroup = rg; }
        inline GenomeLocation getOriginalAlignedLocation() {return originalAlignedLocation;}
        inline unsigned getOriginalMAPQ() {return originalMAPQ;}
        inline unsigned getOriginalSAMFlags() {return originalSAMFlags;}
        inline unsigned getOriginalFrontClipping() {return originalFrontClipping;}
        inline unsigned getOriginalBackClipping() {return originalBackClipping;}
        inline unsigned getOriginalFrontHardClipping() {return originalFrontHardClipping;}
        inline unsigned getOriginalBackHardClipping() {return originalBackHardClipping;}
        inline const char *getOriginalRNEXT() {return originalRNEXT;}
        inline unsigned getOriginalRNEXTLength() {return originalRNEXTLength;}
        inline unsigned getOriginalPNEXT() {return originalPNEXT;}
        inline void setAdditionalFrontClipping(int clipping)
        {
            data += clipping - additionalFrontClipping;
            dataLength -= clipping - additionalFrontClipping;
            quality += clipping - additionalFrontClipping;
            additionalFrontClipping = clipping;
        }

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

        void computeReverseCompliment(char *outputBuffer) { // Caller guarantees that outputBuffer is at least getDataLength() bytes
            for (unsigned i = 0; i < dataLength; i++) {
                outputBuffer[i] = COMPLEMENT[data[dataLength - i - 1]];
            }
        }

        void becomeRC()
        {
            if (RC == currentReadDirection) {
                //
                // We've already RCed ourself.  Switch back.
                //
                if (NULL != upcaseForwardRead) {
                    unclippedData = upcaseForwardRead;
                  } else {
                    unclippedData = externalData;
                }
                unclippedQuality = externalQuality;

                currentReadDirection = FORWARD;
            } else {

                if (rcData != NULL) {
                    //
                    // We've already been RC, just switch back.
                    //
                    unclippedData = rcData;
                    unclippedQuality = rcQuality;
                } else {
                    assureLocalBufferLargeEnough();
                    rcData = localBuffer + localBufferAllocationOffset;
                    localBufferAllocationOffset += unclippedLength;
                    rcQuality = localBuffer + localBufferAllocationOffset;
                    localBufferAllocationOffset += unclippedLength;

                    _ASSERT(localBufferAllocationOffset <= localBufferLength);

                    for (unsigned i = 0; i < unclippedLength; i++) {
                        rcData[i] = COMPLEMENT[unclippedData[unclippedLength - i - 1]];
                        rcQuality[unclippedLength-i-1] = unclippedQuality[i];
                    }

                    unclippedData = rcData;
                    unclippedQuality = rcQuality;
                }

                currentReadDirection = RC;
            }

            //
            // The clipping reverses as we go to/from RC.
            //
            frontClippedLength = unclippedLength - dataLength - frontClippedLength;
            data = unclippedData + frontClippedLength;
            quality = unclippedQuality + frontClippedLength;

            unsigned temp = originalFrontClipping;
            originalFrontClipping = originalBackClipping;
            originalBackClipping = temp;
            temp = originalFrontHardClipping;
            originalFrontHardClipping = originalBackHardClipping;
            originalBackHardClipping = temp;
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
                WriteErrorMessage( "Absurdly long cigar string.\n");
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
        int additionalFrontClipping;

        //
        // Alignment data that was in the read when it was read from a file.  While this should probably also be the place to put
        // information that'll be used by the read writer, for now it's not.  Hence, they're all called "original."
        //
        GenomeLocation originalAlignedLocation;
        unsigned originalMAPQ;
        unsigned originalSAMFlags;
        unsigned originalFrontClipping;
        unsigned originalBackClipping;
        unsigned originalFrontHardClipping;
        unsigned originalBackHardClipping;
        const char *originalRNEXT;
        unsigned originalRNEXTLength;
        unsigned originalPNEXT;


        //
        // Memory that's local to this read and that is used to contain an upcased version of the read as well as 
        // RC read & quality strings.  It survives init() so as to avoid memory allocation overhead.
        //
        char localBuffer[MAX_READ_LENGTH * 3];
        static const unsigned localBufferLength;
        unsigned localBufferAllocationOffset;   // The next location to allocate in the local buffer.
        char *upcaseForwardRead;                // Either NULL or points into localBuffer.  Used when the incoming read isn't all capitalized.  Unclipped.
        char *rcData;                           // Either NULL or points into localBuffer.  Used when we've computed a reverse complement of the read, whether we're using it or not.  Unclipped.
        char *rcQuality;                        // Ditto for quality.
        const char *externalData;               // The data that was passed in at init() time, memory doesn't belong to this.
        const char *externalQuality;            // The quality that was passed in at init() time, memory doens't belong to this.
        Direction currentReadDirection;

        inline void assureLocalBufferLargeEnough()
        {
#if 0   // Always true with static allocation
            if (localBufferLength < 3 * unclippedLength) {
                _ASSERT(0 == localBufferAllocationOffset);  // Can only do this when the buffer is empty
                if (NULL != localBuffer) {
                    BigDealloc(localBuffer);
                }
                localBufferLength = RoundUpToPageSize(3 * unclippedLength);
                localBuffer = (char *)BigAlloc(localBufferLength);
            }

#endif // 0
        }

        // batch for managing lifetime during input
        DataBatch batch;

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
    ReadWithOwnMemory() : Read(), extraBuffer(NULL), dataBuffer(NULL), idBuffer(NULL), qualityBuffer(NULL), auxBuffer(NULL) {}

    ReadWithOwnMemory(const Read &baseRead) {
        set(baseRead);
    }
    
    // must manually call destructor!
    void dispose() {
        if (extraBuffer != NULL) {
            delete [] extraBuffer;
        }
    }

private:

    void set(const Read &baseRead)
    {
        // allocate space in ownBuffer if possible; id/aux might need extraBuffer
        dataBuffer = ownBuffer;
        int ownBufferUsed = baseRead.getUnclippedLength() + 1;
        qualityBuffer = ownBuffer + ownBufferUsed;
        ownBufferUsed += baseRead.getUnclippedLength() + 1;
        unsigned auxLen;
        bool auxSam;
        char* aux = baseRead.getAuxiliaryData(&auxLen, &auxSam);
        if (baseRead.getIdLength() + 1 < sizeof(ownBuffer) - ownBufferUsed) {
            idBuffer = ownBuffer + ownBufferUsed;
            ownBufferUsed += baseRead.getIdLength() + 1;
        } else {
            idBuffer = NULL;
        }
        if (auxLen > 0 && auxLen < sizeof(ownBuffer) - ownBufferUsed) {
            auxBuffer = ownBuffer + ownBufferUsed;
            ownBufferUsed += auxLen;
        } else {
            auxBuffer = NULL;
        }
        if (idBuffer == NULL || (auxLen > 0 && auxBuffer == NULL)) {
            extraBuffer = new char[(idBuffer == NULL ? baseRead.getIdLength() + 1 : 0) + auxLen];
            int extraBufferUsed = 0;
            if (idBuffer == NULL) {
                idBuffer = extraBuffer;
                extraBufferUsed += baseRead.getIdLength() + 1;
            }
            if (auxLen > 0 && auxBuffer == NULL) {
                auxBuffer = extraBuffer + extraBufferUsed;
            }
        } else {
            extraBuffer = NULL;
        }

        // copy data into buffers
        memcpy(idBuffer,baseRead.getId(),baseRead.getIdLength());
        idBuffer[baseRead.getIdLength()] = '\0';    // Even though it doesn't need to be null terminated, it seems like a good idea.

        memcpy(dataBuffer,baseRead.getUnclippedData(),baseRead.getUnclippedLength());
        dataBuffer[baseRead.getUnclippedLength()] = '\0';

        memcpy(qualityBuffer,baseRead.getUnclippedQuality(),baseRead.getUnclippedLength());
        qualityBuffer[baseRead.getUnclippedLength()] = '\0';
    
        init(idBuffer,baseRead.getIdLength(),dataBuffer,qualityBuffer,baseRead.getUnclippedLength());
		clip(baseRead.getClippingState());

        setReadGroup(baseRead.getReadGroup());
        
        if (aux != NULL && auxLen > 0) {
            memcpy(auxBuffer, aux, auxLen);
            setAuxiliaryData(auxBuffer, auxLen);
        } else {
            setAuxiliaryData(NULL, 0);
        }
    }
        
    char ownBuffer[MAX_READ_LENGTH * 2 + 1000]; // internal buffer for copied data
    char* extraBuffer; // extra buffer if internal buffer not big enough

    // should all point into ownBuffer or extraBuffer
    char *idBuffer;
    char *dataBuffer;
    char *qualityBuffer;
    char *auxBuffer;
};

extern const unsigned DEFAULT_MIN_READ_LENGTH;
