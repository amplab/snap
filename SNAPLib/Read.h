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

class Read;

enum ReadClippingType {NoClipping, ClipFront, ClipBack, ClipFrontAndBack};

class ReadReader {
public:
        virtual ~ReadReader() {}
        virtual void readDoneWithBuffer(unsigned *referenceCount) = 0;
        virtual bool getNextRead(Read *readToUpdate) = 0;
        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess) = 0;
};

class PairedReadReader {
public:
        virtual ~PairedReadReader() {}
        virtual bool getNextReadPair(Read *read1, Read *read2) = 0;
        virtual ReadReader *getReaderToInitializeRead(int whichHalfOfPair) = 0;
        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess) = 0;
};

    
class Read {
public:
        Read(ReadReader *i_reader = NULL) : reader(i_reader), id(NULL), data(NULL), quality(NULL), 
                                            localUnclippedDataBuffer(NULL), localUnclippedQualityBuffer(NULL), localBufferSize(0),
                                            originalUnclippedDataBuffer(NULL), originalUnclippedQualityBuffer(NULL), clippingState(NoClipping)
        {
            referenceCounts[0] = referenceCounts[1] = NULL;
        }

        ~Read()
        {
            delete [] localUnclippedDataBuffer;
            delete [] localUnclippedQualityBuffer;
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
                unsigned i_dataLength, 
                unsigned **newReferenceCounts)
        {
            commonInit(i_id,i_idLength,i_data,i_quality,i_dataLength);

            if (NULL != referenceCounts[0]) {
                (*referenceCounts[0])--;
                if (0 == *referenceCounts[0]) {
                    reader->readDoneWithBuffer(referenceCounts[0]);
                }
            }

            if (NULL != referenceCounts[1]) {
                (*referenceCounts[1])--;
                if (0 == *referenceCounts[1]) {
                    reader->readDoneWithBuffer(referenceCounts[1]);
                }
            }

            if (NULL != newReferenceCounts) {
                referenceCounts[0] = newReferenceCounts[0];
                referenceCounts[1] = newReferenceCounts[1];
            }
        }

        void init(
                const char *i_id, 
                unsigned i_idLength,
                const char *i_data, 
                const char *i_quality, 
                unsigned i_dataLength)
        {
            commonInit(i_id,i_idLength,i_data,i_quality,i_dataLength);

            referenceCounts[0] = referenceCounts[1] = NULL;
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

        void clip(ReadClippingType clipping) {
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
                while (dataLength > 0 && quality[dataLength - 1] == '#') {
                    dataLength--;
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
            }

            if (dataLength-frontClippedLength < 50) {
                // There are a lot of quality-2 bases in the read so just use all of it (IS THIS A GOOD IDEA? -BB)
                dataLength = unclippedLength;
                frontClippedLength = 0;
            } else {
                dataLength -= frontClippedLength;
                data += frontClippedLength;
                quality += frontClippedLength;
            }

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
            
            data = localUnclippedDataBuffer + frontClippedLength;
            quality = localUnclippedQualityBuffer + frontClippedLength;
        }

private:

        void commonInit(
                const char *i_id, 
                unsigned i_idLength,
                const char *i_data, 
                const char *i_quality, 
                unsigned i_dataLength)
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
        }

        const char *id;
        const char *data;
        const char *unclippedData;
        const char *unclippedQuality;
        const char *quality;
        unsigned idLength;
        unsigned dataLength;
        unsigned unclippedLength;
        unsigned frontClippedLength;
        ReadClippingType clippingState;

        ReadReader *reader;
        unsigned *referenceCounts[2];

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
};

//
// Reads that copy the memory for their strings.  They're less efficient than the base
// Read class, but you can keep them around without holding references to the IO buffers
// and eventually stopping the IO.
//
class ReadWithOwnMemory : public Read {
public:
        ReadWithOwnMemory(const Read &baseRead) {
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
        }

        ~ReadWithOwnMemory() {
            delete [] dataBuffer;
            delete [] idBuffer;
            delete [] qualityBuffer;
        }

private:

    char *idBuffer;
    char *dataBuffer;
    char *qualityBuffer;
        
};
