/*++

Module Name:

    AlignmentAdjuster.cpp

Abstract:

    Code to adjust alignments to avoid indels at ends, or overhanging contig boundaries

Authors:

    Bill Bolosky, May, 2016

Environment:

    User mode service.

Revision History:

    Pulled out of other places in SNAP

--*/

#include "stdafx.h"
#include "AlignmentAdjuster.h"

AlignmentAdjuster::AlignmentAdjuster(const Genome *_genome)
{
    genome = _genome;
}

    void 
AlignmentAdjuster::AdjustAlignment(Read *read, SingleAlignmentResult *result)
{
    int additionalFrontClipping = 0;

    if (result->status == NotFound) {
        //
        // No adjustment required on unaligned reads.
        //
        return;
    }

    read->setAdditionalFrontClipping(0);
    int cumulativeAddFrontClipping = 0;

    unsigned nAdjustments = 0;

    const int cigarBufLen = 1000;
    char cigarBuf[cigarBufLen + 1];
    char dataBuffer[MAX_READ_LENGTH];
    const char *data = dataBuffer;
    const char *clippedData;
    int netIndel = 0;
    GenomeDistance dataLength;
    const Genome::Contig *contig;
    GenomeDistance extraBasesClippedAfter;
    const char *reference;
    int cigarBufUsed = 0;

    for (;;) {
        // Write the data and quality strings. If the read is reverse complemented, these need to
        // be backwards from the original read. Also, both need to be unclipped.
        GenomeDistance fullLength = read->getUnclippedLength();
        GenomeDistance basesClippedBefore;
        GenomeDistance basesClippedAfter;
        dataLength = read->getDataLength();

        if (fullLength > MAX_READ_LENGTH) {
            WriteErrorMessage("AlignmentAdjuster: read too long (how did it get past the reader?), %d > %d\n", fullLength, MAX_READ_LENGTH);
            soft_exit(1);
        }


        if (result->direction == RC) {
            for (unsigned i = 0; i < fullLength; i++) {
                dataBuffer[fullLength - 1 - i] = COMPLEMENT[read->getUnclippedData()[i]];
            }
            clippedData = &data[fullLength - dataLength - read->getFrontClippedLength()];
            basesClippedBefore = fullLength - dataLength - read->getFrontClippedLength();
            basesClippedAfter = read->getFrontClippedLength();
        } else {
            memcpy(dataBuffer, read->getUnclippedData(), read->getUnclippedLength());
            clippedData = read->getData();
            basesClippedBefore = read->getFrontClippedLength();
            basesClippedAfter = fullLength - dataLength - basesClippedBefore;
        }

        GenomeDistance extraBasesClippedBefore = 0;

        contig = genome->getContigForRead(result->location, read->getDataLength(), &extraBasesClippedBefore);
        _ASSERT(NULL != contig && contig->length > genome->getChromosomePadding());

        GenomeLocation genomeLocation = result->location + extraBasesClippedBefore;

        clippedData += extraBasesClippedBefore;
        dataLength -= extraBasesClippedBefore;


        if (genomeLocation + dataLength > contig->beginningLocation + contig->length - genome->getChromosomePadding()) {
            //
            // The read hangs off the end of the contig.  Soft clip it at the end.  This is a tentative amount that assumes no net indels in the
            // mapping, we'll refine it later if needed.
            //
            extraBasesClippedAfter = genomeLocation + dataLength - (contig->beginningLocation + contig->length - genome->getChromosomePadding());
        } else {
            extraBasesClippedAfter = 0;
        }

        reference = genome->getSubstring(genomeLocation, dataLength - extraBasesClippedAfter);
        _ASSERT(NULL != reference); // Since we just clipped it to fit.
        if (NULL == reference) {    // This shouldn't happen, but the old code checked it, so I'll leave it in.
            result->status = NotFound;
            result->location = InvalidGenomeLocation;
            return;
        }


        result->score = lv.computeEditDistanceNormalized(
            reference,
            (int)(dataLength - extraBasesClippedAfter + MAX_K), // Add space incase of indels.  We know there's enough, because the reference is padded.
            clippedData,
            (int)(dataLength - extraBasesClippedAfter),
            MAX_K - 1,
            cigarBuf,
            cigarBufLen,
            true,   // useM (doesn't matter, since we're throwing the cigar string away anyway)
            COMPACT_CIGAR_STRING,
            &cigarBufUsed,
            &additionalFrontClipping,
            &netIndel);

        if (0 == additionalFrontClipping) {
            break;
        }

        nAdjustments++;

        // redo if read modified (e.g. to add soft clipping, or move alignment for a leading I.
        const Genome::Contig *originalContig = genome->getContigAtLocation(result->location);
        const Genome::Contig *newContig = genome->getContigAtLocation(result->location + additionalFrontClipping);
        if (newContig == NULL || newContig != originalContig || result->location + additionalFrontClipping > originalContig->beginningLocation + originalContig->length - genome->getChromosomePadding() ||
            nAdjustments > read->getDataLength()) {
            //
            // Altering this would push us over a contig boundary, or we're stuck in a loop.  Just give up on the read.
            //
            result->status = NotFound;
            result->location = InvalidGenomeLocation;
            return;
        } else {
            cumulativeAddFrontClipping += additionalFrontClipping;
            read->setAdditionalFrontClipping(__max(0,cumulativeAddFrontClipping));
            result->clippingForReadAdjustment = __max(0, cumulativeAddFrontClipping);
            result->location = result->location + additionalFrontClipping;
        }
    }  // for ever

    //
    // Code copied from SAM.cpp:computeCigar()
    //
    // Normally, we'd be done.  However, if the amount that we would clip at the end of the read because of hanging off of the end
    // of the contig changed, then we need to recompute.  In some cases this is an iterative processess as we add or remove bits
    // of read.  
    //
    GenomeDistance newExtraBasesClippedAfter = __max(0, result->location + dataLength + netIndel - (contig->beginningLocation + contig->length - genome->getChromosomePadding()));
    for (GenomeDistance pass = 0; pass < dataLength; pass++) {
        if (newExtraBasesClippedAfter == extraBasesClippedAfter) {
            return;
        }

        extraBasesClippedAfter = newExtraBasesClippedAfter;

        result->score = lv.computeEditDistanceNormalized(
            reference,
            (int)(dataLength - extraBasesClippedAfter + MAX_K), // Add space incase of indels.  We know there's enough, because the reference is padded.
            data,
            (int)(dataLength - extraBasesClippedAfter),
            MAX_K - 1,
            cigarBuf,
            cigarBufLen,
            true, // useM (doesn't matter, since we're throwing the cigar string away anyway)
            COMPACT_CIGAR_STRING,
            &cigarBufUsed,
            &additionalFrontClipping,
            &netIndel);

        newExtraBasesClippedAfter = __max(0, result->location + dataLength + netIndel - (contig->beginningLocation + contig->length - genome->getChromosomePadding()));
    }

    _ASSERT(!"cigar computation didn't converge");
}

    
    void 
AlignmentAdjuster::AdjustAlignments(Read **reads, PairedAlignmentResult *result)
{
    for (int i = 0; i < NUM_READS_PER_PAIR; i++) {
        SingleAlignmentResult singleResult;
        singleResult.direction = result->direction[i];
        singleResult.location = result->location[i];
        singleResult.mapq = result->mapq[i];
        singleResult.score = result->score[i];
        singleResult.status = result->status[i];
        singleResult.clippingForReadAdjustment = 0;

        AdjustAlignment(reads[i], &singleResult);

        result->direction[i] = singleResult.direction;
        result->location[i] = singleResult.location;
        result->mapq[i] = singleResult.mapq;
        result->score[i] = singleResult.score;
        result->status[i] = singleResult.status;
        result->clippingForReadAdjustment[i] = singleResult.clippingForReadAdjustment;
    }
}

