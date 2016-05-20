/*++

Module Name:

    AlignmentAdjuster.h

Abstract:

    Header for code to adjust alignments to avoid indels at ends, or overhanging contig boundaries

Authors:

    Bill Bolosky, May, 2016

Environment:

User mode service.


Revision History:

    Pulled out of other places in SNAP

--*/

#pragma once
#include "GenomeIndex.h"
#include "LandauVishkin.h"
#include "Read.h"


class AlignmentAdjuster {
public:
    AlignmentAdjuster(const Genome *genome);

    void AdjustAlignment(Read *read, SingleAlignmentResult *result);
    void AdjustAlignments(Read **reads, PairedAlignmentResult *result);

private:

    LandauVishkinWithCigar lv;
    const Genome *genome;
};