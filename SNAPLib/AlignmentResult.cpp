/*++

Module Name:

AlignmentResult.cpp

Abstract:

Code for SNAP genome alignment results

Authors:

Bill Bolosky, March, 2015

Environment:

Revision History:


--*/

#include "stdafx.h"
#include "AlignmentResult.h"
#include "GenomeIndex.h"



    int 
SingleAlignmentResult::compareByContigAndScore(const void *first_, const void *second_)
{
    extern GenomeIndex *g_index;    // Sorry, but no easy way to get this into here

    const SingleAlignmentResult *first = (SingleAlignmentResult *)first_;
    const SingleAlignmentResult *second = (SingleAlignmentResult *)second_;

    int firstContig = g_index->getGenome()->getContigNumAtLocation(first->location);
    int secondContig = g_index->getGenome()->getContigNumAtLocation(second->location);

    if (firstContig < secondContig) {
        return -1;
    } else if (firstContig > secondContig) {
        return 1;
    } else if (first->score < second->score) {
        return -1;
    } else if (first->score > second->score) {
        return 1;
    } else {
        return 0;
    }
 }

int
    SingleAlignmentResult::compareByScore(const void *first_, const void *second_)
{
    const SingleAlignmentResult *first = (SingleAlignmentResult *)first_;
    const SingleAlignmentResult *second = (SingleAlignmentResult *)second_;

    if (first->score < second->score) {
        return -1;
    } else if (first->score > second->score) {
        return 1;
    } else {
        return 0;
    }
}

    int
PairedAlignmentResult::compareByContigAndScore(const void *first_, const void *second_)
{
    extern GenomeIndex *g_index;    // Sorry, but no easy way to get this into here

    const PairedAlignmentResult *first = (PairedAlignmentResult *)first_;
    const PairedAlignmentResult *second = (PairedAlignmentResult *)second_;

    int firstContig = g_index->getGenome()->getContigNumAtLocation(first->location[0]);
    int secondContig = g_index->getGenome()->getContigNumAtLocation(second->location[0]);

    if (firstContig < secondContig) {
        return -1;
    } else if (firstContig > secondContig) {
        return 1;
    } else if (first->score < second->score) {
        return -1;
    } else if (first->score > second->score) {
        return 1;
    } else {
        return 0;
    }
}

int
PairedAlignmentResult::compareByScore(const void *first_, const void *second_)
{
    const PairedAlignmentResult *first = (PairedAlignmentResult *)first_;
    const PairedAlignmentResult *second = (PairedAlignmentResult *)second_;

    int firstScore = first->score[0] + first->score[1];
    int secondScore = second->score[0] + second->score[1];

    if (firstScore < secondScore) {
        return -1;
    } else if (firstScore > secondScore) {
        return 1;
    } else {
        return 0;
    }
}