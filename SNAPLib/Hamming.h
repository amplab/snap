/*++

Module Name:

    Hamming.h

Abstract:

    Hamming distance computation

Authors:

    Bill Bolosky, January, 2013

Environment:

    User mode service.

Revision History:

--*/

#pragma once

int ComputeHammingDistance(const char *text, const char *pattern, size_t len, const char *quality, int scoreLimit, double *mapProbability);
void InitializeHammingDistance(double snpProbability);