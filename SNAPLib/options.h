/*++

Module Name:

    options.h

Abstract:

    Compile time options for cSNAP

Authors:

    Bill Bolosky, September, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once

#ifndef USE_DEVTEAM_OPTIONS
#define USE_DEVTEAM_OPTIONS         0   // Options only for the development team, not for the release version
#endif  // USE_DEVTEAM_OPTIONS

#ifndef INSTRUMENTATION_FOR_PAPER
#define INSTRUMENTATION_FOR_PAPER   0   // Turn this on to generate raw data about hit sets and their intersections for the paper
#endif  // INSTRUMENTATION_FOR_PAPER

#ifndef HIT_DEPTH_COUNTING
#define HIT_DEPTH_COUNTING          0   // This is a special piece of code that generates information about seeding for an entire genome.
#endif  // HIT_DEPTH_COUNTING

#ifndef TIME_HISTOGRAM
#define TIME_HISTOGRAM              0   // Generate histograms of reads by time to align, by MAPQ and by final edit distance (NM score)
#endif  // TIME_HISTOGRAM