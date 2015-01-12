/*++

Module Name:

    SeedSequencer.h

Abstract:

    Code for determining the order of seeds in a read.

Authors:

    Bill Bolosky, February, 2013

Environment:

    User mode service.

Revision History:

    Factored out of BaseAligner

--*/

#pragma once
#include "exit.h"
#include "Seed.h"
#include "Error.h"


class SeedSequencer {
public:
    SeedSequencer(unsigned i_seedSize);

    inline unsigned SeedOffset(unsigned wrapCount) {
        _ASSERT(wrapCount < seedSize);
        return offsets[wrapCount];
    }

    inline unsigned GetWrappedNextSeedToTest(unsigned wrapCount) {
        _ASSERT(wrapCount < seedSize);
        return(offsets[wrapCount]);
    }

private:

    unsigned seedSize;
    unsigned *offsets;
};



void InitializeSeedSequencers();

unsigned GetWrappedNextSeedToTest(unsigned seedLen, unsigned wrapCount) 
#if     1
        ;   // The definition is in SeedSequencer.cpp
#else   // 1 -- old style with switch/case
{
    if (0 == wrapCount) {
        return 0;
    }
    switch (seedLen) {
        case 32: {
            switch (wrapCount) {
                case 1: return 16;
                case 2: return 8; 
                case 3: return 24;
                case 4: return 4;
                case 5: return 20;
                case 6: return 12;
                case 7: return 28;
                case 8: return 6;
                case 9: return 18;
                case 10: return 10;
                case 11: return 26;
                case 12: return 2;
                case 13: return 14;
                case 14: return 22;
                case 15: return 30;
                case 16: return 3;
                case 17: return 17;
                case 18: return 9;
                case 19: return 23;
                case 20: return 29;
                case 21: return 7;
                case 22: return 19;
                case 23: return 27;
                case 24: return 5;
                case 25: return 25;
                case 26: return 11;
                case 27: return 21;
                case 28: return 31;
                case 29: return 15;
                case 30: return 1;
                case 31: return 13;
            }
        }        
        case 31: {
            switch (wrapCount) {
                case 1: return 15;
                case 2: return 23; 
                case 3: return 8;
                case 4: return 19;
                case 5: return 4;
                case 6: return 27;
                case 7: return 11;
                case 8: return 17;
                case 9: return 6;
                case 10: return 25;
                case 11: return 2;
                case 12: return 29;
                case 13: return 13;
                case 14: return 21;
                case 15: return 9;
                case 16: return 24;
                case 17: return 3;
                case 18: return 18;
                case 19: return 10;
                case 20: return 26;
                case 21: return 5;
                case 22: return 30;
                case 23: return 14;
                case 24: return 7;
                case 25: return 22;
                case 26: return 1;
                case 27: return 16;
                case 28: return 28;
                case 29: return 12;
                case 30: return 20;
             }
        }
        case 25: {
            switch (wrapCount) {
                case 1: return 13;
                case 2: return 6; 
                case 3: return 19;
                case 4: return 3;
                case 5: return 16;
                case 6: return 22;
                case 7: return 9;
                case 8: return 11;
                case 9: return 1;
                case 10: return 14;
                case 11: return 7;
                case 12: return 20;
                case 13: return 4;
                case 14: return 17;
                case 15: return 23;
                case 16: return 2;
                case 17: return 15;
                case 18: return 5;
                case 19: return 21;
                case 20: return 8;
                case 21: return 24;
                case 22: return 10;
                case 23: return 18;
                case 24: return 12;
            }
        }

        case 24:{
            switch (wrapCount) {
                case 1: return 12;
                case 2: return 6; 
                case 3: return 18;
                case 4: return 3;
                case 5: return 15;
                case 6: return 21;
                case 7: return 9;
                case 8: return 1;
                case 9: return 13;
                case 10: return 19;
                case 11: return 7;
                case 12: return 16;
                case 13: return 4;
                case 14: return 22;
                case 15: return 10;
                case 16: return 2;
                case 17: return 14;
                case 18: return 20;
                case 19: return 5;
                case 20: return 17;
                case 21: return 8;
                case 22: return 23;
                case 23: return 11;
            }
        }

        case 23: {
            switch (wrapCount) {
                case 1: return 12;
                case 2: return 6;
                case 3: return 17;
                case 4: return 3;
                case 5: return 9;
                case 6: return 20;
                case 7: return 14;
                case 8: return 1;
                case 9: return 4;
                case 10: return 7;
                case 11: return 10;
                case 12: return 15;
                case 13: return 18;
                case 14: return 21;
                case 15: return 4;
                case 16: return 2;
                case 17: return 5;
                case 18: return 11;
                case 19: return 16;
                case 20: return 19;
                case 21: return 22;
                case 22: return 8;
            }
        }

        case 22: {
            switch (wrapCount) {
                case 1: return 11;
                case 2: return 6;
                case 3: return 16;
                case 4: return 3;
                case 5: return 9;
                case 6: return 14;
                case 7: return 19;
                case 8: return 2;
                case 9: return 7;
                case 10: return 12;
                case 11: return 17;
                case 12: return 20;
                case 13: return 4;
                case 14: return 1;
                case 15: return 10;
                case 16: return 13;
                case 17: return 15;
                case 18: return 18;
                case 19: return 21;
                case 20: return 5;
                case 21: return 8;
                default: _ASSERT(!"NOTREACHED");
            }
        }

        case 21: {
            switch (wrapCount) {
                case 1: return 11;
                case 2: return 6;
                case 3: return 16;
                case 4: return 3;
                case 5: return 9;
                case 6: return 13;
                case 7: return 17;
                case 8: return 18;
                case 9: return 2;
                case 10: return 5;
                case 11: return 8;
                case 12: return 15;
                case 13: return 20;
                case 14: return 1;
                case 15: return 4;
                case 16: return 7;
                case 17: return 10;
                case 18: return 12;
                case 19: return 14;
                case 20: return 19;
                default: _ASSERT(!"NOTREACHED");
            }
        }

        case 20: {
            switch (wrapCount) {
                case 1: return 10;
                case 2: return 5;
                case 3: return 15;
                case 4: return 2;
                case 5: return 7;
                case 6: return 12;
                case 7: return 17;
                case 8: return 3;
                case 9: return 9;
                case 10: return 11;
                case 11: return 13;
                case 12: return 19;
                case 13: return 1;
                case 14: return 4;
                case 15: return 6;
                case 16: return 8;
                case 17: return 14;
                case 18: return 18;
                case 19: return 16;
                default: _ASSERT(!"NOTREACHED");
            }
        }

        case 19: {
            switch (wrapCount) {
                case 1: return 10;
                case 2: return 4;
                case 3: return 14;
                case 4: return 2;
                case 5: return 6;
                case 6: return 8;
                case 7: return 12;
                case 8: return 16;
                case 9: return 18;
                case 10: return 1;
                case 11: return 3;
                case 12: return 5;
                case 13: return 7;
                case 14: return 9;
                case 15: return 11;
                case 16: return 13;
                case 17: return 15;
                case 18: return 17;
                default: _ASSERT(!"NOTREACHED");
            }
        }

        case 18: {
            switch (wrapCount) {
                case 1: return 9;
                case 2: return 4;
                case 3: return 13;
                case 4: return 2;
                case 5: return 6;
                case 6: return 11;
                case 7: return 15;
                case 8: return 1;
                case 9: return 3;
                case 10: return 5;
                case 11: return 7;
                case 12: return 8;
                case 13: return 10;
                case 14: return 12;
                case 15: return 14;
                case 16: return 16;
                case 17: return 17;
                default: _ASSERT(!"NOTREACHED");
            }
        }

        case 17: {
            switch (wrapCount) {
                case 1: return 8;
                case 2: return 4;
                case 3: return 12;
                case 4: return 2;
                case 5: return 6;
                case 6: return 10;
                case 7: return 14;
                case 8: return 1;
                case 9: return 3;
                case 10: return 5;
                case 11: return 7;
                case 12: return 9;
                case 13: return 11;
                case 14: return 13;
                case 15: return 15;
                case 16: return 16;
                default: _ASSERT(!"NOTREACHED");
            }
        }

        case 16: {
            switch (wrapCount) {
                case 1: return 8;
                case 2: return 4;
                case 3: return 12;
                case 4: return 2;
                case 5: return 6;
                case 6: return 10;
                case 7: return 14;
                case 8: return 1;
                case 9: return 3;
                case 10: return 5;
                case 11: return 7;
                case 12: return 9;
                case 13: return 11;
                case 14: return 13;
                case 15: return 15;
                default: _ASSERT(!"NOTREACHED");
            }
        } // inner switch
        default: WriteErrorMessage("SeedSequencer: Not set up to run with this seed size\n"); soft_exit(1);
    } // outer switch
    return 0;
} 
#endif  // 1