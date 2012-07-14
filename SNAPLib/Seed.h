/*++

Module Name:

    Seed.h

Abstract:

    Headers for code to handle seeds in the SNAP sequencer.

Authors:

    Bill Bolosky, August, 2011

Environment:
`
    User mode service.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once

#include "Compat.h"
#include "Tables.h"


struct Seed {
    //
    // We exclude seeds with "N" in them.  This checks for that (or other garbage).
    //
    static bool DoesTextRepresentASeed(const char *textBases, unsigned seedLen);

    inline Seed(const char *textBases, unsigned seedLen)
    {

        bases = 0;
        reverseComplement = 0;

        for (unsigned i = 0; i < seedLen; i++) {
            _uint64 encodedBase = BASE_VALUE[textBases[i]];
            _ASSERT(255 != encodedBase);

            bases |= encodedBase << ((seedLen - i - 1) * 2);
            reverseComplement |= (encodedBase ^ 0x3) << (i * 2);
        }
    }

    inline Seed() {}

    inline Seed(_int64 i_bases, _int64 i_reverseComplement)
        : bases(i_bases), reverseComplement(i_reverseComplement)
    {
    }

    inline unsigned getLowBases() const {   // Returns the lowest 16 bases as an unsigned
        return (unsigned) (bases & 0xffffffff);
    }

    inline unsigned getHighBases() const {   // Returns any bases > 16 as an unsigned.  If seedLen <= 16, returns 0.
        return (unsigned) ((bases >> 32) & 0xffffffff);
    }

    inline _int64 getBases() const {
        return bases;
    }

    inline Seed operator~() const
    //
    // Compute the reverse complement of this.  We just copy it from when our constructor ran.
    //
    {
        Seed rc;

        rc.bases = reverseComplement;
        rc.reverseComplement = bases;

        return rc;
    }

    inline bool isBiggerThanItsReverseComplement() const {
        return bases > reverseComplement;
    }

    inline bool isOwnReverseComplement() const {
        return bases == reverseComplement;
    }

    inline bool operator>(Seed &peer) const {
        return bases > peer.bases;
    }

    inline bool operator>=(Seed &peer) const {
        return bases >= peer.bases;
    }

    inline bool operator<(Seed &peer) const {
        return bases < peer.bases;
    }

    inline bool operator<=(Seed &peer) const {
        return bases <= peer.bases;
    }

    inline bool operator==(Seed &peer) const {
        return bases == peer.bases;
    }

    inline bool operator!=(Seed &peer) const {
        return bases != peer.bases;
    }

    inline void setBase(int i, int seedLen, int value) {
        int shift = (seedLen - i - 1) * 2;
        _int64 mask = (_int64) 3 << shift;
        bases = (bases & ~mask) | (((_int64) value << shift) & mask);
        int shift2 = i * 2;
        _int64 mask2 = (_int64) 3 << shift2;
        reverseComplement = (reverseComplement & ~mask2) | (((_int64) (value ^ 3) << shift2) & mask2);
    }

    inline int getBase(int i, int seedLen) {
        return (int) (bases >> ((seedLen - i - 1) * 2)) & 3;
    }

    inline void shiftIn(int b, int seedLen) {
        int shift = (seedLen - 1) * 2;
        bases = ((bases << 2) & ~((_uint64) 3 << (shift + 2))) | (b & 3);
        reverseComplement = (reverseComplement >> 2) | (((_uint64) (b ^ 3) & 3) << shift);
    }

    inline void toString(char* o_bases, int seedLength) {
        for (int i = (seedLength - 1) * 2; i >= 0; i -= 2) {
            *o_bases++ = VALUE_BASE[(bases >> i) & 3];
        }
    }

    inline static Seed fromBases(_int64 bases, int seedLength) {
        _int64 rc = 0;
        _int64 b = bases;
        for (int i = 0; i < seedLength; i++) {
            rc = (rc << 2) | ((b & 3) ^ 3);
            b = b >> 2;
        }
        return Seed(bases, rc);
    }

private:

    _int64   bases;

    //
    // Since we pretty much always compute the reverse complement of a seed, we just keep it
    // here.  That way we only execute the loop once: when the constructor runs.
    //
    _int64   reverseComplement;
};
