#include "stdafx.h"
#include "Tables.h"


static const Tables tables;

const char *COMPLEMENT = tables.getComplement();
const char *IS_N = tables.getIsN();
const int  *BASE_VALUE = tables.getBaseValue();
const int  *BASE_VALUE_NO_N = tables.getBaseValueNoN();
const char *VALUE_BASE = tables.getValueBase();
const unsigned char *VALUE4_RC = tables.getValue4RC();
const char *PACKED_BASE_VALUE = tables.getPackedBaseValue();
const char *PACKED_QUALITY_MASK = tables.getPackedQualityMask();
const char *PACKED_VALUE_BASE = tables.getPackedValueBase();
const unsigned *IS_LOWER_CASE_OR_DOT = tables.getIsLowerCaseOrDot();
const char *TO_UPPER_CASE_DOT_TO_N = tables.getToUpperCaseDotToN();
const char *PACKED_VALUE_BASE_RC = tables.getPackedValueBaseRC();
const char *CIGAR_QUAL_TO_SAM = tables.getCigarQualToSam();

Tables::Tables()
{
    memset(complement, 0, sizeof(complement));
    memset(isN, 0, sizeof(isN));

    complement['A'] = 'T';
    complement['C'] = 'G';
    complement['G'] = 'C';
    complement['T'] = 'A';
    complement['N'] = 'N';
    complement['n'] = 'n';

    isN['N'] = 1;
    isN['n'] = 1;

    // Base values chosen so that complements are bitwise opposites.
    for (unsigned i = 0; i < 256; i++) {
        baseValue[i] = 4;// Everything's an N unless it's not
    }
    baseValue['A'] = 0;
    baseValue['G'] = 1;
    baseValue['C'] = 2;
    baseValue['T'] = 3;

    // inverse of BASE_VALUE
    valueBase[0] = 'A';
    valueBase[1] = 'G';
    valueBase[2] = 'C';
    valueBase[3] = 'T';
    valueBase[4] = 'N';

    // Version that maps N's value to 0 instead of 4
    memset(baseValueNoN, 0, sizeof(baseValueNoN));
    baseValueNoN['A'] = 0;
    baseValueNoN['G'] = 1;
    baseValueNoN['C'] = 2;
    baseValueNoN['T'] = 3;

    // reverse complement of a byte of 4x2-bit values
    for (int i = 0; i < 256; i++) {
        value4RC[i] = 0xff ^ (((i & 0x03) << 6) | ((i & 0x0c) << 2) | ((i & 0x30) >> 2) | ((i & 0xc0) >> 6));
    }

    // packed base tables
    for (int i = 0; i < 256; i++) {
        packedValueBase[i] = i < 4 ? 'N' : "AGCT"[i >> 6];
        packedValueBaseRC[i] = i < 4 ? 'N' : "TCGA"[i >> 6];
    }

    memset(packedBaseValue, 0, sizeof(packedBaseValue));
    packedBaseValue['A'] = packedBaseValue['a'] = 0x00;
    packedBaseValue['G'] = packedBaseValue['g'] = 0x40;
    packedBaseValue['C'] = packedBaseValue['c'] = (char) 0x80;
    packedBaseValue['T'] = packedBaseValue['t'] = (char) 0xc0;
    
    memset(packedQualityMask, 0, 4);
    memset(packedQualityMask + 4, 0x3f, sizeof(packedQualityMask) - 4);

    for (unsigned i = 0; i < 256; i++) {
        isLowerCaseOrDot[i] = 0;
        toUpperCaseDotToN[i] = i;
    }
    for (unsigned i = 0x61; i <= 0x7a; i++) {
        isLowerCaseOrDot[i] = 1;
        toUpperCaseDotToN[i] = i - 0x20;
    }

    isLowerCaseOrDot['.'] = 1;
    toUpperCaseDotToN['.'] = 'N';

    for (unsigned i = 0; i < 256; i++) {
        cigarQualToSam[i] = i > ('~' - '!') ? '!' : '!' + i;
    }
}
