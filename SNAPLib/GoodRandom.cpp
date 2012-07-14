/*++


Module Name:

	GoodRandom.cpp

Abstract:

    Cryptographically secure random number generator that's fast because it batches getting its random
    source from CryptGenRand().  It's careful not to bias results when the request range isn't a power
    of two.

    In Linux, uses the Mersenne Twister to get the raw numbers and mostly just applies the anti-biasing
    related to random number ranges that aren't powers of two.

Authors:

    Bill Bolosky, April, 2011

Environment:

--*/
#include "stdafx.h"
#include "GoodRandom.h"
#include "mt64.h"

unsigned MinBytesToStore(_uint64 maxValue)
{
    if (0 == maxValue) {
        return 0;
    } else if (maxValue < 0x100) {
        return 1;
    } else if (maxValue < 0x10000) {
        return 2;
    } else if (maxValue < 0x1000000) {
        return 3;
    } else if (maxValue < 0x100000000) {
        return 4;
    } else if (maxValue < 0x10000000000) {
        return 5;
    } else if (maxValue < 0x1000000000000) {
        return 6;
    } else if (maxValue < 0x100000000000000) {
        return 7;
    } else {
        // Wow
        return 8;
    }
}

#ifdef _MSC_VER
const unsigned GFRandomBufferSize = 10 * 1024 * 1024;

struct GFRandomState {
    char               *buffer;
    unsigned            bufferUsed;
};

__declspec(thread) GFRandomState *g_randomState = NULL;

__declspec(thread) HCRYPTPROV g_hRandomCryptProv;
#endif

_uint64 GoodFastRandom(_uint64 maxValue)
{
    if (0 == maxValue) {
        return 0;
    }
#ifdef  _MSC_VER
    if (maxValue > 0xffffffffffffff && maxValue != 0xffffffffffffffff) {
        fprintf(stderr,"GoodFastRandom: writeme\n");
        exit(1);
        // This case requires different math in the roundoff error check below, because as written
        // it would use 1 << 64, which isn't representable.  I'm too lazy to write it now.
    }

    if (NULL == g_randomState) {
        g_randomState= new GFRandomState;
        g_randomState->buffer = new char[GFRandomBufferSize];
        g_randomState->bufferUsed = GFRandomBufferSize; // Forces us to get new random data the first time through.

        BOOL worked = CryptAcquireContext(
                                             &g_hRandomCryptProv,
                                             "BillKeyContainer",
                                             NULL,             // default provider
                                             PROV_RSA_FULL,
                                             CRYPT_MACHINE_KEYSET);

        if (!worked && NTE_BAD_KEYSET == GetLastError()) {
            worked = CryptAcquireContext(
                                                 &g_hRandomCryptProv,
                                                 "BillKeyContainer",
                                                 NULL,             // default provider
                                                 PROV_RSA_FULL,
                                                 CRYPT_MACHINE_KEYSET | CRYPT_NEWKEYSET);
                
        }

        if (!worked) {
            fprintf(stderr,"Unable to get crypt provider, %d\n",GetLastError());
            exit(1);
        }
    }
#endif   // _MSC_VER

    unsigned bytesToGet = MinBytesToStore(maxValue);

    for (;;) {
        _uint64 rawValue = 0;

#ifdef  _MSC_VER
        if (g_randomState->bufferUsed + bytesToGet > GFRandomBufferSize) {
            if (!CryptGenRandom(g_hRandomCryptProv,GFRandomBufferSize,(PBYTE)g_randomState->buffer)) {
                fprintf(stderr,"CryptGenRandom failed, %d\b\n",GetLastError());
                exit(1);
            }
            g_randomState->bufferUsed = 0;
        }

        memcpy(&rawValue,g_randomState->buffer + g_randomState->bufferUsed,bytesToGet);
        g_randomState->bufferUsed += bytesToGet;
#else   // _MSC_VER
        #ifdef RELEASE
        #if RAND_MAX % 0x100 != 0xff
        fprintf(stderr,"Jesse was too lazy to correct random bias on your platform.\n");
        exit(1);
        #endif
        for (unsigned b = 0; b < bytesToGet; ++b) {
            // RAND_MAX could in theory be less than 0xffff.
            _uint64 randomByte = rand() & 0xff;
            rawValue |= randomByte << (b * 8);
        }
        #else
        rawValue = genrand64_int64() >> (8 - bytesToGet) * 8;
        #endif // RELEASE
#endif  // _MSC_VER

        //
        // Be careful here not to bias the result.  If maxValue + 1 doesn't go evenly into 256^bytesToGet
        // then there would be a bias toward the lower values, since they get one extra representation.
        // So, if the value is in the last part, throw it away and try again.
        //
        if (0xffffffffffffffff == maxValue) {
            //
            // Special case for full range, because it wouldn't work with the code below.
            //
            return rawValue;
        }

        _uint64 maxRawValuePlusOne = ((_uint64)1) << (bytesToGet * 8);
        if (rawValue < maxRawValuePlusOne - maxRawValuePlusOne % (maxValue+1)) {
            return rawValue % (maxValue+1);
        }
    }

    /*NOTREACHED*/
}

