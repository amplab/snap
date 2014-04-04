/*++

Module Name:

    Util.cpp

Abstract:

    Generic support routines that don't seem to belong elsewhere.

Authors:

    Bill Bolosky, March, 2013

Environment:

    User mode service.

Revision History:

    Factored from other places

--*/

#include "stdafx.h"
#include "Util.h"


int FirstPowerOf2GreaterThanOrEqualTo(int value)
{
    int highestBitSet;
    for (highestBitSet = 0; highestBitSet <= 30; highestBitSet++) { // Only go to 31, since this is signed
        if (!(value & ~((1 << highestBitSet) - 1))) {
            highestBitSet -= 1;
            break;
        }
    }
    if (1 << highestBitSet == value) return value;
    return 1 << (highestBitSet + 1);
}

int cheezyLogBase2(_int64 value)
{
    int retVal = 0;
    value /= 2; // Since 2^0 = 1; we'll also define cheezyLogBase2(x) = 0 where x<= 0.
    while (value > 0) {
        retVal++;
        value /= 2;
    }
    return retVal;
}

    void
util::memrevcpy(
    void* dst,
    const void* src,
    size_t bytes)
{
    size_t dwords = bytes >> 3;
    _uint64* p = (_uint64*) dst;
    const _uint64* q = (const _uint64*) ((const char*)src + bytes - 8);
    for (size_t i = 0; i < dwords; i++) {
        *p++ = ByteSwapUI64(*q--);
    }
    int left = (int) (bytes & 7);
    if (left > 0) {
        char* p2 = (char*) p;
        const char* q2 = (left - 1) + (const char*) src;
        for (int i = 0; i < left; i++) {
            *p2++ = *q2--;
        }
    }
}

NWaiter::NWaiter(size_t n)
{
	_signalsRequired = n;
	_signalsReceived = 0;
	InitializeExclusiveLock(&_lock);
	CreateEventObject(&_waiter);
}

NWaiter::~NWaiter()
{
	DestroyExclusiveLock(&_lock);
	DestroyEventObject(&_waiter);
}


void NWaiter::wait()
{
	while (true) {
		bool done;

		AcquireExclusiveLock(&_lock);
		done = (_signalsReceived >= _signalsRequired);
		ReleaseExclusiveLock(&_lock);

		if (done)
			return;
		else {
			WaitForEvent(&_waiter);
		}
	}
}

void NWaiter::signal()
{
	AcquireExclusiveLock(&_lock);
	_signalsReceived += 1;
	ReleaseExclusiveLock(&_lock);

	AllowEventWaitersToProceed(&_waiter);
}
