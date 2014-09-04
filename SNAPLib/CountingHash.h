/*++

Module Name:

    CountingHash.h

Abstract:

    Approximate counting hash table

Authors:

    Ravi Pandya, March 2014

Environment:
    User mode service.

    Uses interlocked operations so it is thread-safe

Revision History:

--*/

#pragma once 
#include "stdafx.h"
#include "Compat.h"
#include "Tables.h"
#include "Util.h"
#include "Error.h"

// approximate hash table 

template < typename Key >
class SaltedNumericHash
{
public:
    inline _uint64 operator() (const Key& key, _uint64 salt) const
    { return util::hash64(key ^ salt); }
};

// probabilistic hash table for storing incremental counts
// key must have _uint64 hash64(_uint64 salt) that gives independent values for each salt value
// bits per value must be factor of 64 (i.e. 1, 2, 4, 8, 16, 32, 64)
// when Bits=1 this is just a Bloom filter
template < typename Key >
class CountingHashBase
{
public:

    int BitsForMaxValue(_uint64 value)
    {
        if (value <= 1) {
            return 1;
        }
        unsigned long n;
        CountLeadingZeroes(value, n);
        return 64 - n;
    }

    // calculation from http://en.wikipedia.org/wiki/Bloom_filter
    CountingHashBase(_uint64 i_keyCount, _uint64 i_maxValue, double pError = 1e-4, bool i_big = false) :
		maxKeyCount(i_keyCount),
		keyCount(0),
		capacity((_uint64) (0.5 - (i_keyCount * log(pError) / (log(2.0) * log(2.0))))),
		probes((int) (0.5 - log(pError) / log(2.0))),
        bits(BitsForMaxValue(i_maxValue)),
        maxValue(i_maxValue),
        packing(64 / BitsForMaxValue(i_maxValue)),
        mask((((_uint64) 1) << BitsForMaxValue(i_maxValue)) - 1),
        big(i_big)
    {
        if (big) {
            data = (_uint64*) BigAlloc(8 * (1 + capacity / packing));
        } else {
            data = new _uint64[1 + capacity / packing];
            memset(data, 0, 8 * (1 + capacity / packing));
        }
    }

    ~CountingHashBase()
    {
        if (big) {
            BigDealloc(data);
        } else {
            delete[] data;
        }
    }

    // get count for a key
    _uint64 get(const Key& key)
    {
        _uint64 result = maxValue;;
        for (int i = 0; i < probes; i++) {
            _uint64 index = hash64(key, i) % capacity;
            _uint64 dw = data[index / packing];
            _uint64 val = (dw >> (bits * (index % packing))) & mask;
            if (val < result) {
                result = val;
            }
        }
        return result;
    }

    // increment count for a key, return PREVIOUS count
    _uint64 increment(const Key& key)
    {
        _uint64 result = maxValue;;
		keyCount++;  //increment to keep track of how many keys are in here. 
        for (int i = 0; i < probes; i++) {
            _uint64 index = hash64(key, i) % capacity;
            _uint64 slot = index / packing;
            int shift = bits * (index % packing);
            while (true) {
                _uint64 dw = data[slot];
                _uint64 val = (dw >> shift) & mask;
                if (val < maxValue) {
                    _uint64 dw2 = (dw & ~(mask << shift)) | ((val + 1) << shift);
                    if (InterlockedCompareExchange64AndReturnOldValue(&data[slot], dw2, dw) != dw) {
                        continue; // try again
                    }
                    if (val < result) {
                        result = val;
                    }
                }
                break;
            }
        }
        return result;
    }

	_uint64 getBitsSet(){
		_uint64 bitsSet = 0; 
		for (_uint64 slot = 0; slot < capacity/packing; slot++){
			_uint64 dw = data[slot]; 
			for (int i = 0; i< packing; i++){
				_uint64 val = (dw >> i) & mask;
				if (val >0){
					bitsSet++;
				}
			}
		}
		return bitsSet; 
	}
    
	bool isTooFull(){
		return keyCount > maxKeyCount; }
	
	_uint64 getKeyCount(){
		return keyCount; }

	double calculatePerror(){
		return pow((1 - exp(-(double)probes*keyCount / capacity)), probes); }

    virtual _uint64 hash64(const Key& key, _uint64 salt) = 0;

private:
    const _uint64 capacity; // #slots in hashtable
	const _uint64 maxKeyCount; // #maximum number of keys you can put in the table and maintain your fp rate (p-value)
	const _uint64 maxValue; // maximum possible value, will be < 2^bits
    const int bits; // bits per value
    const _uint64 mask; // mask to AND out value (2^bits-1)
    const int packing; // #values packed per DWORD
    const int probes; // #hash probes for each value
    const bool big; // BigAlloc vs. operator new for data
    _uint64* data;
	_uint64  keyCount;  //keep track of the number of keys in the table.  If it's greater than max keycount, the bloom filter is useless (always returns positive) 
};

// counter for integer keys
class CountingHash : public CountingHashBase<_uint64>
{
public:
    CountingHash(_uint64 i_keyCount, _uint64 i_maxValue, double pError = 1e-4, bool i_big = false)
        : CountingHashBase(i_keyCount, i_maxValue, pError, i_big)
    {}

    ~CountingHash() {}

    _uint64 hash64(const _uint64& key, _uint64 salt)
    { return util::hash64(key ^ salt); }
};

