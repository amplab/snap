#pragma once

#include "Compat.h"
#include "BigAlloc.h"
#include "exit.h"
#include "Error.h"


//
// A hash function for numeric types.
//
template<typename T>
class NumericHash
{
public:
    inline _uint64 operator() (T value) {
        return (_uint64) (value * 131);
    }
};


//
// A fixed-size hash map that allows for efficient clearing and reuse through epochs
// and does not perform any memory allocation.
//
// This class only allows the capacity to be a power of 2.
//
// Use epoch + 1 as a tombstone for deleted values
//
// K must be a numeric type that supports shift, mask and xor operators.
//
template< typename K, typename V, typename Hash = NumericHash<K> >
class FixedSizeMap
{
public:
    FixedSizeMap(unsigned capacity_ = 16): entries(NULL), size(0) {
        reserve(capacity_);
    }

    ~FixedSizeMap() {
        delete[] entries;
    }

    void reserve(unsigned capacity) {
        if (!isPowerOf2(capacity)) {
            WriteErrorMessage("FixedSizeMap capacity must be a power of 2\n");
            soft_exit(1);
        }
        if (entries != NULL) {
            if (size > 0) {
                WriteErrorMessage("reserve() called on a non-empty FixedSizeMap\n");
                soft_exit(1);
            }
            delete[] entries;
        }
        this->capacity = capacity;
        this->mask = capacity - 1;
        entries = new Entry[capacity];
        for (unsigned i = 0; i < capacity; i++) {
            entries[i].epoch = 0;
        }
        epoch = 1;

        clearBloomFilter();
    }

    void clear() {
        size = 0;
        epoch += 2;
        if (epoch > 100000000) {
            // Reset the epoch of every bucket to 0 and the current epoch to 1
            for (unsigned i = 0; i < capacity; i++) {
                entries[i].epoch = 0;
            }
            epoch = 1;
        }

        clearBloomFilter();
    }

    void resize(unsigned size)
    {
        // Do something here to limit the size of the hash table to reduce cache missing.
        _ASSERT(size <= capacity);
    }

    static const unsigned MaxQuadraticProbes = 4;

    inline V get(K key) {
        unsigned pos = hash(key) & mask;
#if     0
        //
        // Prefetch the data.  If it hits in the Bloom Filter then we can overlap the cache fetch
        // with the Bloom Filter computation, making the latter essentially free.  If it's not in the
        // Bloom Filter, then this will bring the cache line in for the add that's doubtless coming
        // soon after.
        //
        _mm_prefetch((const char *)(&entries[pos]),_MM_HINT_T2);

        if (!checkBloomFilter(key)) {
            //
            // Not in the Bloom Filter means not in the cache.
            //
            return V();
        }
#endif  // 0

        unsigned i = 1;
        while (true) {
            if (entries[pos].epoch < epoch) {
                return V();
            } else if (entries[pos].key == key && entries[pos].epoch == epoch) {
                return entries[pos].value;
            } else {
                pos = (pos + (i <= MaxQuadraticProbes ? i : 1)) & mask;
                i++;
                if (i > capacity + MaxQuadraticProbes) {
                    return V();
                }
            }
        }
    }

    inline void put(K key, V value) {
        _ASSERT(size < capacity);
//        addToBloomFilter(key);

        unsigned pos = hash(key) & mask;
        unsigned i = 1;

        while (true) {
            if (entries[pos].epoch != epoch) {
                entries[pos].key = key;
                entries[pos].value = value;
                entries[pos].epoch = epoch;
                size++;
                return;
            } else if (entries[pos].key == key) {
                entries[pos].value = value;
                return;
            } else {
                pos = (pos + (i <= MaxQuadraticProbes ? i : 1)) & mask;
                i++;
                _ASSERT(i <= capacity + MaxQuadraticProbes); // todo: overlow condition?
            }
        }
    }

    inline void erase(K key) {
        _ASSERT(size <= capacity);
        unsigned pos = hash(key) & mask;
        unsigned i = 1;
        while (true) {
            if (entries[pos].epoch < epoch) {
                return;
            } else if (entries[pos].key == key && entries[pos].epoch == epoch) {
                entries[pos].epoch = epoch + 1; // mark with tombstone
                size--;
                return;
            } else {
                pos = (pos + (i <= MaxQuadraticProbes ? i : 1)) & mask;
                i++;
                _ASSERT(i <= capacity + MaxQuadraticProbes); // todo: overlow condition?
            }
        }
    }

    inline int getSize() { return size; }

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

    typedef void* iterator;

    iterator begin()
    {
        return next(&entries[-1]);
    }

    iterator next(iterator i)
    {
        Entry* final = &entries[capacity];
        Entry* x = (Entry*) i;
        if (x < final) {
            do {
                x++;
            } while (x < final && x->epoch != epoch);
        }
        return x;
    }

    iterator end()
    {
        return &entries[capacity];
    }

    K key(iterator i)
    {
        return ((Entry*)i)->key;
    }

    V& value(iterator i)
    {
        return ((Entry*)i)->value;
    }

private:
    //
    // To avoid cache misses on failed lookups, we have a cheezy Bloom filter.  It's fixed at 512 bits (which is 64 bytes, typically
    // a cache line), and two features.
    //
    static const unsigned bloomFilterFeatureSizeInBits = 9; // Must be >=3.  Using 9 results in 64 bytes of Bloom Filter, which is cache-line sized (though not necessarily aligned)
    static const unsigned bloomFilterSizeInChar = (1 << (bloomFilterFeatureSizeInBits - 3));
    static const _uint64 bloomFilterFeatureMask = (1 << bloomFilterFeatureSizeInBits) - 1;

    unsigned char bloomFilter[bloomFilterSizeInChar];

    static inline void getBloomFilterFeatures(K key, unsigned *feature0Word, unsigned *feature0Bit, unsigned *feature1Word, unsigned *feature1Bit)
    {
        //
        // We know the bloom filter is 2^bloomFilterFeatureSizeInBits bits wide.  Use alternating bloomFilterFeatureSizeInBits bit chunks from the key to build up each of the features.
        //
        _uint64 feature[2] = {0, 0};

        for (int i = 0; i < sizeof(K) * 8; i += bloomFilterFeatureSizeInBits * 2) {
            feature[0] ^= ((key >> i) & bloomFilterFeatureMask);
            feature[1] ^= ((key >> (i+bloomFilterFeatureSizeInBits)) & bloomFilterFeatureMask);
        }

        *feature0Word = feature[0] / 8;
        *feature0Bit = feature[0] % 8;
        *feature1Word = feature[1] / 8;
        *feature1Bit = feature[1] % 8;

        _ASSERT(*feature0Word < bloomFilterSizeInChar && *feature1Word < bloomFilterSizeInChar);
    }

    //
    // false means that this entry is NOT in the cache.  true means we can't say for sure.
    //
    inline bool checkBloomFilter(K key)
    {
        unsigned feature0Word, feature0Bit, feature1Word, feature1Bit;

        getBloomFilterFeatures(key, &feature0Word, &feature0Bit, &feature1Word, &feature1Bit);

        return (bloomFilter[feature0Word] & (1 << feature0Bit)) && (bloomFilter[feature1Word] & (1 << feature1Bit));
    }

    inline void addToBloomFilter(K key)
    {
        unsigned feature0Word, feature0Bit, feature1Word, feature1Bit;

        getBloomFilterFeatures(key, &feature0Word, &feature0Bit, &feature1Word, &feature1Bit);

        bloomFilter[feature0Word] |= 1 << feature0Bit;
        bloomFilter[feature1Word] |= 1 << feature1Bit;
    }

    void clearBloomFilter()
    {
       memset(bloomFilter, 0, bloomFilterSizeInChar * sizeof(bloomFilter[0]));
    }

    struct Entry {
        K key;
        V value;
        int epoch;

        void *operator new[](size_t size) {return BigAlloc(size);}
        void operator delete[](void *ptr) {BigDealloc(ptr);}
    };

    Entry *entries;
    unsigned capacity;
    unsigned size;
    int mask;
    int epoch;
    Hash hash;



    bool isPowerOf2(int n) {
        while (n > 0) {
            if (n == 1) {
                return true;
            } else if (n % 2 == 1) {
                return false;
            } else {
                n /= 2;
            }
        }
        return false;
    }
};

