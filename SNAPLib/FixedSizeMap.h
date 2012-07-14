#pragma once

#include "Compat.h"
#include "BigAlloc.h"


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
template< typename K, typename V, typename Hash = NumericHash<K> >
class FixedSizeMap
{
public:
    FixedSizeMap(int capacity_ = 16): entries(NULL), size(0) {
        reserve(capacity_);
    }

    ~FixedSizeMap() {
        delete[] entries;
    }
    
    void reserve(int capacity) {
        if (!isPowerOf2(capacity)) {
            fprintf(stderr, "FixedSizeMap capacity must be a power of 2\n");
            exit(1);
        }
        if (entries != NULL) {
            if (size > 0) {
                fprintf(stderr, "reserve() called on a non-empty FixedSizeMap\n");
                exit(1);
            }
            delete[] entries;
        }
        this->capacity = capacity;
        this->mask = capacity - 1;
        entries = new Entry[capacity];
        for (int i = 0; i < capacity; i++) {
            entries[i].epoch = 0;
        }
        epoch = 1;
    }
    
    void clear() {
        size = 0;
        epoch++;
        if (epoch > 100000000) {
            // Reset the epoch of every bucket to 0 and the current epoch to 1
            for (int i = 0; i < capacity; i++) {
                entries[i].epoch = 0;
            }
            epoch = 1;
        }
    }
    
    inline V get(K key) {
        unsigned pos = hash(key) & mask;
        int i = 1;
        while (true) {
            if (entries[pos].epoch != epoch) {
                return V();
            } else if (entries[pos].key == key) {
                return entries[pos].value;
            } else {
                pos = (pos + i) & mask;
                i++;
            }
        }
    }

    inline void put(K key, V value) {
        _ASSERT(size < capacity);
        unsigned pos = hash(key) & mask;
        int i = 1;
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
                pos = (pos + i) & mask;
                i++;
            }
        }
    }
    
    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

private:
    struct Entry {
        K key;
        V value;
        int epoch;

        void *operator new[](size_t size) {return BigAlloc(size);}
        void operator delete[](void *ptr) {BigDealloc(ptr);}
    };
    
    Entry *entries;
    int capacity;
    int size;
    int maxSize;
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
