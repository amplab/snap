#pragma once

#include "Compat.h"
#include "BigAlloc.h"
#include "FixedSizeMap.h"


// 
// A fixed-capacity hash set that allows for efficient clearing and reuse through epochs
// and does not perform any memory allocation.
//
// This class only allows the capacity to be a power of 2.
//
template< typename K, typename Hash = NumericHash<K> >
class FixedSizeSet
{
public:
    FixedSizeSet(int capacity_ = 16): entries(NULL), size(0) {
        reserve(capacity_);
    }

    ~FixedSizeSet() {
        delete[] entries;
    }
    
    void reserve(int capacity) {
        if (!isPowerOf2(capacity)) {
            fprintf(stderr, "FixedSizeSet capacity must be a power of 2\n");
            exit(1);
        }
        if (entries != NULL) {
            if (size > 0) {
                fprintf(stderr, "reserve() called on a non-empty FixedSizeSet\n");
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
    
    inline bool contains(K key) {
        unsigned pos = hash(key) & mask;
        int i = 1;
        while (true) {
            if (entries[pos].epoch != epoch) {
                return false;
            } else if (entries[pos].key == key) {
                return true;
            } else {
                pos = (pos + i) & mask;
                i++;
            }
        }
    }

    inline void add(K key) {
        _ASSERT(size < capacity);
        unsigned pos = hash(key) & mask;
        int i = 1;
        while (true) {
            if (entries[pos].epoch != epoch) {
                entries[pos].key = key;
                entries[pos].epoch = epoch;
                size++;
                return;
            } else if (entries[pos].key == key) {
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
