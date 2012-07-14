/*++

Module Name:

	HashTable.h

Abstract:

    Headers for large closed hash table that used to be general, but now just handles
    seeds for SNAP.

Authors:

    Bill Bolosky, March, 2011

Environment:


--*/

#pragma once

#include "Compat.h"


class SNAPHashTable {
    public:

        SNAPHashTable(
            unsigned      i_tableSize,
            bool          i_useBigAlloc = true);

        //
        // Load from file.
        //
        SNAPHashTable(char *loadFileName);

        ~SNAPHashTable();

        bool saveToFile(const char *saveFileName);

        //
        // Fails if either the table is full or key already exists.
        //
        bool Insert(unsigned key, const unsigned *data);

        size_t GetUsedElementCount() const {return usedElementCount;}
        size_t GetTableSize() const {return tableSize;}

        _uint64 GetHashTableMemorySize();

        unsigned GetKeySizeInBytes() const {return keySizeInBytes;}
        unsigned GetDataSizeInBytes() const {return dataSizeInBytes;}

        static inline unsigned hash(unsigned key) {
            //
            // Hash the key.  Use the hash finalizer from the 64 bit MurmurHash3, http://code.google.com/p/smhasher/wiki/MurmurHash3,
            // which is public domain code.
            //
    
            key ^= key >> 16; 
            key *= 0x85ebca6b; 
            key ^= key >> 13; 
            key *= 0xc2b2ae35; 
            key ^= key >> 16;
            return key;
        }

        inline unsigned *Lookup(unsigned key) const {
            unsigned tableIndex = hash(key) % tableSize;
            if (table[tableIndex].key == key && table[tableIndex].value1 != 0xffffffff) {
                return &(table[tableIndex].value1);
            } else {
                unsigned nProbes = 0;
                Entry* entry;
                unsigned value1;
                do {
                    nProbes++;
                    if (nProbes > tableSize + QUADRATIC_CHAINING_DEPTH) {
                        return NULL;
                    }
                    if (nProbes < QUADRATIC_CHAINING_DEPTH) {
                        tableIndex = (tableIndex + nProbes * nProbes) % tableSize;
                    } else {
                        tableIndex = (tableIndex + 1) % tableSize;
                    }
                    entry = &table[tableIndex];
                    value1 = entry->value1;
                } while (entry->key != key && value1 != 0xffffffff);

                extern _int64 nProbesInGetEntryForKey;
                nProbesInGetEntryForKey += nProbes;

                if (value1 == 0xffffffff) {
                    return NULL;
                } else {
                    return &(entry->value1);
                }
            }
        }

        //
        // A version of Lookup that works properly when the table is (nearly) full and the key being looked up isn't
        // there.  It's, as you might imagine, slower than Lookup.
        //
        unsigned *SlowLookup(unsigned key);

private:

        static const unsigned QUADRATIC_CHAINING_DEPTH = 5; // Chain quadratically for this long, then linerarly  Set to 0 for linear chaining

        struct Entry {
            unsigned        key;
            unsigned        value1;
            unsigned        value2;
        };

        // Free Entries have value1 == 0xffffffff


        Entry *table;
        size_t tableSize;
        static const unsigned keySizeInBytes = 4;
        static const unsigned dataSizeInBytes = 8;
        static const unsigned elementSize = keySizeInBytes + dataSizeInBytes;
        size_t usedElementCount;

        size_t virtualAllocSize;

        bool useBigAlloc;

    
        unsigned bytesToCheckForUnusedEntry;

        //
        // Returns either the entry for this key, or else the entry where the key would be
        // inserted if it's not in the table.
        //
        Entry* getEntryForKey(__in unsigned key) const;

        friend class SeedCountIterator;
};
