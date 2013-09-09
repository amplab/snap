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
#include "Genome.h"


class SNAPHashTable {
    public:

        SNAPHashTable(
            unsigned      i_tableSize,
            unsigned      i_keySizeInBytes);

        //
        // Load from file.
        //
        static SNAPHashTable *loadFromFile(char *loadFileName);

        static SNAPHashTable *loadFromFile(FILE *loadFile);

        ~SNAPHashTable();

        bool saveToFile(const char *saveFileName);

        bool saveToFile(FILE *saveFile);

        //
        // Fails if either the table is full or key already exists.
        //
        bool Insert(_uint64 key, const unsigned *data);

        size_t GetUsedElementCount() const {return usedElementCount;}
        size_t GetTableSize() const {return tableSize;}

        _uint64 GetHashTableMemorySize();

        unsigned GetKeySizeInBytes() const {return keySizeInBytes;}
        unsigned GetDataSizeInBytes() const {return dataSizeInBytes;}

        static inline _uint64 hash(_uint64 key) {
            //
            // Hash the key.  Use the hash finalizer from the 64 bit MurmurHash3, http://code.google.com/p/smhasher/wiki/MurmurHash3,
            // which is public domain code.
            //
           
            key ^= key >> 33;
            key *= 0xff51afd7ed558ccd;
            key ^= key >> 33;
            key *= 0xc4ceb9fe1a85ec53;
            key ^= key >> 33;
            
            return key;
        }

        inline unsigned *Lookup(_uint64 key) const {
            _ASSERT(keySizeInBytes == 8 || (key & ~((((_uint64)1) << (keySizeInBytes * 8)) - 1)) == 0);    // High bits of the key aren't set.
            _uint64 tableIndex = hash(key) % tableSize;
            Entry *entry = getEntry(tableIndex);
            if (isKeyEqual(entry, key) && entry->value1 != InvalidGenomeLocation) {
                return &(entry->value1);
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
                    entry = getEntry(tableIndex);
                    value1 = entry->value1;
                } while (!isKeyEqual(entry, key) && value1 != InvalidGenomeLocation);

                extern _int64 nProbesInGetEntryForKey;
                nProbesInGetEntryForKey += nProbes;

                if (value1 == InvalidGenomeLocation) {
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
        unsigned *SlowLookup(_uint64 key);

private:

        SNAPHashTable() {}

        static const unsigned QUADRATIC_CHAINING_DEPTH = 5; // Chain quadratically for this long, then linerarly  Set to 0 for linear chaining

        struct Entry {
            unsigned        value1;
            unsigned        value2;
            unsigned char   key[1]; // Actual size of key determined by keySizeInBytes
        };

        // Free Entries have value1 == InvalidGenomeLocation
        
        inline Entry *getEntry(_uint64 whichEntry) const {
            return (Entry *) ((char *)Table + elementSize * whichEntry);
        }


        inline bool isKeyEqual(const Entry *entry, _uint64 key) const 
        {
            return !memcmp(entry->key, &key, keySizeInBytes);
        }

        inline void clearKey(Entry *entry)
        {
            memset(entry->key, 0 , keySizeInBytes);
        }

        inline void setKey(Entry *entry, _uint64 key)
        {
            memcpy(entry->key, &key, keySizeInBytes);
        }
 
        Entry *Table;
        size_t tableSize;
        unsigned keySizeInBytes;
        static const unsigned dataSizeInBytes;
        unsigned elementSize;
        size_t usedElementCount;

        size_t virtualAllocSize;

        //
        // Returns either the entry for this key, or else the entry where the key would be
        // inserted if it's not in the table.
        //
        Entry* getEntryForKey(__in _uint64 key) const;

        friend class SeedCountIterator;

        static const unsigned magic;
};
