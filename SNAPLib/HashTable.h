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
#include "GenericFile_Blob.h"
#include "Genome.h"


class SNAPHashTable {
    public:

        typedef _uint64 ValueType;  // Values can be smaller than this, but they're expanded in the interface
        typedef _uint64 KeyType;    // Likewise for keys.

        SNAPHashTable(
            _int64		i_tableSize,
            unsigned    i_keySizeInBytes,
            unsigned    i_valueSizeInBytes,
            unsigned    i_valueCount,
            _uint64		i_invalidValueValue);

        //
        // Load from file.
        //
        static SNAPHashTable *loadFromBlob(GenericFile_Blob *loadFile);
		static SNAPHashTable *loadFromGenericFile(GenericFile *loadFile);

        ~SNAPHashTable();

        bool saveToFile(const char *saveFileName, size_t *bytesWritten);

        bool saveToFile(FILE *saveFile, size_t *bytesWritten);

        //
        // Fails if either the table is full or key already exists.  Inserts ALL
        // values for a key.
        //
        bool Insert(KeyType key, ValueType *data);

        size_t GetUsedElementCount() const {return usedElementCount;}
        size_t GetTableSize() const {return tableSize;}

        unsigned GetKeySizeInBytes() const {return keySizeInBytes;}
        unsigned GetValueSizeInBytes() const {return valueSizeInBytes;}
        unsigned GetValueCount() const {return valueCount;}

		void *getEntryValues(_uint64 whichEntry) 
		{
			_ASSERT(whichEntry < GetTableSize());
			return getEntry(whichEntry);
		}

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

        inline ValueType *GetFirstValueForKey(KeyType key) const {
            _ASSERT(keySizeInBytes == 8 || (key & ~((((_uint64)1) << (keySizeInBytes * 8)) - 1)) == 0);    // High bits of the key aren't set.
            _uint64 tableIndex = hash(key) % tableSize;
            void *entry = getEntry(tableIndex);
            if (isKeyEqual(entry, key) && !doesEntryHaveInvalidValue(entry)) {
                return (ValueType *)entry;
            } else {
                unsigned nProbes = 0;
                void* entry;
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
                } while (!isKeyEqual(entry, key) && !doesEntryHaveInvalidValue(entry));

                extern _int64 nProbesInGetEntryForKey;
                nProbesInGetEntryForKey += nProbes;

                if (doesEntryHaveInvalidValue(entry)) {
                    return NULL;
                } else {
                    return (ValueType *)entry;
                }
            }
        }


        inline bool Lookup(KeyType key, unsigned nValuesToFill, ValueType *values) const {
            _ASSERT(nValuesToFill <= valueCount);
            char *entry = (char *)GetFirstValueForKey(key);
            if (NULL == entry) {
                return false;
            }

            for (unsigned i = 0; i < nValuesToFill; i++) {
                values[i] = 0;
                memcpy(values + i, entry + i * valueSizeInBytes, valueSizeInBytes);
            }

            return true;
        }
        //
        // A version of Lookup that works properly when the table is (nearly) full and the key being looked up isn't
        // there.  It's, as you might imagine, slower than Lookup.
        //
        ValueType *SlowLookup(KeyType key);

private:

        SNAPHashTable() {}

        static const unsigned QUADRATIC_CHAINING_DEPTH = 5; // Chain quadratically for this long, then linerarly  Set to 0 for linear chaining
		static SNAPHashTable *loadCommon(GenericFile *loadFile);

        //
        // A hash table entry consists of a set of valueCount values, each of valueSizeInBytes bytes, followed by
        // a key of keySizeInBytes bytes.  The key size must be between 4 and 8 bytes, inclusive.
        //
        // Because the size of the fields and the count of values is variable, we can't use a struct and are stuck
        // with ugly memory manipulation to access the hash table entries.
        //
        // The format is 1 or 2 (valueCount) values of size valueSize, followed by keySize bytes of key.
        //
#if 0
        struct Entry {
            unsigned        value1;
            unsigned        value2;
            unsigned char   key[1]; // Actual size of key determined by keySizeInBytes
        };
#endif // 0

        // Free Entries have the first valueSizeInByes bytes of value == invalidValueValue.  The following methods
        // understand the format and try to make it less opaque to use them.
        
        inline void *getEntry(_uint64 whichEntry) const {
            return ((char *)Table + elementSize * whichEntry);
        }

        inline bool doesEntryHaveInvalidValue(void *entry) const
        {
            return !memcmp(entry, &invalidValueValue, valueSizeInBytes);
        }

        inline ValueType getValueFromEntry(void *entry, unsigned whichValue) const
        {
            _ASSERT(whichValue < valueCount);
            ValueType value = 0;    // Need =0 because valueSizeInBytes might be < sizeeof(ValueType)
            memcpy(&value, (char *)entry + whichValue * valueSizeInBytes, valueSizeInBytes);    // Assumes little-endian
            return value;
        }

        inline bool isKeyEqual(const void *entry, KeyType key) const 
        {
            return !memcmp((const char *)entry + valueSizeInBytes * valueCount, &key, keySizeInBytes);
        }

        inline void clearKey(void *entry)
        {
            memset((char *)entry + valueSizeInBytes * valueCount, 0 , keySizeInBytes);
        }

        inline void setKey(void *entry, KeyType key)
        {
            memcpy((char *)entry + valueSizeInBytes * valueCount, &key, keySizeInBytes);
        }
 
        void *Table;
        size_t tableSize;
        unsigned keySizeInBytes;
        unsigned elementSize;
        size_t usedElementCount;
        bool ownsMemoryForTable;
        unsigned valueSizeInBytes;
        unsigned valueCount;
        ValueType invalidValueValue;
 
        //
        // Returns either the entry for this key, or else the entry where the key would be
        // inserted if it's not in the table.
        //
        void* getEntryForKey(__in KeyType key) const;

        friend class SeedCountIterator;

        static const unsigned magic;
};
