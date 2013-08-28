/*++

Module Name:

	HashTable.cpp

Abstract:

    Large closed hash table, specialized for SNAP

Authors:

    Bill Bolosky, March, 2011

Environment:


--*/
#include "stdafx.h"
#include "HashTable.h"
#include "BigAlloc.h"
#include "exit.h"
#include "Genome.h"

SNAPHashTable::SNAPHashTable(
    unsigned i_tableSize,
    bool     i_useBigAlloc)
/*++

Routine Description:

    Constructor for a new, empty closed hash table.

Arguments:
    tableSize           - How many slots should the table have.
--*/
{
    tableSize = i_tableSize;
    useBigAlloc = i_useBigAlloc;
    usedElementCount = 0;
    table = NULL;

    if (tableSize <= 0) {
        tableSize = 0;
        return;
    }

#ifdef _MSC_VER
    if (useBigAlloc) {
        table = (Entry *)BigAlloc(tableSize * elementSize,&virtualAllocSize);
    } else {
        SYSTEM_INFO systemInfo[1];
        GetSystemInfo(systemInfo);

        virtualAllocSize = ((tableSize * elementSize + systemInfo->dwPageSize - 1) / systemInfo->dwPageSize) * systemInfo->dwPageSize;

        table = (Entry *)VirtualAlloc(0,virtualAllocSize,MEM_COMMIT|MEM_RESERVE,PAGE_READWRITE);
    }
#else
    table = (Entry *)BigAlloc(tableSize * elementSize,&virtualAllocSize);
#endif
    if (NULL == table) {
        fprintf(stderr,"SNAPHashTable: Unable to allocate table of size %lld",(_int64)tableSize * elementSize);
        soft_exit(1);
    }


    //
    // Run through the table and set all of the value1s to InvalidGenomeLocation, which means
    // unused.
    //

    for (unsigned i = 0; i < tableSize; i++) {
        table[i].key = 0;
        table[i].value1 = InvalidGenomeLocation;
    }
}

SNAPHashTable *SNAPHashTable::loadFromFile(char *loadFileName)
/*++

Routine Description:

    Constructor to load a previously saved hash table.

Arguments:

    loadFileName    - the base file name you want to load (the saver actually makes two files and both are needed to load)

--*/
{
    FILE *loadFile = fopen(loadFileName,"rb");
    if (loadFile == NULL) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable(%s) fopen failed\n",loadFileName);
        soft_exit(1);
    }

    SNAPHashTable *table = loadFromFile(loadFile);
    fclose (loadFile);

    return table;
}

SNAPHashTable *SNAPHashTable::loadFromFile(FILE *loadFile)
{
    SNAPHashTable *table = new SNAPHashTable();

    table->useBigAlloc = true;

    unsigned fileMagic;
    if (1 != fread(&fileMagic, sizeof(magic), 1, loadFile)) {
        fprintf(stderr,"Magic number mismatch on hash table load.  %d != %d\n", fileMagic, magic);
        soft_exit(1);
    }
 
    if (1 != fread(&table->tableSize, sizeof(table->tableSize), 1, loadFile)) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable fread table size failed\n");
        soft_exit(1);
    }


    if (1 != fread(&table->usedElementCount, sizeof(table->usedElementCount), 1, loadFile)) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable fread data size failed\n");
        soft_exit(1);

    }

    if (table->tableSize <= 0) {
        fprintf(stderr,"Zero or negative hash table size\n");
        soft_exit(1);
    }

    table->table = (Entry *)BigAlloc(table->tableSize * elementSize, &table->virtualAllocSize);

    size_t maxReadSize = 100 * 1024 * 1024;
    size_t readOffset = 0;
    while (readOffset < table->tableSize * elementSize) {
        size_t amountToRead = __min(table->tableSize * table->elementSize - readOffset,
			__min(maxReadSize,table->virtualAllocSize - readOffset));
        size_t bytesRead = fread((char*)table->table + readOffset, 1, amountToRead, loadFile);
        if (bytesRead < amountToRead) {
            fprintf(stderr,"SNAPHashTable::SNAPHashTable: fread failed, %d, %lu, %lu\n",errno,bytesRead,amountToRead);
            soft_exit(1);
        }
 
        readOffset += bytesRead;
    }

    return table;
}

SNAPHashTable::~SNAPHashTable()
{
#ifdef _MSC_VER
    if (useBigAlloc) {
        BigDealloc(table);
    } else {
        VirtualFree(table,0,MEM_RELEASE);
    }
#else
    BigDealloc(table);
#endif
}

    bool
SNAPHashTable::saveToFile(const char *saveFileName)
{
    FILE *saveFile = fopen(saveFileName,"wb");
    if (saveFile == NULL) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable(%s) fopen failed\n",saveFileName);
        return false;
    }

    bool worked = saveToFile(saveFile);
    fclose(saveFile);

    return worked;
}

bool
SNAPHashTable::saveToFile(FILE *saveFile) 
{
    if (1 != fwrite(&magic,sizeof(magic),1,saveFile)) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable fwrite magic number failed\n");
        return false;

    }    
    
    if (1 != fwrite(&tableSize,sizeof(tableSize),1,saveFile)) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable fwrite table size failed\n");
        return false;

    }

    if (1 != fwrite(&usedElementCount,sizeof(usedElementCount),1,saveFile)) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable fwrite data size failed\n");
        return false;

    }

    size_t maxWriteSize = 100 * 1024 * 1024;
    size_t writeOffset = 0;
    while (writeOffset < tableSize * elementSize) {
        size_t amountToWrite = __min(maxWriteSize,tableSize * elementSize - writeOffset);
        size_t bytesWritten = fwrite((char*)table + writeOffset, 1, amountToWrite, saveFile);
        if (bytesWritten < amountToWrite) {
            fprintf(stderr,"SNAPHashTable::saveToFile: fwrite failed, %d\n",errno);
            fprintf(stderr,"handle %p, addr %p, atr: %lu, &bw %p\n",saveFile,(char*)table + writeOffset,amountToWrite,&bytesWritten);
            return false;
        }
        writeOffset += bytesWritten;
    }

    return true;
}
    
_int64 nCallsToGetEntryForKey = 0;
_int64 nProbesInGetEntryForKey = 0;

SNAPHashTable::Entry *
SNAPHashTable::getEntryForKey(__in unsigned key) const
{
    nCallsToGetEntryForKey++;

    unsigned tableIndex = hash(key) % tableSize;

    bool wrapped = false;
    unsigned nProbes = 1;

    //
    // Chain through the table until we hit either a match on the key or an unused element
    //
    while (table[tableIndex].key != key && table[tableIndex].value1 != InvalidGenomeLocation) {
        nProbesInGetEntryForKey++;

        if (nProbes < QUADRATIC_CHAINING_DEPTH) {
            tableIndex += nProbes * nProbes;
        } else {
            tableIndex++;
        }

        nProbes++;

        if (tableIndex >= tableSize) {
            if (wrapped) {
                return NULL;
            }
            wrapped = true;
            tableIndex = tableIndex % tableSize;
        }
    }

    nProbesInGetEntryForKey++;

    return &table[tableIndex];
}

bool
SNAPHashTable::Insert(unsigned key, const unsigned *data)
{
    _ASSERT(data[0] != InvalidGenomeLocation); // This is the unused value that represents an empty hash table.  You can't use it.

    Entry *entry = getEntryForKey(key);
    if (NULL == entry) {
        return false;
    }

    if (entry->key != key) {
        _ASSERT(entry->key == 0);
        entry->key = key;
        usedElementCount++;
    }

    entry->value1 = data[0];
    entry->value2 = data[1];

    return true;
}

_uint64
SNAPHashTable::GetHashTableMemorySize()
{
    return virtualAllocSize;
}

unsigned *
SNAPHashTable::SlowLookup(unsigned key)
{
    Entry *entry = getEntryForKey(key);

    if (NULL == entry || entry->value1 == InvalidGenomeLocation) {
        return NULL;
    }

    return &entry->value1;
}

const unsigned SNAPHashTable::magic = 0xb111b010;
