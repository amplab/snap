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
        exit(1);
    }


    //
    // Run through the table and set all of the value1s to 0xffffffff, which means
    // unused.
    //

    for (unsigned i = 0; i < tableSize; i++) {
        table[i].key = 0;
        table[i].value1 = 0xffffffff;
    }
}

SNAPHashTable::SNAPHashTable(char *loadFileName)
/*++

Routine Description:

    Constructor to load a previously saved hash table.

Arguments:

    loadFileName    - the base file name you want to load (the saver actually makes two files and both are needed to load)

--*/
{
    useBigAlloc = true;

    FILE *loadFile = fopen(loadFileName,"rb");
    if (loadFile == NULL) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable(%s) fopen failed\n",loadFileName);
        exit(1);
    }

    if (1 != fread(&tableSize,sizeof(tableSize),1,loadFile)) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable(%s) fread table size failed\n",loadFileName);
        exit(1);
    }

    if (1 != fread(&bytesToCheckForUnusedEntry, sizeof(bytesToCheckForUnusedEntry),1,loadFile)) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable(%s) fread unused entry size failed\n",loadFileName);
        exit(1);
    }

    if (1 != fread(&usedElementCount,sizeof(usedElementCount),1,loadFile)) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable(%s) fread data size failed\n",loadFileName);
        exit(1);

    }

    if (tableSize <= 0) {
        table = NULL;
        tableSize = 0;
        fclose(loadFile);
        return;
    }

    table = (Entry *)BigAlloc(tableSize * elementSize, &virtualAllocSize);

    if (NULL == table) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable(%s) unable to allocate table memory\n",loadFileName);
        exit(1);
    }

    const char *tableFileExtension = ".table";
    size_t tableFileNameLength = strlen(loadFileName) + strlen(tableFileExtension) + 1;    
    char *tableFileName  = new char[tableFileNameLength];
    snprintf(tableFileName,tableFileNameLength,"%s%s",loadFileName,tableFileExtension);
    FILE* tableFile = fopen(tableFileName,"rb");
    if (tableFile == NULL) {
        fprintf(stderr,"CloseHashTable::SNAPHashTable: Unable to open table file %s, %d\n",tableFileName,errno);
        exit(1);
    }

    size_t maxReadSize = 100 * 1024 * 1024;
    size_t readOffset = 0;
    while (readOffset < tableSize * elementSize) {
        size_t amountToRead = __min(tableSize * elementSize - readOffset,
			__min(maxReadSize,virtualAllocSize - readOffset));
        size_t bytesRead = fread((char*)table + readOffset, 1, amountToRead, tableFile);
        if (bytesRead < amountToRead) {
            fprintf(stderr,"SNAPHashTable::SNAPHashTable: fread failed, %d, %lu, %lu\n",errno,bytesRead,amountToRead);
            exit(1);
        }
        // MATEI: Not sure this is needed
        /*
        if (0 == bytesRead) {
            fprintf(stderr,"SNAPHashTable::SNAPHashTable: fread read no data\n");
            exit(1);
        }
        */
        readOffset += bytesRead;
    }

    fclose(loadFile);
    fclose(tableFile);
    delete [] tableFileName;

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

    if (1 != fwrite(&tableSize,sizeof(tableSize),1,saveFile)) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable(%s) fwrite table size failed\n",saveFileName);
        return false;

    }

    if (1 != fwrite(&bytesToCheckForUnusedEntry,sizeof(bytesToCheckForUnusedEntry),1,saveFile)) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable(%s) fwrite unused entry size failed\n",saveFileName);
        return false;
    }

    if (1 != fwrite(&usedElementCount,sizeof(usedElementCount),1,saveFile)) {
        fprintf(stderr,"SNAPHashTable::SNAPHashTable(%s) fwrite data size failed\n",saveFileName);
        return false;

    }

    const char *tableFileExtension = ".table";
    size_t tableFileNameLength = strlen(saveFileName) + strlen(tableFileExtension) + 1;    
    char *tableFileName  = new char[tableFileNameLength];
    snprintf(tableFileName,tableFileNameLength,"%s%s",saveFileName,tableFileExtension);
    FILE* tableFile = fopen(tableFileName,"wb");
    if (tableFile == NULL) {
        fprintf(stderr,"CloseHashTable::saveToFile: Unable to open table file %s, %d\n",tableFileName,errno);
        exit(1);
    }

    size_t maxWriteSize = 100 * 1024 * 1024;
    size_t writeOffset = 0;
    while (writeOffset < tableSize * elementSize) {
        size_t amountToWrite = __min(maxWriteSize,virtualAllocSize - writeOffset);
        size_t bytesWritten = fwrite((char*)table + writeOffset, 1, amountToWrite, tableFile);
        if (bytesWritten < amountToWrite) {
            fprintf(stderr,"SNAPHashTable::saveToFile: fwrite failed, %d\n",errno);
            fprintf(stderr,"handle %p, addr %p, atr: %lu, &bw %p\n",tableFile,(char*)table + writeOffset,amountToWrite,&bytesWritten);
            return false;
        }
        // MATEI: Not sure this is needed
        /*
        if (0 == bytesWritten) {
            fprintf(stderr,"SNAPHashTable::saveToFile: WriteFile wrote no data\n");
            return false;
        }
        */
        writeOffset += bytesWritten;
    }

    fclose(tableFile);
    delete [] tableFileName;
    

    fclose(saveFile);
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
    while (table[tableIndex].key != key && table[tableIndex].value1 != 0xffffffff) {
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
    _ASSERT(data[0] != 0xffffffff); // This is the unused value that represents an empty hash table.  You can't use it.

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

    if (NULL == entry || entry->value1 == 0xffffffff) {
        return NULL;
    }

    return &entry->value1;
}
