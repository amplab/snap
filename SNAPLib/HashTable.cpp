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
#include "Error.h"
#include "GenericFile_Blob.h"

SNAPHashTable::SNAPHashTable(
    _int64      i_tableSize,
    unsigned    i_keySizeInBytes,
    unsigned    i_valueSizeInBytes,
    unsigned    i_valueCount,
    _uint64     i_invalidValueValue)
/*++

Routine Description:

    Constructor for a new, empty closed hash table.

Arguments:
    tableSize           - How many slots should the table have.
--*/
{
    keySizeInBytes = i_keySizeInBytes;
    valueSizeInBytes = i_valueSizeInBytes;
    valueCount = i_valueCount;
    invalidValueValue = i_invalidValueValue;
    elementSize = keySizeInBytes + valueSizeInBytes * valueCount;
    tableSize = i_tableSize;
    usedElementCount = 0;
    Table = NULL;

    if (tableSize <= 0) {
        tableSize = 0;
        return;
    }

	Table = BigAlloc(tableSize * elementSize);
    ownsMemoryForTable = true;

    //
    // Run through the table and set all of the first values to invalidValueValue, which means
    // unused.
    //

    for (size_t i = 0; i < tableSize; i++) {
        void *entry = getEntry(i);
		_ASSERT(entry >= Table && entry <= (char *)Table + tableSize * elementSize);
        clearKey(entry);
        memcpy(getEntry(i), &invalidValueValue, valueSizeInBytes);
    }
}

SNAPHashTable *SNAPHashTable::loadFromBlob(GenericFile_Blob *loadFile)
{
	SNAPHashTable *table = loadCommon(loadFile);

	size_t bytesMapped;
	table->Table = loadFile->mapAndAdvance(table->tableSize * table->elementSize, &bytesMapped);
	if (bytesMapped != table->tableSize * table->elementSize) {
		WriteErrorMessage("SNAPHashTable: unable to map table\n");
		soft_exit(1);
	}
	table->ownsMemoryForTable = false;

	return table;
}

SNAPHashTable *SNAPHashTable::loadFromGenericFile(GenericFile *loadFile)
{
	SNAPHashTable *table = loadCommon(loadFile);
	table->Table = BigAlloc(table->tableSize * table->elementSize);
	loadFile->read(table->Table, table->tableSize * table->elementSize);
	table->ownsMemoryForTable = true;

	return table;
}

SNAPHashTable *SNAPHashTable::loadCommon(GenericFile *loadFile)
{
    SNAPHashTable *table = new SNAPHashTable();

    unsigned fileMagic;
    if (sizeof(magic) != loadFile->read(&fileMagic, sizeof(magic))) {
        WriteErrorMessage("Magic number mismatch on hash table load.  %d != %d\n", fileMagic, magic);
        soft_exit(1);
    }

    if (fileMagic != magic) {
        WriteErrorMessage("SNAPHashTable: magic number mismatch.  Perhaps you have a corruped index.  %d != %d\n", fileMagic, magic);
        soft_exit(1);
    }
 
    if (sizeof(table->tableSize) != loadFile->read(&table->tableSize, sizeof(table->tableSize))) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable fread table size failed\n");
        soft_exit(1);
    }

    if (sizeof(table->usedElementCount) != loadFile->read(&table->usedElementCount, sizeof(table->usedElementCount))) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable fread used element count failed\n");
        soft_exit(1);
    }

    if (sizeof(table->keySizeInBytes) != loadFile->read(&table->keySizeInBytes, sizeof(table->keySizeInBytes))) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable fread keySizeInBytes size failed.  Perhaps this is an old format hash table and needs to be rebuilt.\n");
        soft_exit(1);
    }

    if (table->keySizeInBytes < 4 || table->keySizeInBytes > 8) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable Key size must be between 4 and 8 inclusive.  Perhaps this is an old format hash table and needs to be rebuilt.\n");
        soft_exit(1);
    }

    if (sizeof(table->valueSizeInBytes) != loadFile->read(&table->valueSizeInBytes, sizeof(table->valueSizeInBytes))) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable fread dataSizeInBytes size failed.  Perhaps this is an old format hash table and needs to be rebuilt.\n");
        soft_exit(1);
    }

    if (table->valueSizeInBytes == 0 || table->valueSizeInBytes > sizeof(_uint64)) {
        //
        // It must be at least one byte, because we need that much for the unused value value. The code stuffs
        // values into _uint64, so it can't be bigger than that.
        //
        WriteErrorMessage(
            "SNAPHashTable::SNAPHashTable value size in bytes (%d) must be between 1 and 8.  Perhaps you have a hash table from a future version of SNAP?  Or else it's corrupt.\n", table->valueSizeInBytes);
        soft_exit(1);
    }

    if (sizeof(table->valueCount) != loadFile->read(&table->valueCount, sizeof(table->valueCount))) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable: value count failed to read.\n");
        soft_exit(1);
    }

    if (table->valueCount == 0 || table->valueCount > 2) {
        // Technically, > 2 would work fine with the code, but SNAP doesn't use it, so the check is here to detect corruption.
        WriteErrorMessage("SNAPHashTable::SNAPHashTable: invalid value count (%d), possible corruption or bad file format.\n", table->valueCount);
        soft_exit(1);
    }

    table->invalidValueValue = 0;   // Need this in case valueSizeInBytes < sizeof(ValueType)
    if (table->valueSizeInBytes != loadFile->read(&table->invalidValueValue, table->valueSizeInBytes)) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable: unable to read invalid value value\n");
        soft_exit(1);
    }

    if (table->tableSize <= 0) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable Zero or negative hash table size\n");
        soft_exit(1);
    }

    table->elementSize = table->keySizeInBytes + table->valueSizeInBytes * table->valueCount;



    return table;
}

SNAPHashTable::~SNAPHashTable()
{
    if (ownsMemoryForTable) {
        BigDealloc(Table);
    }
}

    bool
SNAPHashTable::saveToFile(const char *saveFileName, size_t *bytesWritten)
{
    FILE *saveFile = fopen(saveFileName,"wb");
    if (saveFile == NULL) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable(%s) fopen failed\n",saveFileName);
        return false;
    }

    bool worked = saveToFile(saveFile, bytesWritten);
    fclose(saveFile);

    return worked;
}

bool
SNAPHashTable::saveToFile(FILE *saveFile, size_t *bytesWritten) 
{
    *bytesWritten = 0;
    if (1 != fwrite(&magic,sizeof(magic), 1, saveFile)) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable fwrite magic number failed\n");
        return false;
    }    
    (*bytesWritten) += sizeof(magic);
    
    if (1 != fwrite(&tableSize,sizeof(tableSize), 1, saveFile)) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable fwrite table size failed\n");
        return false;
    }
    (*bytesWritten) += sizeof(tableSize);

    if (1 != fwrite(&usedElementCount,sizeof(usedElementCount), 1, saveFile)) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable fwrite used element count size failed\n");
        return false;
    }
    (*bytesWritten) += sizeof(usedElementCount);

    if (1 != fwrite(&keySizeInBytes, sizeof(keySizeInBytes), 1, saveFile)) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable fwrite key size failed\n");
        return false;
    }
    (*bytesWritten) += sizeof(keySizeInBytes);

    if (1 != fwrite(&valueSizeInBytes, sizeof(valueSizeInBytes), 1, saveFile)) {
        WriteErrorMessage("SNAPHashTable::SNAPHashTable fwrite data size failed\n");
        return false;
    }
    (*bytesWritten) += sizeof(valueSizeInBytes);

    if (1 != fwrite(&valueCount, sizeof(valueCount), 1, saveFile)) {
        WriteErrorMessage("SNAPHashTable: fwrite value count failed\n");
        return false;
    }
    (*bytesWritten) += sizeof(valueCount);

    if (1 != fwrite(&invalidValueValue, valueSizeInBytes, 1, saveFile)) {
        WriteErrorMessage("SNAPHashTable: fwrite invalid value value failed\n");
        return false;
    }
    (*bytesWritten) += valueSizeInBytes;

    size_t maxWriteSize = 100 * 1024 * 1024;
    size_t writeOffset = 0;
    while (writeOffset < tableSize * elementSize) {
        size_t amountToWrite = __min(maxWriteSize,tableSize * elementSize - writeOffset);
        size_t thisWrite = fwrite((char*)Table + writeOffset, 1, amountToWrite, saveFile);
        if (thisWrite < amountToWrite) {
            WriteErrorMessage("SNAPHashTable::saveToFile: fwrite failed, %d\n"
                              "handle %p, addr %p, atr: %lu, &bw %p\n",errno, saveFile,(char*)Table + writeOffset, amountToWrite, &bytesWritten);
            return false;
        }
        writeOffset += thisWrite;
        (*bytesWritten) += thisWrite;
    }

    return true;
}
    
_int64 nCallsToGetEntryForKey = 0;
_int64 nProbesInGetEntryForKey = 0;

void *
SNAPHashTable::getEntryForKey(KeyType key) const
{
    nCallsToGetEntryForKey++;

    _uint64 tableIndex = hash(key) % tableSize;

    bool wrapped = false;
    _uint64 nProbes = 1;

    //
    // Chain through the table until we hit either a match on the key or an unused element
    //
    void *entry = getEntry(tableIndex);
    while (!isKeyEqual(entry, key) && !doesEntryHaveInvalidValue(entry)) {
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

        entry = getEntry(tableIndex);
    }

    nProbesInGetEntryForKey++;

    return entry;
}



    SNAPHashTable::ValueType * 
SNAPHashTable::SlowLookup(KeyType key)
{
    void *entry = getEntryForKey(key);

    if (NULL == entry || doesEntryHaveInvalidValue(entry)) {
        return NULL;
    }

    return (ValueType *)entry;
}

    bool 
SNAPHashTable::Insert(KeyType key, ValueType *data)
{
    void *entry = getEntryForKey(key);

    if (NULL == entry) {
        return false;
    }

	if (!isKeyEqual(entry, key)) {
		setKey(entry, key);
		usedElementCount++;
	}

    for (unsigned i = 0; i < valueCount; i++) {
        memcpy((char *)entry + i * valueSizeInBytes, &data[i], valueSizeInBytes);   // Assumes little endian
    }

    return true;
}



const unsigned SNAPHashTable::magic = 0xb111b010;
