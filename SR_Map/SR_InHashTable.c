/*
 * =====================================================================================
 *
 *       Filename:  SR_InHashTable.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/13/2011 00:49:19
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>

#include "SR_Error.h"
#include "SR_InHashTable.h"

SR_InHashTable* SR_InHashTableAlloc(unsigned char hashSize)
{
    SR_InHashTable* pNewTable = (SR_InHashTable*) malloc(sizeof(SR_InHashTable));
    if (pNewTable == NULL)
        SR_ErrSys("ERROR: Not enough memory for a reference hash table object.\n");
    
    pNewTable->id = 0;
    pNewTable->hashSize = hashSize;

    pNewTable->highEndMask = GET_HIGH_END_MASK(hashSize);
    pNewTable->numHashes = (uint32_t) 1 << (2 * hashSize);
    pNewTable->indices = (uint32_t*) malloc(sizeof(uint32_t) * pNewTable->numHashes);
    if (pNewTable->indices == NULL)
        SR_ErrSys("ERROR: Not enough memory for the hash index array in a hash table object.\n");

    pNewTable->numPos = 0;
    pNewTable->hashPos = NULL;

    return pNewTable;
}

void SR_InHashTableFree(SR_InHashTable* pHashTable)
{
    if (pHashTable != NULL)
    {
        free(pHashTable->hashPos);
        free(pHashTable->indices);

        free(pHashTable);
    }
}

int64_t SR_InHashTableReadStart(unsigned char* pHashSize, FILE* htInput)
{
    size_t readSize = 0;
    int64_t refHeaderPos = 0;

    readSize = fread(&refHeaderPos, sizeof(refHeaderPos), 1, htInput);
    if (readSize != 1)
        SR_ErrSys("ERROR: Cannot read the reference header position from the hash table file.\n");

    readSize = fread(pHashSize, sizeof(*pHashSize), 1, htInput);
    if (readSize != 1)
        SR_ErrSys("ERROR: Cannot read the hash size from the hash table file.\n");

    return refHeaderPos;
}

SR_Bool SR_InHashTableJump(FILE* htInput, const SR_RefHeader* pRefHeader, int32_t refID)
{
    int64_t jumpPos = pRefHeader->htFilePos[refID];

    if (fseeko(htInput, jumpPos, SEEK_SET) != 0)
        return FALSE;

    return TRUE;
}

SR_Bool SR_InHashTableRead(SR_InHashTable* pHashTable, FILE* htInput)
{
    size_t readSize = 0;

    readSize = fread(&(pHashTable->id), sizeof(pHashTable->id), 1, htInput);
    if (readSize != 1)
    {
        if (feof(htInput))
            return FALSE;
        else
            SR_ErrSys("ERROR: Cannot read the chromsome number from the hash table file.\n");
    }

    readSize = fread(pHashTable->indices, sizeof(uint32_t), pHashTable->numHashes, htInput);
    if (readSize != pHashTable->numHashes)
        SR_ErrSys("ERROR: Cannot read the indices from the hash table file.\n");

    readSize = fread(&(pHashTable->numPos), sizeof(uint32_t), 1, htInput);
    if (readSize != 1)
        SR_ErrSys("ERROR: Cannot read the total number of hash positions from the hash table file.\n");

    free(pHashTable->hashPos);
    pHashTable->hashPos = (uint32_t*) malloc(sizeof(uint32_t) * pHashTable->numPos);
    if (pHashTable->hashPos == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of hash positions in the hash table object.\n");

    readSize = fread(pHashTable->hashPos, sizeof(uint32_t), pHashTable->numPos, htInput);
    if (readSize != pHashTable->numPos)
        SR_ErrSys("ERROR: Cannot read the hash positions from the hash table.\n");

    return TRUE;
}


SR_Bool SR_InHashTableSearch(HashPosView* hashPosView, const SR_InHashTable* pHashTable, uint32_t hashKey)
{
    if(hashKey >= pHashTable->numHashes)
        SR_ErrSys("ERROR: Invalid hash key.\n");

    uint32_t index = pHashTable->indices[hashKey];
    uint32_t nextIndex = hashKey == (pHashTable->numHashes - 1) ? pHashTable->numPos : pHashTable->indices[hashKey + 1];

    if (index == nextIndex)
        return FALSE;
    
    hashPosView->size = nextIndex - index;
    hashPosView->data = pHashTable->hashPos + index;

    return TRUE;
}
