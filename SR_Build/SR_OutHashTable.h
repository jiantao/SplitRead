/*
 * =====================================================================================
 *
 *       Filename:  SR_OutHashTable.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/01/2011 08:26:38 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  REFHASHTABLE_H
#define  REFHASHTABLE_H

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>

#include "SR_Types.h"


#define DEFAULT_HASH_SIZE 7

#define DEFAULT_POS_ARR_CAPACITY 2000

#define SR_HashPosArrayGet(hashPosArray, index) (((hashPosArray)->hashPos)[(index)])

#define SR_HashPosArrayReset(hashPosArray) ((hashPosArray)->size = 0)

#define WriteHashSize(outputStream, hashSize)                                                \
    do                                                                                       \
    {                                                                                        \
        size_t writeSize = 0;                                                                \
        writeSize = fwrite(&(hashSize), sizeof(unsigned char), 1, (outputStream));           \
        if (writeSize != 1)                                                                  \
            SR_ErrQuit("ERROR: Cannot write the hash size into output file.\n");             \
    }while(0)

// generate a mask to clear the highest 2 bits in a hash key (the leftmost base pair)
#define GetHighEndMask(hashSize) ((uint32_t) 0xffffffff >> (34 - (2 * (hashSize))))

typedef struct SR_HashPosArray
{
    uint32_t* hashPos;      // an array to hold all positions found in the reference for a certain hash

    unsigned int size;      // number of positions

    unsigned int capacity;  // maximum number of positions that can be hold in the array
}SR_HashPosArray;

typedef struct SR_OutHashTable
{
    unsigned char chr;           // chromsome
    
    unsigned char hashSize;      // size of hash

    char md5[MD5_STR_LEN + 1];   // md5 checksum string

    SR_HashPosArray* hashPosTable;  // table holds all the position of hashes found in reference

    uint32_t  numPos;            // total number of hash positions found in reference

    uint32_t  numHashes;         // total number of different hashes
}SR_OutHashTable;


SR_HashPosArray* SR_HashPosArrayAlloc(unsigned int capacity);

void SR_HashPosArrayInit(SR_HashPosArray* hashPosArray, unsigned int capacity);

void SR_HashPosArrayFree(SR_HashPosArray* hashPosArray);

void SR_HashPosArrayPushBack(SR_HashPosArray* hashPosArray, uint32_t pos);


SR_OutHashTable* SR_OutHashTableAlloc(unsigned char hashSize);

void SR_OutHashTableFree(SR_OutHashTable* pHashTable);

void SR_OutHashTableLoad(SR_OutHashTable* pHashTable, const char* refSeq, const char* md5String, uint32_t refLen, unsigned char chr);

off_t SR_OutHashTableWrite(FILE* hashTableOutput, const SR_OutHashTable* pHashTable);

void SR_OutHashTableReset(SR_OutHashTable* pHashTable);

#endif  /*REFHASHTABLE_H*/
