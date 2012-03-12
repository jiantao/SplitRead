/*
 * =====================================================================================
 *
 *       Filename:  SR_ReadPairBuild.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/08/2011 11:44:06 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_READPAIRBUILD_H
#define  SR_READPAIRBUILD_H

#include "SR_LibInfo.h"

typedef struct
{
    double cutoff;

    double trimRate;

    unsigned int binLen;

    uint32_t detectSet;

    unsigned char minMQ;

    FILE* fileListInput;

    char* workingDir;

}SR_ReadPairBuildPars;

typedef struct SR_LocalPairInfo
{
    int32_t readGrpID;

    int32_t refID:16, upMapQ:8, downMapQ:8;

    int32_t upPos;

    int32_t upEnd;

    int32_t downPos;

    int32_t fragLen;

    int32_t upNumMM:16, downNumMM:16;

    int32_t fragLenQual:16, readPairType:8, pairMode:8;

}SR_LocalPair;

typedef struct SR_CrossPair
{
    int32_t readGrpID;

    int32_t upRefID:16, downRefID:16;
            
    int32_t upMapQ:16, downMapQ:16;

    int32_t upPos;

    int32_t upEnd;

    int32_t downPos;

    int32_t downEnd;

    int32_t upNumMM:16, downNumMM:16;

    int32_t fragLenQual:16, readPairType:8, pairMode:8;

}SR_CrossPair;

typedef struct SR_SpecialPair
{
    int32_t readGrpID;

    int16_t refID[2];

    int32_t pos[2];

    int32_t end[2];

    int32_t numSpeicalHits;

    int16_t numMM[2];
    
    uint8_t bestMQ[2];

    uint8_t secMQ[2];

    int32_t fragLenQual:16, readPairType:8, pairMode:8;

    uint32_t specialID;

}SR_SpecialPair;

typedef struct SR_LocalPairArray
{
    SR_LocalPair* data;

    uint64_t* chrCount;

    uint64_t size;

    uint64_t capacity;

}SR_LocalPairArray;

typedef struct SR_CrossPairArray
{
    SR_CrossPair* data;

    uint64_t* chrCount;

    uint64_t size;

    uint64_t capacity;

}SR_CrossPairArray;

typedef struct SR_SpecialPairArray
{
    SR_SpecialPair* data;

    uint64_t* chrCount;

    uint64_t size;

    uint64_t capacity;

}SR_SpecialPairArray;

typedef struct SR_SpecialPairTable
{
    char (*names)[3];

    void* nameHash;

    SR_SpecialPairArray array;

    SR_SpecialPairArray crossArray;

    uint32_t size;

    uint32_t capacity;

}SR_SpecialPairTable;

typedef struct SR_ReadPairTable
{
    SR_LocalPairArray* pLongPairArray;

    SR_LocalPairArray* pShortPairArray;

    SR_LocalPairArray* pReversedPairArray;

    SR_LocalPairArray* pInvertedPairArray;

    SR_CrossPairArray* pCrossPairArray;

    SR_SpecialPairTable* pSpecialPairTable;

    uint32_t detectSet;

    uint32_t numChr;

    uint64_t numPairs;

}SR_ReadPairTable;

SR_SpecialPairTable* SR_SpecialPairTableAlloc(unsigned int numChr);

void SR_SpecialPairTableFree(SR_SpecialPairTable* pSpecialPairTable);

SR_ReadPairTable* SR_ReadPairTableAlloc(uint32_t numChr, uint32_t detectSet);

void SR_ReadPairTableFree(SR_ReadPairTable* pInfoTable);


void SR_ReadPairBuild(const SR_ReadPairBuildPars* pBuildPars);

void SR_ReadPairTableClear(SR_ReadPairTable* pReadPairTable);

void SR_ReadPairTableUpdate(SR_ReadPairTable* pReadPairTable, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const SR_ZAtag* pZAtag,
                            const SR_PairStats* pPairStats, const SR_LibInfoTable* pLibTable, const SR_FragLenHistArray* pHistArray, uint8_t minMQ);


void* SR_ReadPairFilesOpen(const SR_LibInfoTable* pLibTable, uint32_t detectSet, const char* workingDir);

void SR_ReadPairFilesClose(void* pFileHash);

void SR_ReadPairTableWrite(const SR_ReadPairTable* pReadPairTable, void* pHash);

void SR_SpecialPairTableWirteID(const SR_SpecialPairTable* pSpecialTable, FILE* libOutput);

SR_Status SR_CheckWorkingDir(const char* workingDir);

#endif  /*SR_READPAIRBUILD_H*/
