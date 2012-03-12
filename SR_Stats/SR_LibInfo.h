/*
 * =====================================================================================
 *
 *       Filename:  SR_LibInfo.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/03/2012 02:17:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_LIBINFO_H
#define  SR_LIBINFO_H

#include <stdint.h>
#include <stdio.h>

#include "SR_Types.h"
#include "SR_BamHeader.h"
#include "SR_FragLenHist.h"

#define SR_NUM_READ_PAIR_TYPES 8

#define SR_NUM_SV_TYPES 5


// a map used to map the pair mode into its corresponding number
// negative value means invalid mode
static const int8_t SR_PairModeMap[16] = 
{ 
     -1, -1, 0, 1,

     -1, -1, 2, 3,

     4, 5, -1, -1,

     6, 7, -1, -1
};  

static const int8_t SR_SeqTechMap[4][8] = 
{
    {0, 1, 2, 3, 4, 5, 6, 7},

    {2, 3, 0, 1, 5, 4, 7, 6},

    {1, 0, 3, 2, 6, 7, 4, 5},

    {3, 2, 1, 0, 7, 6, 5, 4}
};

typedef enum
{
    PT_UNKNOWN    = -1,

    PT_NORMAL     = 0,

    PT_LONG       = 1,

    PT_SHORT      = 2,

    PT_REVERSED   = 3,

    PT_INVERTED3  = 4,

    PT_INVERTED5  = 5,

    PT_CROSS      = 6,

    PT_SPECIAL3   = 7,

    PT_SPECIAL5   = 8

}SV_ReadPairType;

static const int8_t SV_ReadPairTypeMap[2][8] = 
{
    {PT_INVERTED3, PT_NORMAL, PT_UNKNOWN, PT_INVERTED5, PT_UNKNOWN, PT_UNKNOWN, PT_REVERSED, PT_UNKNOWN},

    {PT_UNKNOWN, PT_UNKNOWN, PT_REVERSED, PT_UNKNOWN, PT_INVERTED3, PT_NORMAL, PT_UNKNOWN, PT_INVERTED5}
};

typedef enum
{
    SV_NORMAL             = 0,

    SV_DELETION           = 1,

    SV_TANDEM_DUP         = 2,

    SV_INVERSION          = 3,

    SV_INTER_CHR_TRNSLCTN = 4,

    SV_SPECIAL            = 5

}SV_EventType;

typedef enum
{
    ST_ILLUMINA = 0,

    ST_454 = 1,

    ST_SOLID = 2,

    ST_ILLUMINA_LONG = 3

}SR_SeqTech;

// the object used to hold the basic statistics of a pair of alignments
typedef struct SR_PairStats
{
    int32_t readGrpID;        // name of the read group

    int fragLen;              // fragment length of the pair

    SR_PairMode pairMode;     // orientation mode of the pair

}SR_PairStats;

typedef struct SR_ZAtag
{
    uint8_t bestMQ[2];

    uint8_t secMQ[2];

    int16_t numMM[2];

    int32_t numMappings[2];

    char spRef[2][3];

}SR_ZAtag;

typedef struct SR_AnchorInfo
{
    char** pAnchors;

    void* pAnchorHash;

    char* pMd5s;

    int32_t* pLength;

    uint32_t size;

    uint32_t capacity;

}SR_AnchorInfo;

typedef struct SR_SampleInfo
{
    char** pSamples;

    void* pSampleHash;

    double* pReadFraction;

    uint32_t size;

    uint32_t capacity;

}SR_SampleInfo;

typedef struct SR_LibInfo
{
    int32_t fragLenMedian;

    int32_t fragLenHigh;

    int32_t fragLenLow;

}SR_LibInfo;

typedef struct SR_LibInfoTable
{
    SR_AnchorInfo* pAnchorInfo;

    SR_SampleInfo* pSampleInfo;

    char** pReadGrps;

    void* pReadGrpHash;

    SR_LibInfo* pLibInfo;

    int32_t* pSampleMap;

    int8_t* pSeqTech;

    uint32_t fragLenMax;

    uint32_t size;

    uint32_t capacity;

    double cutoff;

    double trimRate;

}SR_LibInfoTable;


SR_AnchorInfo* SR_AnchorInfoAlloc(uint32_t capacity);

void SR_AnchorInfoFree(SR_AnchorInfo* pInfo);

SR_SampleInfo* SR_SampleInfoAlloc(uint32_t capacity);

void SR_SampleInfoFree(SR_SampleInfo* pInfo);

SR_LibInfoTable* SR_LibInfoTableAlloc(uint32_t capAnchor, uint32_t capSample, uint32_t capReadGrp);

void SR_LibInfoTableFree(SR_LibInfoTable* pTable);

int SR_GetNumMismatchFromBam(const bam1_t* pAlgn);

int SR_GetNumMismatchFromZA(const char* cigarStr, unsigned int cigarLen, const char* mdStr, unsigned int mdLen);

SR_Status SR_LoadZAtag(SR_ZAtag* pZAtag, const bam1_t* pAlignment);

SR_Status SR_LoadPairStats(SR_PairStats* pPairStats, const bam1_t* pAlignment, const SR_LibInfoTable* pTable);

//====================================================================
// function:
//      check if a pair of read is normal
//
// args:
//      1. ppUpAlgn: a pointer of the pointer to a bam node object
//                   for the alignment with smaller coordinate
//      1. ppDownAlgn: a pointer of the pointer to a bam node object
//                   for the alignment with greater coordinate
//
// return:
//      if the pair of read is normal, return TRUE; else, return
//      FALSE
//=====================================================================
SR_Bool SR_IsNormalPair(SR_PairStats* pPairStats, unsigned int* pBackHistIndex,
                        const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const SR_LibInfoTable* pTable, unsigned short minMQ);

static inline void SR_LibInfoTableSetCutoff(SR_LibInfoTable* pTable, double cutoff)
{
    pTable->cutoff = cutoff;
}

static inline void SR_LibInfoTableSetTrimRate(SR_LibInfoTable* pTable, double trimRate)
{
    pTable->trimRate = trimRate;
}

SR_Status SR_LibInfoTableSetRG(SR_LibInfoTable* pTable, unsigned int* oldSize, const SR_BamHeader* pBamHeader);

SR_Status SR_LibInfoTableGetRGIndex(int32_t* pReadGrpIndex, const SR_LibInfoTable* pTable, const char* pRreadGrpName);

SR_Status SR_LibInfoTableCheckPair(unsigned int* pHistIndex, const SR_LibInfoTable* pTable, SR_PairStats* pPairStats);

void SR_LibInfoTableUpdate(SR_LibInfoTable* pTable, const SR_FragLenHistArray* pHistArray, unsigned int oldSize);

void SR_LibInfoTableWrite(const SR_LibInfoTable* pTable, FILE* libFile);

SR_LibInfoTable* SR_LibInfoTableRead(FILE* libFile);

#endif  /*SR_LIBINFO_H*/
