/*
 * =====================================================================================
 *
 *       Filename:  SR_ReadPairBuild.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/08/2011 11:45:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <math.h>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>

#include "khash.h"
#include "SR_Error.h"
#include "SR_Utilities.h"
#include "SR_BamPairAux.h"
#include "SR_BamInStream.h"
#include "SR_ReadPairBuild.h"

#define DEFAULT_RP_INFO_CAPACITY 50

#define SR_NUM_RP_FILE 6

#define DEFAULT_SP_TABLE_CAPACITY 10

const char* SR_LibTableFileName = "lib_table.dat";

const char* SR_HistFileName = "hist.dat";

const char* SR_READ_PAIR_FILE_NAME_TEMPLATE[] = 
{
    "refXXX_long_pairs.dat",

    "refXXX_short_pairs.dat",

    "refXXX_reversed_pairs.dat",

    "refXXX_inverted_pairs.dat",

    "refXXX_cross_pairs.dat",

    "refXXX_special_pairs.dat",
};

enum
{
    SR_LONG_PAIR_FILE = 0,

    SR_SHORT_PAIR_FILE = 1,

    SR_REVERSED_PAIR_FILE = 2,

    SR_INVERTED_PAIR_FILE = 3,

    SR_CROSS_PAIR_FILE = 4,

    SR_SPEICAL_PAIR_FILE = 5
};

typedef struct SR_ReadPairOutStream
{
    FILE* output;

    int64_t numPairs;

}SR_ReadPairOutStream;

KHASH_MAP_INIT_STR(name, uint32_t);

KHASH_MAP_INIT_INT(file, SR_ReadPairOutStream);

static SV_ReadPairType SR_CheckReadPairType(const SR_ZAtag* pZAtag, const SR_PairStats* pPairStats, const SR_LibInfoTable* pLibTable)
{
    if (pZAtag != NULL)
    {
        if (pZAtag->spRef[0][0] != ' ' && pZAtag->spRef[1][0] != ' ')
            return PT_UNKNOWN;
        else if (pZAtag->spRef[0][0] == ' ' && pZAtag->spRef[1][0] != ' ')
        {
            if (pZAtag->numMappings[0] == 1)
                return PT_SPECIAL3;
            else
                return PT_UNKNOWN;
        }
        else if (pZAtag->spRef[0][0] != ' ' && pZAtag->spRef[1][0] == ' ')
        {
            if (pZAtag->numMappings[1] == 1)
                return PT_SPECIAL5;
            else
                return PT_UNKNOWN;
        }

        if (pZAtag->numMappings[0] != 1 || pZAtag->numMappings[1] != 1)
            return PT_UNKNOWN;
    }

    // two mates of a read pair aligned to different references
    if (pPairStats->fragLen == -1)
        return PT_CROSS;

    SV_ReadPairType type1 = SV_ReadPairTypeMap[0][pPairStats->pairMode];
    SV_ReadPairType type2 = SV_ReadPairTypeMap[1][pPairStats->pairMode];

    if (type1 == PT_NORMAL || type2 == PT_NORMAL)
    {
        if (pPairStats->fragLen > pLibTable->pLibInfo[pPairStats->readGrpID].fragLenHigh)
            return PT_LONG;
        else if (pPairStats->fragLen < pLibTable->pLibInfo[pPairStats->readGrpID].fragLenLow)
            return PT_SHORT;
        else
            return PT_NORMAL;
    }
    else if ((type1 == PT_UNKNOWN && type2 == PT_UNKNOWN) || (type1 != PT_UNKNOWN && type2 != PT_UNKNOWN))
    {
        return PT_UNKNOWN;
    }
    else
    {
        return (type1 == PT_UNKNOWN ? type2 : type1);
    }
}

static void SR_LocalPairArrayUpdate(SR_LocalPairArray* pLocalPairArray, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const SR_LibInfoTable* pLibTable, 
        const SR_FragLenHistArray* pHistArray, const SR_PairStats* pPairStats, SV_ReadPairType readPairType)
{
    if (pLocalPairArray->size == pLocalPairArray->capacity)
        SR_ARRAY_RESIZE(pLocalPairArray, pLocalPairArray->capacity * 2, SR_LocalPair);

    SR_LocalPair* pLocalPair = pLocalPairArray->data + pLocalPairArray->size;

    pLocalPair->readGrpID = pPairStats->readGrpID;
    pLocalPair->refID = pUpAlgn->core.tid;

    pLocalPair->upPos = pUpAlgn->core.pos;
    pLocalPair->downPos = pDownAlgn->core.pos;

    pLocalPair->fragLen = pPairStats->fragLen;
    pLocalPair->upEnd = bam_calend(&(pUpAlgn->core), bam1_cigar(pUpAlgn));

    pLocalPair->upNumMM = SR_GetNumMismatchFromBam(pUpAlgn);
    pLocalPair->downNumMM = SR_GetNumMismatchFromBam(pDownAlgn);

    pLocalPair->upMapQ = pUpAlgn->core.qual;
    pLocalPair->downMapQ = pDownAlgn->core.qual;

    pLocalPair->pairMode = pPairStats->pairMode;
    pLocalPair->readPairType = readPairType;

    pLocalPair->fragLenQual = SR_FragLenHistArrayGetFragLenQual(pHistArray, pLibTable->size - pPairStats->readGrpID, pPairStats->fragLen);

    ++(pLocalPairArray->chrCount[pLocalPair->refID]);
    ++(pLocalPairArray->size);
}

static void SR_InvertedPairArrayUpdate(SR_LocalPairArray* pInvertedPairArray, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, 
        const SR_PairStats* pPairStats, SV_ReadPairType readPairType)
{
    if (pInvertedPairArray->size == pInvertedPairArray->capacity)
        SR_ARRAY_RESIZE(pInvertedPairArray, pInvertedPairArray->capacity * 2, SR_LocalPair);

    SR_LocalPair* pInvertedPair = pInvertedPairArray->data + pInvertedPairArray->size;

    pInvertedPair->readGrpID = pPairStats->readGrpID;
    pInvertedPair->refID = pUpAlgn->core.tid;

    pInvertedPair->upPos = pUpAlgn->core.pos;
    pInvertedPair->downPos = pDownAlgn->core.pos;

    pInvertedPair->fragLen = pPairStats->fragLen;
    pInvertedPair->upEnd = bam_calend(&(pUpAlgn->core), bam1_cigar(pUpAlgn));

    pInvertedPair->upNumMM = SR_GetNumMismatchFromBam(pUpAlgn);
    pInvertedPair->downNumMM = SR_GetNumMismatchFromBam(pDownAlgn);

    pInvertedPair->upMapQ = pUpAlgn->core.qual;
    pInvertedPair->downMapQ = pDownAlgn->core.qual;

    pInvertedPair->pairMode = pPairStats->pairMode;
    pInvertedPair->readPairType = readPairType;

    pInvertedPair->fragLenQual = INVALID_FRAG_LEN_QUAL;

    ++(pInvertedPairArray->chrCount[pInvertedPair->refID]);
    ++(pInvertedPairArray->size);
}

static void SR_CrossPairArrayUpdate(SR_CrossPairArray* pCrossPairArray, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const SR_PairStats*  pPairStats)
{
    if (pCrossPairArray->size == pCrossPairArray->capacity)
        SR_ARRAY_RESIZE(pCrossPairArray, pCrossPairArray->capacity * 2, SR_CrossPair);

    SR_CrossPair* pCrossPair = pCrossPairArray->data + pCrossPairArray->size;

    pCrossPair->readGrpID = pPairStats->readGrpID;

    pCrossPair->upRefID = pUpAlgn->core.tid;
    pCrossPair->downRefID = pDownAlgn->core.tid;

    pCrossPair->upPos = pUpAlgn->core.pos;
    pCrossPair->downPos = pDownAlgn->core.pos;

    pCrossPair->upEnd = bam_calend(&(pUpAlgn->core), bam1_cigar(pUpAlgn));
    pCrossPair->downEnd = bam_calend(&(pDownAlgn->core), bam1_cigar(pDownAlgn));

    pCrossPair->upNumMM = SR_GetNumMismatchFromBam(pUpAlgn);
    pCrossPair->downNumMM = SR_GetNumMismatchFromBam(pDownAlgn);

    pCrossPair->upMapQ = pUpAlgn->core.qual;
    pCrossPair->downMapQ = pDownAlgn->core.qual;

    pCrossPair->pairMode = pPairStats->pairMode;
    pCrossPair->readPairType = PT_CROSS;

    pCrossPair->fragLenQual = INVALID_FRAG_LEN_QUAL;

    ++(pCrossPairArray->chrCount[pCrossPair->upRefID]);
    ++(pCrossPairArray->size);
}

static void SR_SpecialPairTableUpdate(SR_SpecialPairTable* pSpecialPairTable, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const SR_LibInfoTable* pLibTable, 
        const SR_FragLenHistArray* pHistArray, const SR_PairStats* pPairStats, const SR_ZAtag* pZAtag, SV_ReadPairType readPairType)
{
    if (pSpecialPairTable->size == pSpecialPairTable->capacity)
    {
        pSpecialPairTable->capacity *= 2;
        pSpecialPairTable->names = (char (*)[3]) realloc(pSpecialPairTable->names, sizeof(char) * pSpecialPairTable->capacity);
        if (pSpecialPairTable->names == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the special reference names.\n");

        for (unsigned int i = pSpecialPairTable->size; i != pSpecialPairTable->capacity; ++i)
            pSpecialPairTable->names[i][2] = '\0';
    }

    int anchorIndex = readPairType - PT_SPECIAL3;
    int specialIndex = anchorIndex ^ 1;

    pSpecialPairTable->names[pSpecialPairTable->size][0] = pZAtag->spRef[specialIndex][0];
    pSpecialPairTable->names[pSpecialPairTable->size][1] = pZAtag->spRef[specialIndex][1];

    int ret = 0;
    int spRefID = 0;
    khiter_t khIter = kh_put(name, pSpecialPairTable->nameHash, pSpecialPairTable->names[pSpecialPairTable->size], &ret);

    if (ret != 0)
    {
        kh_value((khash_t(name)*) pSpecialPairTable->nameHash, khIter) = pSpecialPairTable->size;
        spRefID = pSpecialPairTable->size;
        ++(pSpecialPairTable->size);
    }
    else
        spRefID = kh_value((khash_t(name)*) pSpecialPairTable->nameHash, khIter);

    const bam1_t* pAlgns[2] = {pUpAlgn, pDownAlgn};
    SR_SpecialPairArray* pSpecialPairArray = NULL;

    if (pAlgns[anchorIndex]->core.tid == pUpAlgn->core.tid)
        pSpecialPairArray = &(pSpecialPairTable->array);
    else
        pSpecialPairArray = &(pSpecialPairTable->crossArray);

    if (pSpecialPairArray->capacity == pSpecialPairArray->size)
        SR_ARRAY_RESIZE(pSpecialPairArray, pSpecialPairArray->capacity * 2, SR_SpecialPair);

    SR_SpecialPair* pSpecialPair = pSpecialPairArray->data + pSpecialPairArray->size;

    pSpecialPair->readGrpID = pPairStats->readGrpID;

    pSpecialPair->refID[0] = pAlgns[anchorIndex]->core.tid;
    pSpecialPair->refID[1] = pAlgns[specialIndex]->core.tid;

    pSpecialPair->pos[0] = pAlgns[anchorIndex]->core.pos;
    pSpecialPair->pos[1] = pAlgns[specialIndex]->core.pos;

    pSpecialPair->end[0] = bam_calend(&(pAlgns[anchorIndex]->core), bam1_cigar(pAlgns[anchorIndex]));
    pSpecialPair->end[1] = bam_calend(&(pAlgns[specialIndex]->core), bam1_cigar(pAlgns[specialIndex]));

    pSpecialPair->numMM[0] = pZAtag->numMM[anchorIndex];
    pSpecialPair->numMM[1] = pZAtag->numMM[specialIndex];

    // FIXME: the mapping quality in the ZA tag is useless for current version.
    // use the mapping quality from the alignment instead
    pSpecialPair->bestMQ[0] = pAlgns[anchorIndex]->core.qual;
    pSpecialPair->bestMQ[1] = pAlgns[specialIndex]->core.qual;

    pSpecialPair->secMQ[0] = pZAtag->secMQ[anchorIndex];
    pSpecialPair->secMQ[1] = pZAtag->secMQ[specialIndex];

    pSpecialPair->numSpeicalHits = pZAtag->numMappings[specialIndex];
    pSpecialPair->pairMode = pPairStats->pairMode;
    pSpecialPair->readPairType = readPairType;
    pSpecialPair->specialID = spRefID;

    if (SV_ReadPairTypeMap[0][pPairStats->pairMode] == PT_NORMAL || SV_ReadPairTypeMap[1][pPairStats->pairMode] == PT_NORMAL)
        pSpecialPair->fragLenQual = SR_FragLenHistArrayGetFragLenQual(pHistArray, pLibTable->size - pPairStats->readGrpID, pPairStats->fragLen);
    else
        pSpecialPair->fragLenQual = INVALID_FRAG_LEN_QUAL;

    if (pAlgns[anchorIndex]->core.tid == pUpAlgn->core.tid)
        ++(pSpecialPairArray->chrCount[pSpecialPair->refID[0]]);

    ++(pSpecialPairArray->size);
}

void SR_SpecialPairTableClear(SR_SpecialPairTable* pSpecialPairTable, unsigned int numChr)
{
    pSpecialPairTable->array.size = 0;
    memset(pSpecialPairTable->array.chrCount, 0, sizeof(uint64_t) * numChr);

    pSpecialPairTable->crossArray.size = 0;
}

SR_SpecialPairTable* SR_SpecialPairTableAlloc(unsigned int numChr)
{
    SR_SpecialPairTable* pSpecialPairTable = (SR_SpecialPairTable*) malloc(sizeof(SR_SpecialPairTable));
    if (pSpecialPairTable == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a special pair table object.\n");

    pSpecialPairTable->size = 0;
    pSpecialPairTable->capacity = DEFAULT_SP_TABLE_CAPACITY;

    pSpecialPairTable->names = (char (*)[3]) calloc(sizeof(char), DEFAULT_SP_TABLE_CAPACITY * 3);
    if (pSpecialPairTable->names == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the special reference namese.\n");

    pSpecialPairTable->nameHash = kh_init(name);

    SR_SpecialPairArray* pArray = &(pSpecialPairTable->array);
    SR_ARRAY_INIT(pArray, DEFAULT_RP_INFO_CAPACITY, SR_SpecialPair);

    SR_SpecialPairArray* pCrossArray = &(pSpecialPairTable->crossArray);
    SR_ARRAY_INIT(pCrossArray, DEFAULT_RP_INFO_CAPACITY, SR_SpecialPair);

    pSpecialPairTable->array.chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
    if (pSpecialPairTable->array.chrCount == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");

    return pSpecialPairTable;
}

void SR_SpecialPairTableFree(SR_SpecialPairTable* pSpecialPairTable)
{
    if (pSpecialPairTable != NULL)
    {
        free(pSpecialPairTable->array.chrCount);
        SR_ARRAY_FREE(&(pSpecialPairTable->array), FALSE);

        SR_ARRAY_FREE(&(pSpecialPairTable->crossArray), FALSE);

        free(pSpecialPairTable->names);
        kh_destroy(name, pSpecialPairTable->nameHash);

        free(pSpecialPairTable);
    }
}

SR_ReadPairTable* SR_ReadPairTableAlloc(uint32_t numChr, uint32_t detectSet)
{
    SR_ReadPairTable* pNewInfoTable = (SR_ReadPairTable*) calloc(sizeof(SR_ReadPairTable), 1);
    if (pNewInfoTable == NULL)
        SR_ErrQuit("ERROR: Not enough memory for read pair information table.\n");

    pNewInfoTable->detectSet = detectSet;

    if ((detectSet & (1 << SV_DELETION)) != 0)
    {
        SR_ARRAY_ALLOC(pNewInfoTable->pLongPairArray, DEFAULT_RP_INFO_CAPACITY, SR_LocalPairArray, SR_LocalPair);
        pNewInfoTable->pLongPairArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
        if (pNewInfoTable->pLongPairArray->chrCount == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");
    }

    if ((detectSet & (1 << SV_TANDEM_DUP)) != 0)
    {
        SR_ARRAY_ALLOC(pNewInfoTable->pShortPairArray, DEFAULT_RP_INFO_CAPACITY, SR_LocalPairArray, SR_LocalPair);
        pNewInfoTable->pShortPairArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
        if (pNewInfoTable->pShortPairArray->chrCount == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");

        SR_ARRAY_ALLOC(pNewInfoTable->pReversedPairArray, DEFAULT_RP_INFO_CAPACITY, SR_LocalPairArray, SR_LocalPair);
        pNewInfoTable->pReversedPairArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
        if (pNewInfoTable->pReversedPairArray->chrCount == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");
    }

    if ((detectSet & (1 << SV_INVERSION)) != 0)
    {
        SR_ARRAY_ALLOC(pNewInfoTable->pInvertedPairArray, DEFAULT_RP_INFO_CAPACITY, SR_LocalPairArray, SR_LocalPair);
        pNewInfoTable->pInvertedPairArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
        if (pNewInfoTable->pInvertedPairArray->chrCount == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");
    }

    if ((detectSet & (1 << SV_INVERSION)) != 0)
    {
        SR_ARRAY_ALLOC(pNewInfoTable->pCrossPairArray, DEFAULT_RP_INFO_CAPACITY, SR_CrossPairArray, SR_CrossPair);
        pNewInfoTable->pCrossPairArray->chrCount = (uint64_t*) calloc(sizeof(uint64_t), numChr);
        if (pNewInfoTable->pCrossPairArray->chrCount == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the reference index.\n");
    }

    if ((detectSet & (1 << SV_SPECIAL)) != 0)
        pNewInfoTable->pSpecialPairTable = SR_SpecialPairTableAlloc(numChr);

    pNewInfoTable->numChr = numChr;
    pNewInfoTable->numPairs = 0;
    return pNewInfoTable;
}

void SR_ReadPairTableFree(SR_ReadPairTable* pReadPairTable)
{
    if (pReadPairTable != NULL)
    {
        if (pReadPairTable->pLongPairArray != NULL)
        {
            free(pReadPairTable->pLongPairArray->chrCount);
            SR_ARRAY_FREE(pReadPairTable->pLongPairArray, TRUE);
        }

        if (pReadPairTable->pShortPairArray != NULL)
        {
            free(pReadPairTable->pShortPairArray->chrCount);
            SR_ARRAY_FREE(pReadPairTable->pShortPairArray, TRUE);
        }

        if (pReadPairTable->pReversedPairArray != NULL)
        {
            free(pReadPairTable->pReversedPairArray->chrCount);
            SR_ARRAY_FREE(pReadPairTable->pReversedPairArray, TRUE);
        }

        if (pReadPairTable->pInvertedPairArray != NULL)
        {
            free(pReadPairTable->pInvertedPairArray->chrCount);
            SR_ARRAY_FREE(pReadPairTable->pInvertedPairArray, TRUE);
        }

        if (pReadPairTable->pCrossPairArray != NULL)
        {
            free(pReadPairTable->pCrossPairArray->chrCount);
            SR_ARRAY_FREE(pReadPairTable->pCrossPairArray, TRUE);
        }

        SR_SpecialPairTableFree(pReadPairTable->pSpecialPairTable);

        free(pReadPairTable);
    }
}

void SR_ReadPairBuild(const SR_ReadPairBuildPars* pBuildPars)
{
    // some default capacity of the containers
    unsigned int buffCapacity = 100;
    unsigned int capAnchor = 150;
    unsigned int capSample = 20;
    unsigned int capReadGrp = 20;
    unsigned int capHist = 10;

    // phony parameters for the bam input stream
    unsigned int reportSize = 0;
    unsigned int numThread = 0;

    // initialize the library information table
    SR_LibInfoTable* pLibTable = SR_LibInfoTableAlloc(capAnchor, capSample, capReadGrp);
    SR_LibInfoTableSetCutoff(pLibTable, pBuildPars->cutoff);
    SR_LibInfoTableSetTrimRate(pLibTable, pBuildPars->trimRate);

    // this is the data used to filter the read pairs in the bam
    SR_FilterDataRP* pFilterData = SR_FilterDataRPAlloc(pLibTable->pAnchorInfo, pBuildPars->binLen);

    // set the stream mode to read pair filter
    SR_StreamMode streamMode;
    SR_SetStreamMode(&streamMode, SR_ReadPairFilter, pFilterData, SR_READ_PAIR_MODE);

    // structure initialization
    SR_BamInStream* pBamInStream = SR_BamInStreamAlloc(pBuildPars->binLen, numThread, buffCapacity, reportSize, &streamMode);
    SR_FragLenHistArray* pHistArray = SR_FragLenHistArrayAlloc(capHist);
    SR_BamHeader* pBamHeader = NULL;

    // a flag to indicate if we have created the read pair table
    SR_Bool hasReadPairTable = FALSE;
    SR_ReadPairTable* pReadPairTable = NULL;

    // file hash
    khash_t(file)* pFileHash = NULL;

    // required arguments for bam in stream structure
    char* histOutputFile = SR_CreateFileName(pBuildPars->workingDir, SR_HistFileName);

    // open the fragment length histogram output file
    FILE* histOutput = fopen(histOutputFile, "w");
    if (histOutput == NULL)
        SR_ErrQuit("ERROR: Cannot open fragment length histogram file: %s\n", histOutputFile);

    char bamFileName[SR_MAX_LINE];
    while (SR_GetNextLine(bamFileName, SR_MAX_LINE, pBuildPars->fileListInput) == SR_OK)
    {
        // update the fitler data with new bam file
        // this will affect the bam in stream
        SR_FilterDataRPInit(pFilterData, bamFileName);
        SR_FilterDataRPTurnOffCross(pFilterData);

        // open the bam file
        SR_BamInStreamOpen(pBamInStream, bamFileName);

        // load the bam header before read any alignments
        pBamHeader = SR_BamInStreamLoadHeader(pBamInStream);

        // get the file position of the bam aligments
        int64_t bamPos = SR_BamInStreamTell(pBamInStream);

        // process the header information
        unsigned int oldSize = 0;
        if (SR_LibInfoTableSetRG(pLibTable, &oldSize, pBamHeader) != SR_OK)
            SR_ErrQuit("ERROR: Found an error when loading the bam file.\n");

        // initialize the fragment length histogram array with the number of newly added libraries in the bam file
        SR_FragLenHistArrayInit(pHistArray, pLibTable->size - oldSize);

        SR_BamNode* pUpNode = NULL;
        SR_BamNode* pDownNode = NULL;

        const bam1_t* pUpAlgn = NULL;
        const bam1_t* pDownAlgn = NULL;

        SR_Status bamStatus = SR_OK;
        while ((bamStatus = SR_BamInStreamLoadPair(&pUpNode, &pDownNode, pBamInStream)) != SR_EOF && bamStatus != SR_ERR)
        {
            // we hit another chromosome
            if (bamStatus == SR_OUT_OF_RANGE)
                continue;

            // load the read pairs
            if (!pFilterData->isFilled)
            {
                // usually we should catch the read pairs here
                // since their fragment length should be smaller than the bin length
                pUpAlgn = &(pUpNode->alignment);
                pDownAlgn = &(pDownNode->alignment);
            }
            else
            {
                // if the fragment length of the read pair is beyond our bin length
                // we should catch them here
                pUpAlgn = pFilterData->pUpAlgn;
                pDownAlgn = pFilterData->pDownAlgn;
            }

            // check if the incoming read pair is normal (unique-unique pair)
            // if yes, then update the corresponding fragment length histogram
            SR_PairStats pairStats;
            unsigned int backHistIndex = 0;
            if (SR_IsNormalPair(&pairStats, &backHistIndex, pUpAlgn, pDownAlgn, pLibTable, pBuildPars->minMQ))
                SR_FragLenHistArrayUpdate(pHistArray, backHistIndex, pairStats.fragLen);

            // recycle those bam nodes that are allocated from the memory pool
            if (!pFilterData->isFilled)
            {
                SR_BamInStreamRecycle(pBamInStream, pUpNode);
                SR_BamInStreamRecycle(pBamInStream, pDownNode);
            }
        }

        // finish the process of the histogram and update the library information table
        SR_FragLenHistArrayFinalize(pHistArray);
        SR_LibInfoTableUpdate(pLibTable, pHistArray, oldSize);

        // write the fragment length histogram into the file
        SR_FragLenHistArrayWrite(pHistArray, histOutput);


        // clear the bam in stream
        SR_BamInStreamClear(pBamInStream);

        // rewind the bam file to the beginning of the alignments (right after the header)
        SR_BamInStreamSeek(pBamInStream, bamPos, SEEK_SET);

        // we only have to create the read pair table and open the read pair files once
        if (!hasReadPairTable)
        {
            pReadPairTable = SR_ReadPairTableAlloc(pLibTable->pAnchorInfo->size, pBuildPars->detectSet); 
            pFileHash = SR_ReadPairFilesOpen(pLibTable, pBuildPars->detectSet, pBuildPars->workingDir);
            hasReadPairTable =TRUE;
        }

        // we need to load the cross pair if we want to detect inter-chromosome translocation
        if ((pBuildPars->detectSet & SV_INTER_CHR_TRNSLCTN) != 0)
            SR_FilterDataRPTurnOnCross(pFilterData);

        while ((bamStatus = SR_BamInStreamLoadPair(&pUpNode, &pDownNode, pBamInStream)) != SR_EOF && bamStatus != SR_ERR)
        {
            // we hit another chromosome
            if (bamStatus == SR_OUT_OF_RANGE)
                continue;

            // load the read pairs
            if (!pFilterData->isFilled)
            {
                // usually we should catch the read pairs here
                // since their fragment length should be smaller than the bin length
                pUpAlgn = &(pUpNode->alignment);
                pDownAlgn = &(pDownNode->alignment);
            }
            else
            {
                // if the fragment length of the read pair is beyond our bin length
                // we should catch them here
                pUpAlgn = pFilterData->pUpAlgn;
                pDownAlgn = pFilterData->pDownAlgn;
            }

            SR_ZAtag zaTag;
            SR_PairStats pairStats;

            SR_Status readStatus = SR_LoadPairStats(&pairStats, pUpAlgn, pLibTable);

            if (readStatus == SR_OK)
            {
                SR_Status ZAstatus = SR_LoadZAtag(&zaTag, pUpAlgn);

                if (ZAstatus == SR_OK)
                    SR_ReadPairTableUpdate(pReadPairTable, pUpAlgn, pDownAlgn, &zaTag, &pairStats, pLibTable, pHistArray, pBuildPars->minMQ);
                else
                    SR_ReadPairTableUpdate(pReadPairTable, pUpAlgn, pDownAlgn, NULL, &pairStats, pLibTable, pHistArray, pBuildPars->minMQ);
            }


            // recycle those bam nodes that are allocated from the memory pool
            if (!pFilterData->isFilled)
            {
                SR_BamInStreamRecycle(pBamInStream, pUpNode);
                SR_BamInStreamRecycle(pBamInStream, pDownNode);
            }
        }

        SR_ReadPairTableWrite(pReadPairTable, pFileHash);
        SR_ReadPairTableClear(pReadPairTable);

        // close the bam file
        SR_BamInStreamClose(pBamInStream);
        SR_BamHeaderFree(pBamHeader);
    }

    // close all the read pair files
    SR_ReadPairFilesClose(pFileHash);

    // get the library table output file name
    char* libTableOutputFile = SR_CreateFileName(pBuildPars->workingDir, SR_LibTableFileName);

    // write the library information table into the file
    FILE* libTableOutput = fopen(libTableOutputFile, "wb");
    if (libTableOutput == NULL)
        SR_ErrQuit("ERROR: Cannot open the library output file: %s\n", libTableOutputFile);

    SR_LibInfoTableWrite(pLibTable, libTableOutput);
    SR_SpecialPairTableWirteID(pReadPairTable->pSpecialPairTable, libTableOutput);

    // clean up
    fclose(histOutput);
    fclose(libTableOutput);

    free(libTableOutputFile);
    free(histOutputFile);

    SR_FragLenHistArrayFree(pHistArray);
    SR_LibInfoTableFree(pLibTable);
    SR_FilterDataRPFree(pFilterData);
    SR_ReadPairTableFree(pReadPairTable);
    SR_BamInStreamFree(pBamInStream);
}

/* 
void SR_ReadPairBuild(const SR_ReadPairBuildPars* pBuildPars)
{
    // some default capacity of the containers
    unsigned int buffCapacity = 100;
    unsigned int capAnchor = 150;
    unsigned int capSample = 20;
    unsigned int capReadGrp = 20;
    unsigned int capHist = 10;

    // phony parameters for the bam input stream
    unsigned int reportSize = 0;
    unsigned int numThread = 0;

    // initialize the library information table
    SR_LibInfoTable* pLibTable = SR_LibInfoTableAlloc(capAnchor, capSample, capReadGrp);
    SR_LibInfoTableSetCutoff(pLibTable, pBuildPars->cutoff);
    SR_LibInfoTableSetTrimRate(pLibTable, pBuildPars->trimRate);

    // this is the data used to filter the read pairs in the bam
    SR_FilterDataRP* pFilterData = SR_FilterDataRPAlloc(pLibTable->pAnchorInfo, pBuildPars->binLen);

    // set the stream mode to read pair filter
    SR_StreamMode streamMode;
    SR_SetStreamMode(&streamMode, SR_ReadPairFilter, pFilterData, SR_READ_PAIR_MODE);

    // structure initialization
    SR_BamInStream* pBamInStream = SR_BamInStreamAlloc(pBuildPars->binLen, numThread, buffCapacity, reportSize, &streamMode);
    SR_FragLenHistArray* pHistArray = SR_FragLenHistArrayAlloc(capHist);
    SR_BamHeader* pBamHeader = NULL;

    // a flag to indicate if we have created the read pair table
    SR_Bool hasReadPairTable = FALSE;
    SR_ReadPairTable* pReadPairTable = NULL;

    // file hash
    khash_t(file)* pFileHash = NULL;

    // required arguments for bam in stream structure
    char* histOutputFile = SR_CreateFileName(pBuildPars->workingDir, SR_HistFileName);

    // open the fragment length histogram output file
    FILE* histOutput = fopen(histOutputFile, "w");
    if (histOutput == NULL)
        SR_ErrQuit("ERROR: Cannot open fragment length histogram file: %s\n", histOutputFile);

    char bamFileName[1024];
    while (SR_GetNextLine(bamFileName, 1024, pBuildPars->fileListInput) == SR_OK)
    {
        // update the fitler data with new bam file
        SR_FilterDataRPInit(pFilterData, bamFileName);

        // open the bam file
        SR_BamInStreamOpen(pBamInStream, bamFileName);

        // load the bam header before read any alignments
        pBamHeader = SR_BamInStreamLoadHeader(pBamInStream);

        // get the file position of the bam aligments
        int64_t bamPos = SR_BamInStreamTell(pBamInStream);

        // process the header information
        unsigned int oldSize = 0;
        SR_LibInfoTableSetRG(pLibTable, &oldSize, pBamHeader);

        // initialize the fragment length histogram array with the number of newly added libraries in the bam file
        SR_FragLenHistArrayInit(pHistArray, pLibTable->size - oldSize);

        SR_BamNode* pUpNode = NULL;
        SR_BamNode* pDownNode = NULL;

        const bam1_t* pUpAlgn = NULL;
        const bam1_t* pDownAlgn = NULL;

        SR_Status bamStatus = SR_OK;
        while ((bamStatus = SR_BamInStreamLoadPair(&pUpNode, &pDownNode, pBamInStream)) != SR_EOF && bamStatus != SR_ERR)
        {
            // we hit another chromosome
            if (bamStatus == SR_OUT_OF_RANGE)
                continue;

            // load the read pairs
            if (!pFilterData->isFilled)
            {
                // usually we should catch the read pairs here
                // since their fragment length should be smaller than the bin length
                pUpAlgn = &(pUpNode->alignment);
                pDownAlgn = &(pDownNode->alignment);
            }
            else
            {
                // if the fragment length of the read pair is beyond our bin length
                // we should catch them here
                pUpAlgn = pFilterData->pUpAlgn;
                pDownAlgn = pFilterData->pDownAlgn;
            }

            // check if the incoming read pair is normal (unique-unique pair)
            // if yes, then update the corresponding fragment length histogram
            SR_PairStats pairStats;
            unsigned int backHistIndex = 0;
            if (SR_IsNormalPair(&pairStats, &backHistIndex, pUpAlgn, pDownAlgn, pLibTable, pBuildPars->minMQ))
                SR_FragLenHistArrayUpdate(pHistArray, backHistIndex, pairStats.fragLen);

            // recycle those bam nodes that are allocated from the memory pool
            if (!pFilterData->isFilled)
            {
                SR_BamInStreamRecycle(pBamInStream, pUpNode);
                SR_BamInStreamRecycle(pBamInStream, pDownNode);
            }
        }

        // finish the process of the histogram and update the library information table
        SR_FragLenHistArrayFinalize(pHistArray);
        SR_LibInfoTableUpdate(pLibTable, pHistArray, oldSize);

        // write the fragment length histogram into the file
        SR_FragLenHistArrayWrite(pHistArray, histOutput);


        // rewind the bam file to the beginning of the alignments (right after the header)
        SR_BamInStreamSeek(pBamInStream, bamPos, SEEK_SET);

        // we only have to create the read pair table and open the read pair files once
        if (!hasReadPairTable)
        {
            pReadPairTable = SR_ReadPairTableAlloc(pLibTable->pAnchorInfo->size, pBuildPars->detectSet); 
            pFileHash = SR_ReadPairFilesOpen(pLibTable, pBuildPars->detectSet, pBuildPars->workingDir);
            hasReadPairTable =TRUE;
        }

        while ((bamStatus = SR_BamInStreamLoadPair(&pUpNode, &pDownNode, pBamInStream)) != SR_EOF && bamStatus != SR_ERR)
        {
            // we hit another chromosome
            if (bamStatus == SR_OUT_OF_RANGE)
                continue;

            // load the read pairs
            if (!pFilterData->isFilled)
            {
                // usually we should catch the read pairs here
                // since their fragment length should be smaller than the bin length
                pUpAlgn = &(pUpNode->alignment);
                pDownAlgn = &(pDownNode->alignment);
            }
            else
            {
                // if the fragment length of the read pair is beyond our bin length
                // we should catch them here
                pUpAlgn = pFilterData->pUpAlgn;
                pDownAlgn = pFilterData->pDownAlgn;
            }

            SR_ZAtag zaTag;
            SR_PairStats pairStats;

            SR_Status readStatus = SR_LoadPairStats(&pairStats, pUpAlgn, pLibTable);
            if (readStatus == SR_OK)
            {
                SR_Status ZAstatus = SR_LoadZAtag(&zaTag, pUpAlgn);
                if (ZAstatus == SR_OK)
                    SR_ReadPairTableUpdate(pReadPairTable, pUpAlgn, pDownAlgn, &zaTag, &pairStats, pLibTable, pHistArray, pBuildPars->minMQ);
                else
                    SR_ReadPairTableUpdate(pReadPairTable, pUpAlgn, pDownAlgn, NULL, &pairStats, pLibTable, pHistArray, pBuildPars->minMQ);
            }

            // recycle those bam nodes that are allocated from the memory pool
            if (!pFilterData->isFilled)
            {
                SR_BamInStreamRecycle(pBamInStream, pUpNode);
                SR_BamInStreamRecycle(pBamInStream, pDownNode);
            }
        }

        SR_ReadPairTableWrite(pReadPairTable, pFileHash);
        SR_ReadPairTableClear(pReadPairTable);

        // close the bam file
        SR_BamInStreamClose(pBamInStream);
        SR_BamHeaderFree(pBamHeader);
    }

    // close all the read pair files
    SR_ReadPairFilesClose(pFileHash);

    // get the library table output file name
    char* libTableOutputFile = SR_CreateFileName(pBuildPars->workingDir, SR_LibTableFileName);

    // write the library information table into the file
    FILE* libTableOutput = fopen(libTableOutputFile, "wb");
    if (libTableOutput == NULL)
        SR_ErrQuit("ERROR: Cannot open the library output file: %s\n", libTableOutputFile);

    SR_LibInfoTableWrite(pLibTable, libTableOutput);
    SR_SpecialPairTableWirteID(pReadPairTable->pSpecialPairTable, libTableOutput);

    // clean up
    fclose(histOutput);
    fclose(libTableOutput);

    free(libTableOutputFile);
    free(histOutputFile);

    SR_FragLenHistArrayFree(pHistArray);
    SR_LibInfoTableFree(pLibTable);
    SR_FilterDataRPFree(pFilterData);
    SR_ReadPairTableFree(pReadPairTable);
    SR_BamInStreamFree(pBamInStream);
}
*/

void SR_ReadPairTableClear(SR_ReadPairTable* pReadPairTable)
{
    if (pReadPairTable->pLongPairArray != NULL)
    {
        memset(pReadPairTable->pLongPairArray->chrCount, 0, sizeof(uint64_t) * pReadPairTable->numChr);
        pReadPairTable->pLongPairArray->size = 0;
    }

    if (pReadPairTable->pShortPairArray != NULL)
    {
        memset(pReadPairTable->pShortPairArray->chrCount, 0, sizeof(uint64_t) * pReadPairTable->numChr);
        pReadPairTable->pShortPairArray->size = 0;
    }

    if (pReadPairTable->pReversedPairArray != NULL)
    {
        memset(pReadPairTable->pReversedPairArray->chrCount, 0, sizeof(uint64_t) * pReadPairTable->numChr);
        pReadPairTable->pReversedPairArray->size = 0;
    }

    if (pReadPairTable->pInvertedPairArray != NULL)
    {
        memset(pReadPairTable->pInvertedPairArray->chrCount, 0, sizeof(uint64_t) * pReadPairTable->numChr);
        pReadPairTable->pInvertedPairArray->size = 0;
    }

    if (pReadPairTable->pCrossPairArray != NULL)
    {
        memset(pReadPairTable->pCrossPairArray->chrCount, 0, sizeof(uint64_t) * pReadPairTable->numChr);
        pReadPairTable->pCrossPairArray->size = 0;
    }

    SR_SpecialPairTableClear(pReadPairTable->pSpecialPairTable, pReadPairTable->numChr);
}

void SR_ReadPairTableUpdate(SR_ReadPairTable* pReadPairTable, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const SR_ZAtag* pZAtag,
        const SR_PairStats* pPairStats, const SR_LibInfoTable* pLibTable, const SR_FragLenHistArray* pHistArray, uint8_t minMQ)
{
    SV_ReadPairType readPairType = PT_NORMAL;
    readPairType = SR_CheckReadPairType(pZAtag, pPairStats, pLibTable);

    if (readPairType != PT_SPECIAL5 && readPairType != PT_SPECIAL3)
    {
        if (pUpAlgn->core.qual < minMQ || pDownAlgn->core.qual < minMQ)
            readPairType = PT_UNKNOWN;
    }

    switch (readPairType)
    {
        case PT_NORMAL:
        case PT_UNKNOWN:
            break;
        case PT_LONG:
            if (pReadPairTable->pLongPairArray != NULL)
                SR_LocalPairArrayUpdate(pReadPairTable->pLongPairArray, pUpAlgn, pDownAlgn, pLibTable, pHistArray, pPairStats, readPairType);
            break;
        case PT_SHORT:
            if (pReadPairTable->pShortPairArray != NULL)
                SR_LocalPairArrayUpdate(pReadPairTable->pShortPairArray, pUpAlgn, pDownAlgn, pLibTable, pHistArray, pPairStats, readPairType);
            break;
        case PT_REVERSED:
            if (pReadPairTable->pReversedPairArray != NULL)
                SR_LocalPairArrayUpdate(pReadPairTable->pReversedPairArray, pUpAlgn, pDownAlgn, pLibTable, pHistArray, pPairStats, readPairType);
            break;
        case PT_INVERTED3:
        case PT_INVERTED5:
            if (pReadPairTable->pInvertedPairArray != NULL)
                SR_InvertedPairArrayUpdate(pReadPairTable->pInvertedPairArray, pUpAlgn, pDownAlgn, pPairStats, readPairType);
            break;
        case PT_CROSS:
            if (pReadPairTable->pCrossPairArray != NULL)
                SR_CrossPairArrayUpdate(pReadPairTable->pCrossPairArray, pUpAlgn, pDownAlgn, pPairStats);
            break;
        case PT_SPECIAL3:
        case PT_SPECIAL5:
            if (pReadPairTable->pSpecialPairTable != NULL) 
                SR_SpecialPairTableUpdate(pReadPairTable->pSpecialPairTable, pUpAlgn, pDownAlgn, pLibTable, pHistArray, pPairStats, pZAtag, readPairType);
            break;
        default:
            break;
    }
}

void* SR_ReadPairFilesOpen(const SR_LibInfoTable* pLibTable, uint32_t detectSet, const char* workingDir)
{
    khash_t(file)* pFileHash = kh_init(file);

    char refBuff[4];
    int dirLen = strlen(workingDir);
    char nameBuff[dirLen + 50];

    strncpy(nameBuff, workingDir, dirLen);
    char* namePos = nameBuff + dirLen;
    char* replacePos = namePos + 3;

    const SR_AnchorInfo* pAnchorInfo = pLibTable->pAnchorInfo;
    if (pAnchorInfo->size > 1000)
        SR_ErrQuit("ERROR: reference ID overflow (0-999).\n");

    // this is the place holder put at the beginning of each read pair file
    // different SV types may have different file headers
    const uint64_t placeHolder = 0;
    SR_ReadPairOutStream stream = {NULL, 0};

    for (unsigned int i = 0; i != pAnchorInfo->size; ++i)
    {
        if (pAnchorInfo->pLength[i] <= 0)
            continue;

        unsigned int offset = 2 - (unsigned int) log10(i);
        sprintf(refBuff, "000");
        sprintf(refBuff + offset, "%d", i);

        int fileID = 0;
        int ret = 0;
        khiter_t khIter = 0;
        unsigned int fileIndex = SR_LONG_PAIR_FILE;


        for (unsigned int j = SV_DELETION; j <= SV_SPECIAL; ++j)
        {
            if ((detectSet & (1 << j)) != 0)
            {
                strcpy(namePos, SR_READ_PAIR_FILE_NAME_TEMPLATE[fileIndex]);
                strncpy(replacePos, refBuff, 3);

                stream.output = fopen(nameBuff, "wb");
                if (stream.output == NULL)
                    SR_ErrQuit("ERROR: Cannot open file: %d\n", nameBuff);

                fileID = i * SR_NUM_RP_FILE + fileIndex; 
                khIter = kh_put(file, pFileHash, fileID, &ret);
                kh_value(pFileHash, khIter) = stream;
                // leave a place for the total number of read pairs
                fwrite(&placeHolder, sizeof(uint64_t), 1, stream.output);

                if (j == SV_TANDEM_DUP)
                {
                    ++fileIndex;

                    strcpy(namePos, SR_READ_PAIR_FILE_NAME_TEMPLATE[fileIndex]);
                    strncpy(replacePos, refBuff, 3);

                    stream.output = fopen(nameBuff, "wb");
                    if (stream.output == NULL)
                        SR_ErrQuit("ERROR: Cannot open file: %d\n", nameBuff);

                    fileID = i * SR_NUM_RP_FILE + fileIndex; 
                    khIter = kh_put(file, pFileHash, fileID, &ret);
                    kh_value(pFileHash, khIter) = stream;

                    // leave a place for the total number of read pairs
                    fwrite(&placeHolder, sizeof(uint64_t), 1, stream.output);
                }
            }

            ++fileIndex;
        }
    }

    return ((void*) pFileHash);
}

void SR_ReadPairFilesClose(void* pFileHash)
{
    khash_t(file)* pHash = (khash_t(file)*) pFileHash;

    for (khiter_t khIter = kh_begin(pHash); khIter != kh_end(pHash); ++khIter)
    {
        if (kh_exist(pHash, khIter))
        {
            fseeko(kh_value(pHash, khIter).output, 0, SEEK_SET);
            fwrite(&(kh_value(pHash, khIter).numPairs), sizeof(uint64_t), 1, kh_value(pHash, khIter).output);
            fclose(kh_value(pHash, khIter).output);
        }
    }

    kh_destroy(file, pHash);
}

void SR_ReadPairTableWrite(const SR_ReadPairTable* pReadPairTable, void* pHash)
{
    khash_t(file)* pFileHash = (khash_t(file)*) pHash;
    unsigned int numChr = pReadPairTable->numChr;

    int fileID = 0;
    khiter_t khIter = 0;

    int currPos[] = {0, 0, 0, 0, 0, 0};

    for (unsigned int i = 0; i != numChr; ++i)
    {
        int numPairs = pReadPairTable->pLongPairArray == NULL ? 0 : pReadPairTable->pLongPairArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * SR_NUM_RP_FILE + SR_LONG_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pLongPairArray->data + currPos[SR_LONG_PAIR_FILE], sizeof(SR_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[SR_LONG_PAIR_FILE] += numPairs;
        }

        numPairs = pReadPairTable->pShortPairArray == NULL ? 0 : pReadPairTable->pShortPairArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * SR_NUM_RP_FILE + SR_SHORT_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pShortPairArray->data + currPos[SR_SHORT_PAIR_FILE], sizeof(SR_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[SR_SHORT_PAIR_FILE] += numPairs;
        }

        numPairs = pReadPairTable->pReversedPairArray == NULL ? 0 : pReadPairTable->pReversedPairArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * SR_NUM_RP_FILE + SR_REVERSED_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pReversedPairArray->data + currPos[SR_REVERSED_PAIR_FILE], sizeof(SR_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[SR_REVERSED_PAIR_FILE] += numPairs;
        }

        numPairs = pReadPairTable->pInvertedPairArray == NULL ? 0 : pReadPairTable->pInvertedPairArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * SR_NUM_RP_FILE + SR_INVERTED_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pInvertedPairArray->data + currPos[SR_INVERTED_PAIR_FILE], sizeof(SR_LocalPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[SR_INVERTED_PAIR_FILE] += numPairs;
        }

        numPairs = pReadPairTable->pCrossPairArray == NULL ? 0 : pReadPairTable->pCrossPairArray->chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * SR_NUM_RP_FILE + SR_CROSS_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pCrossPairArray->data + currPos[SR_CROSS_PAIR_FILE], sizeof(SR_CrossPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[SR_CROSS_PAIR_FILE] += numPairs;
        }

        numPairs = pReadPairTable->pSpecialPairTable == NULL ? 0 : pReadPairTable->pSpecialPairTable->array.chrCount[i];
        if (numPairs > 0)
        {
            fileID = i * SR_NUM_RP_FILE + SR_SPEICAL_PAIR_FILE;
            khIter = kh_get(file, pFileHash, fileID);

            if (khIter != kh_end(pFileHash))
            {
                kh_value(pFileHash, khIter).numPairs += numPairs;
                fwrite(pReadPairTable->pSpecialPairTable->array.data + currPos[SR_SPEICAL_PAIR_FILE], sizeof(SR_SpecialPair), numPairs, kh_value(pFileHash, khIter).output);
            }

            currPos[SR_SPEICAL_PAIR_FILE] += numPairs;
        }
    }

    // take care of those special pairs whose anchor refID is not the same as its mate's refID.
    const SR_SpecialPairArray* pSpecialPairArray = &(pReadPairTable->pSpecialPairTable->crossArray);
    for (unsigned int i = 0; i != pSpecialPairArray->size; ++i)
    {
        fileID = pSpecialPairArray->data[i].refID[0] * SR_NUM_RP_FILE + SR_SPEICAL_PAIR_FILE;
        khIter = kh_get(file, pFileHash, fileID);

        if (khIter != kh_end(pFileHash))
        {
            kh_value(pFileHash, khIter).numPairs += 1;
            fwrite(pSpecialPairArray->data + i, sizeof(SR_SpecialPair), 1, kh_value(pFileHash, khIter).output);
        }
    }
}

void SR_SpecialPairTableWirteID(const SR_SpecialPairTable* pSpecialTable, FILE* libOutput)
{
    // write the name of special reference at the end of the library file
    fwrite(&(pSpecialTable->size), sizeof(uint32_t), 1, libOutput);
    for (unsigned int i = 0; i != pSpecialTable->size; ++i)
        fwrite(pSpecialTable->names + i, sizeof(char), 2, libOutput);

    fflush(libOutput);
}

static SR_Status SR_DoDir(const char* path, mode_t mode)
{
    struct stat     st;
    SR_Status status = SR_OK;

    if (stat(path, &st) != 0)
    {
        /*  Directory does not exist */
        if (mkdir(path, mode) != 0)
            status = SR_ERR;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        status = SR_ERR;
    }

    return status;
}

/* *
 * ** mkpath - ensure all directories in path exist
 * ** Algorithm takes the pessimistic view and works top-down to ensure
 * ** each directory in path exists, rather than optimistically creating
 * ** the last element and working backwards.
 * */
SR_Status SR_MakePath(const char* path, mode_t mode)
{
    char           *pp;
    char           *sp;
    SR_Status      status;
    char           *copypath = strdup(path);

    status = 0;
    pp = copypath;

    while (status == SR_OK && (sp = strchr(pp, '/')) != 0)
    {
        if (sp != pp)
        {
            /*  Neither root nor double slash in path */
            *sp = '\0';
            status = SR_DoDir(copypath, mode);
            *sp = '/';
        }

        pp = sp + 1;
    }

    if (status == 0)
        status = SR_DoDir(path, mode);

    free(copypath);
    return status;
}

// make sure the working dir is empty
SR_Status SR_CheckWorkingDir(const char* workingDir)
{
    SR_Status status = SR_MakePath(workingDir, (S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH));

    if (status == SR_OK)
    {
        DIR* pDir = opendir(workingDir);
        if (pDir != NULL)
        {
            status = SR_ERR;
            struct dirent* pDirRecord = NULL;

            while ((pDirRecord = readdir(pDir)) != NULL)
            {
                status = SR_OK;
                if (strcmp(pDirRecord->d_name, ".") != 0 && strcmp(pDirRecord->d_name, "..") != 0)
                {
                    status = SR_ERR;
                    SR_ErrMsg("ERROR: Working directory \"%s\" is not empty.\n", workingDir);
                    break;
                }
            }

            closedir(pDir);
        }
        else
        {
            status = SR_ERR;
            SR_ErrMsg("ERROR: Cannot open working directory: %s.\n", workingDir);
        }
    }

    return status;
}
