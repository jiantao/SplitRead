/*
 * =====================================================================================
 *
 *       Filename:  SR_ReadPairDetect.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/23/2012 08:14:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <limits.h>

#include "khash.h"
#include "SR_Error.h"
#include "SR_Utilities.h"
#include "SR_ReadPairDetect.h"

static const char* SR_LibTableFileName = "lib_table.dat";

// static const char* SR_HistFileName = "hist.dat";

/*  
static const char* SR_READ_PAIR_FILE_NAME_TEMPLATE[] = 
{
    "refXXX_long_pairs.dat",

    "refXXX_short_pairs.dat",

    "refXXX_reversed_pairs.dat",

    "refXXX_inverted_pairs.dat",

    "refXXX_cross_pairs.dat",

    "refXXX_special_pairs.dat",
};
*/

KHASH_MAP_INIT_STR(name, uint32_t);

static void SR_DelEventMerge(SR_DelArray* pDelArray, SR_Cluster* pDelCluster)
{
    SR_DelEvent* pLastEvent = pDelArray->data + (pDelArray->size - 1);
    SR_DelEvent* pNewEvent = pDelArray->data + pDelArray->size;

    pLastEvent->pos5[0] = (pLastEvent->pos5[0] < pNewEvent->pos5[0] ? pLastEvent->pos5[0] : pNewEvent->pos5[0]);
    pLastEvent->pos5[1] = (pLastEvent->pos5[1] > pNewEvent->pos5[1] ? pLastEvent->pos5[1] : pNewEvent->pos5[1]);
    pLastEvent->pos5[2] = (pLastEvent->pos5[2] > pNewEvent->pos5[2] ? pLastEvent->pos5[2] : pNewEvent->pos5[2]);

    pLastEvent->pos3[0] = (pLastEvent->pos3[0] < pNewEvent->pos3[0] ? pLastEvent->pos3[0] : pNewEvent->pos3[0]);
    pLastEvent->pos3[1] = (pLastEvent->pos3[1] > pNewEvent->pos3[1] ? pLastEvent->pos3[1] : pNewEvent->pos3[1]);
    pLastEvent->pos3[2] = (pLastEvent->pos3[2] > pNewEvent->pos3[2] ? pLastEvent->pos3[2] : pNewEvent->pos3[2]);

    int numReadPair = SR_ClusterMerge(pDelCluster, pLastEvent->clusterID, pNewEvent->clusterID);

    pLastEvent->pos = pLastEvent->pos5[1] + 1;
    pLastEvent->end = pLastEvent->pos3[0];
    pLastEvent->length = DoubleRoundToInt(pDelCluster->pElmntArray->data[pLastEvent->clusterID].mean[1]);
    pLastEvent->quality = (int) ((numReadPair * 100.0) / (numReadPair + 10.0));

    pLastEvent->CIpos[0] = -(pLastEvent->pos5[2] - pLastEvent->pos5[0]) / numReadPair;
    pLastEvent->CIpos[1] = (pLastEvent->pos3[2] - pLastEvent->pos3[0]) / numReadPair;

    pLastEvent->CIend[0] = -(pLastEvent->pos5[2] - pLastEvent->pos5[0]) / numReadPair;
    pLastEvent->CIend[1] = (pLastEvent->pos3[2] - pLastEvent->pos3[0]) / numReadPair;

    pLastEvent->CIlen[0] = 0;
    pLastEvent->CIlen[1] = 0;
}

static int CompareDelEvents(const void* pEvent1, const void* pEvent2)
{
    const SR_DelEvent* pE1 = pEvent1;
    const SR_DelEvent* pE2 = pEvent2;

    if (pE1->refID < pE2->refID)
        return -1;
    else if (pE1->refID > pE2->refID)
        return 1;
    else
    {
        if (pE1->pos < pE2->pos)
            return -1;
        else if (pE1->pos > pE2->pos)
            return 1;
        else
        {
            if (pE1->length < pE2->length)
                return -1;
            else if (pE1->length > pE2->length)
                return 1;
            else
            {
                if (pE1->quality < pE2->quality)
                    return -1;
                else if (pE1->quality > pE2->quality)
                    return 1;
            }
        }
    }

    return 0;
}


void SR_ReadPairDetect(const SR_ReadPairDetectPars* pDetectPars)
{
    // create a buffer to store the file name 
    int dirLen = strlen(pDetectPars->workingDir);
    char nameBuff[dirLen + 60];

    // create the library information file name
    strncpy(nameBuff, pDetectPars->workingDir, dirLen);
    strcpy(nameBuff + dirLen, SR_LibTableFileName);

    // read the library information into memory
    FILE* pLibInput = fopen(nameBuff, "rb");
    if (pLibInput == NULL)
        SR_ErrQuit("ERROR: Cannot open library information file: \"%s\".\n", nameBuff);

    // SR_LibInfoTable* pLibTable = SR_LibInfoTableRead(pLibInput);

    // read the detect set
    uint32_t detectSet = 0;
    fread(&detectSet, sizeof(uint32_t), 1, pLibInput);
}


void SR_LocalPairArrayRead(SR_LocalPairArray* pLocalPairArray, FILE* input)
{
    int64_t size = 0;
    fread(&(size), sizeof(int64_t), 1, input);

    if (size > pLocalPairArray->capacity)
        SR_ARRAY_RESIZE_NO_COPY(pLocalPairArray, size, SR_LocalPair);
    else if (size <= 0)
        return;

    pLocalPairArray->size = size;
    fread(pLocalPairArray->data, sizeof(SR_LocalPair), size, input);
}

void SR_CrossPairArrayRead(SR_CrossPairArray* pCrossPairArray, FILE* input)
{
    int64_t size = 0;
    fread(&(size), sizeof(int64_t), 1, input);

    if (size > pCrossPairArray->capacity)
        SR_ARRAY_RESIZE_NO_COPY(pCrossPairArray, size, SR_CrossPair);
    else if (size <= 0)
        return;

    pCrossPairArray->size = size;
    fread(pCrossPairArray->data, sizeof(SR_LocalPair), size, input);
}

void SR_SpecialPairArrayRead(SR_SpecialPairArray* pSpecialPairArray, FILE* input)
{
    int64_t size = 0;
    fread(&(size), sizeof(int64_t), 1, input);

    if (size > pSpecialPairArray->capacity)
        SR_ARRAY_RESIZE_NO_COPY(pSpecialPairArray, size, SR_SpecialPair);
    else if (size <= 0)
        return;

    pSpecialPairArray->size = size;
    fread(pSpecialPairArray->data, sizeof(SR_LocalPair), size, input);
}

void SR_SpecialPairTableReadID(SR_SpecialPairTable* pSpeicalPairTable, FILE* libInput)
{
    fread(&(pSpeicalPairTable->size), sizeof(uint32_t), 1, libInput);

    if (pSpeicalPairTable->size > pSpeicalPairTable->capacity)
    {
        free(pSpeicalPairTable->names);
        pSpeicalPairTable->names = (char (*)[3]) malloc(sizeof(char) * 3 * pSpeicalPairTable->size);
        pSpeicalPairTable->capacity = pSpeicalPairTable->size;
    }
    else if (pSpeicalPairTable->size == 0)
        return;

    int ret = 0;
    khiter_t khIter = 0;
    for (unsigned int i = 0; i != pSpeicalPairTable->size; ++i)
    {
        fread(pSpeicalPairTable->names[i], sizeof(char), 2, libInput);
        pSpeicalPairTable->names[i][2] = '\0';

        khash_t(name)* pHash = pSpeicalPairTable->nameHash;
        khIter = kh_put(name, pHash, pSpeicalPairTable->names[i], &ret);
        kh_value(pHash, khIter) = i;
    }
}

void SR_ReadPairFindDel(SR_DelArray* pDelArray, SV_AssistArray* pAssistArray, const SR_LocalPairArray* pLongPairArray,
                        SR_Cluster* pDelCluster, const SR_LibInfoTable* pLibTable, const SR_ReadPairDetectPars* pPars)
{
    SR_ARRAY_RESIZE_NO_COPY(pDelArray, pDelCluster->pElmntArray->size, SR_DelEvent);

    unsigned int i = 0;
    unsigned int counter = 0;
    unsigned int numClusters = pDelCluster->pElmntArray->size;

    while (counter != numClusters)
    {
        if (pDelCluster->pElmntArray->data[i].numReadPair == 0)
        {
            ++i;
            continue;
        }

        int posMin5 = INT_MAX;
        int posMax5 = 0;

        int posMin3 = INT_MAX;
        int posMax3 = 0;

        int endMax5 = 0;
        int endMax3 = 0;

        int fragLenDiffMin = INT_MAX;
        int fragLenDiffMax = 0;

        int fragLenMax = 0;

        int numReadPair = pDelCluster->pElmntArray->data[i].numReadPair;
        // SV_AssistArrayResize(pAssistArray, numReadPair);

        unsigned int j = pDelCluster->pElmntArray->data[i].startIndex;
        unsigned int startIndex = pDelCluster->pElmntArray->data[i].startIndex;
        const SR_LocalPair* pLongPair = NULL;

        do
        {
            unsigned int origIndex = pDelCluster->pAttrbtArray->data[j].origIndex;
            pLongPair = (pLongPairArray->data + origIndex);

            if (pLongPair->upPos < posMin5)
                posMin5 = pLongPair->upPos;

            if (pLongPair->upPos > posMax5)
                posMax5 = pLongPair->upPos;

            if (pLongPair->upEnd > endMax5)
                endMax5 = pLongPair->upEnd;

            if (pLongPair->downPos < posMin3)
                posMin3 = pLongPair->downPos;

            if (pLongPair->downPos > posMax3)
                posMax3 = pLongPair->downPos;

            unsigned int downEnd = pLongPair->upPos + pLongPair->fragLen;
            if (downEnd > endMax3)
                endMax3 = downEnd;

            unsigned int upperFragLen = pLibTable->pLibInfo[pLongPair->readGrpID].fragLenHigh;
            unsigned int medianFragLen = pLibTable->pLibInfo[pLongPair->readGrpID].fragLenMedian;

            if (upperFragLen > fragLenMax)
                fragLenMax = upperFragLen;

            pAssistArray->pFragLenDiff[j] = pLongPair->fragLen - medianFragLen;
            pAssistArray->pMapQ5[j] = pLongPair->upMapQ;
            pAssistArray->pMapQ3[j] = pLongPair->downMapQ;

            if (pAssistArray->pFragLenDiff[j] < fragLenDiffMin)
                fragLenDiffMin = pAssistArray->pFragLenDiff[j];

            if (pAssistArray->pFragLenDiff[j] > fragLenDiffMax)
                fragLenDiffMax = pAssistArray->pFragLenDiff[j];

            j = pDelCluster->pNext[j];

        }while(j != startIndex);

        int eventLength = FindMedianInt(pAssistArray->pFragLenDiff, pAssistArray->size);
        if (eventLength < 1)
            continue;

        int k = pDelArray->size;

        pDelArray->data[k].clusterID = i;

        pDelArray->data[k].refID = pLongPair->refID;
        pDelArray->data[k].pos = endMax5 + 1;
        pDelArray->data[k].end = posMin3;
        pDelArray->data[k].length = FindMedianInt(pAssistArray->pFragLenDiff, pAssistArray->size);
        pDelArray->data[k].quality = (int) ((numReadPair * 100.0) / (numReadPair + 10.0));

        pDelArray->data[k].pos5[0] = posMin5;
        pDelArray->data[k].pos5[1] = endMax5;
        pDelArray->data[k].pos5[2] = posMax5;

        pDelArray->data[k].pos3[0] = posMin3;
        pDelArray->data[k].pos3[1] = endMax3;
        pDelArray->data[k].pos3[2] = posMax3;

        pDelArray->data[k].CIpos[0] = -(posMax5 - posMin5) / numReadPair;
        pDelArray->data[k].CIpos[1] = (posMax3 - posMin3) / numReadPair;

        pDelArray->data[k].CIend[0] = -(posMax5 - posMin5) / numReadPair;
        pDelArray->data[k].CIend[1] = (posMax3 - posMin3) / numReadPair;

        pDelArray->data[k].CIlen[0] = -(eventLength - fragLenDiffMin) / numReadPair;
        pDelArray->data[k].CIlen[1] = (fragLenDiffMax - eventLength) / numReadPair;

        pDelArray->data[k].mapQ5 = FindMedianInt(pAssistArray->pMapQ5, pAssistArray->size);
        pDelArray->data[k].mapQ3 = FindMedianInt(pAssistArray->pMapQ3, pAssistArray->size);

        SR_Bool isOverlap = FALSE;
        if (pDelArray->size > 0)
        {
            int lastIndex = pDelArray->size - 1;

            if (pDelArray->data[lastIndex].end > pDelArray->data[k].pos 
                && abs(pDelArray->data[lastIndex].pos - pDelArray->data[k].pos) < fragLenMax)
            {
                isOverlap = TRUE;
            }
        }

        // correction to pos and len for deletions bracketed by repeat region 
        int fragLenDiffError = eventLength - DoubleRoundToInt(pDelCluster->pElmntArray->data[i].mean[1]);
        if (fragLenDiffError > 4 * (pDelArray->data[k].CIpos[1]))
        {
            pDelArray->data[k].pos += (int) (fragLenDiffError / 2);
            pDelArray->data[k].CIpos[0] = (int) (-fragLenDiffError / 2);
            pDelArray->data[k].CIpos[1] = (int) (fragLenDiffError / 2);

            pDelArray->data[k].length = DoubleRoundToInt(pDelCluster->pElmntArray->data[i].mean[1]);
            pDelArray->data[k].CIlen[0] = DoubleRoundToInt(-pDelCluster->pElmntArray->data[i].std[1] / sqrt(numReadPair));
            pDelArray->data[k].CIlen[1] = DoubleRoundToInt(pDelCluster->pElmntArray->data[i].std[1] / sqrt(numReadPair));
        }

        if (isOverlap)
            SR_DelEventMerge(pDelArray, pDelCluster);
        else
        {
            if (numReadPair >= pPars->minNumClustered && pDelArray->data[k].length >= pPars->minEventLength)
                ++(pDelArray->size);
        }

        ++i;
        ++counter;
    }

    qsort(pDelArray->data, pDelArray->size, sizeof(SR_DelEvent), CompareDelEvents);

    //SR_DelEventGenotype(pDelArray, pLibTable);
}
