/*
 * =====================================================================================
 *
 *       Filename:  SR_FragLenDstrb.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/15/2011 03:43:08 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "SR_FragLenDstrb.h"

SR_FragLenDstrb* SR_FragLenDstrbAlloc(uint8_t drctField, uint32_t capacity)
{
    SR_FragLenDstrb* pNewDstrb = (SR_FragLenDstrb*) malloc(sizeof(*pNewDstrb));
    if (pNewDstrb == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the fragment length distribution object.\n");

    if (capacity == 0)
        capacity = DEFAULT_FRAG_DSTRB_CAP;

    pNewDstrb->pReadGrpNames = (char**) calloc(capacity, sizeof(char*));
    if (pNewDstrb->pReadGrpNames == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of read group names in the fragment length distribution object.\n");

    pNewDstrb->pReadGrpHash = kh_init(readGrpName);
    pNewDstrb->pFragLenHists = NULL;

    pNewDstrb->drctField = drctField;
    pNewDstrb->size = 0;
    pNewDstrb->capacity = capacity;

    return pNewDstrb;
}

void SR_FragLenDstrbFree(SR_FragLenDstrb* pDstrb)
{
    if (pDstrb != NULL)
    {
        if (pDstrb->pFragLenHists != NULL)
        {
            for (unsigned int i = 0; i != pDstrb->size; ++i)
                gsl_histogram_free(pDstrb->pFragLenHists[i]);

            free(pDstrb->pFragLenHists);
        }

        if (pDstrb->pReadGrpNames != NULL)
        {
            for (unsigned int i = 0; i != pDstrb->size; ++i)
                free(pDstrb->pReadGrpNames[i]);

            free(pDstrb->pReadGrpNames);
        }

        kh_destroy(readGrpName, pDstrb->pReadGrpHash);

        free(pDstrb);
    }
}

SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader)
{
    for (const char* readGrpPos = pBamHeader->pOrigHeader->text; ; ++readGrpPos)
    {
        if ((readGrpPos = strstr(readGrpPos, "@RG")) != NULL 
            && (readGrpPos = strstr(readGrpPos, "ID:")) != NULL)
        {
            const char* readGrpEnd = strpbrk(readGrpPos, " \t\n\0");
            if (readGrpEnd == NULL)
                return SR_ERR;

            if (pDstrb->size == pDstrb->capacity)
            {
                pDstrb->capacity *= 2;
                pDstrb->pReadGrpNames = (char**) realloc(pDstrb->pReadGrpNames, pDstrb->capacity * sizeof(char*));
                if (pDstrb->pReadGrpNames == NULL)
                    SR_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");
            }

            size_t readGrpLen = readGrpEnd - readGrpPos + 1;
            pDstrb->pReadGrpNames[pDstrb->size] = (char*) calloc(readGrpLen + 1, sizeof(char));
            if (pDstrb->pReadGrpNames[pDstrb->size] == NULL)
                SR_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");

            strncpy(pDstrb->pReadGrpNames[pDstrb->size], readGrpPos, readGrpLen);
            ++(pDstrb->size);
        }
        else
            break;
    }

    if (pDstrb->size == 0)
    {
        SR_FragLenDstrbFree(pDstrb);
        return SR_NOT_FOUND;
    }
    else
    {
        int khRet = 0;
        khiter_t khIter;
        for (unsigned int i = 0; i != pDstrb->size; ++i)
        {
            khIter = kh_put(readGrpName, pDstrb->pReadGrpHash, pDstrb->pReadGrpNames[i], &khRet);
            if (khRet == 0)
            {
                SR_ErrMsg("ERROR: Found a non-unique read group name.\n");
                return SR_ERR;
            }

            kh_value((khash_t(readGrpName)*) pDstrb->pReadGrpHash, khIter) = i;
        }
    }

    return SR_OK;
}

void SR_FragLenDstrbSetHist(SR_FragLenDstrb* pDstrb, size_t numBins)
{
    pDstrb->pFragLenHists = (gsl_histogram**) malloc(pDstrb->size * sizeof(gsl_histogram*));
    if (pDstrb->pFragLenHists == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storageof the histogram in the fragment length distribution object.\n");

    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        pDstrb->pFragLenHists[i] = gsl_histogram_alloc(numBins);
        if (pDstrb->pFragLenHists[i] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storageof the histogram in the fragment length distribution object.\n");
    }
}

