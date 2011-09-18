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

#include "khash.h"
#include "SR_Error.h"
#include "SR_FragLenDstrb.h"

// default capacity of fragment length distribution
#define DEFAULT_FRAG_DSTRB_CAP 10

// default fragment length distribution
#define DEFAULT_FRAG_LEN 10000.0

// read group name hash
KHASH_MAP_INIT_STR(readGrpName, uint32_t);

SR_FragLenDstrb* SR_FragLenDstrbAlloc(unsigned short minMQ, uint32_t capacity)
{
    SR_FragLenDstrb* pNewDstrb = (SR_FragLenDstrb*) calloc(1, sizeof(*pNewDstrb));
    if (pNewDstrb == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the fragment length distribution object.\n");

    if (capacity == 0)
        capacity = DEFAULT_FRAG_DSTRB_CAP;

    pNewDstrb->pReadGrpNames = (char**) calloc(capacity, sizeof(char*));
    if (pNewDstrb->pReadGrpNames == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of read group names in the fragment length distribution object.\n");


    pNewDstrb->pReadGrpHash = NULL;
    pNewDstrb->pFragLenHists = NULL;
    pNewDstrb->pFragLenRange = NULL;

    pNewDstrb->size = 0;
    pNewDstrb->capacity = capacity;

    pNewDstrb->minMQ = minMQ
    pNewDstrb->hasRG = TRUE;

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

        free(pDstrb->pFragLenRange);
        kh_destroy(readGrpName, pDstrb->pReadGrpHash);

        free(pDstrb);
    }
}

SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader)
{
    const char* readGrpPos = pBamHeader->pOrigHeader->text;
    while ((readGrpPos = strstr(readGrpPos, "@RG")) != NULL)
    {
        // expand the array if necessary
        if (pDstrb->size == pDstrb->capacity)
        {
            pDstrb->capacity *= 2;
            pDstrb->pReadGrpNames = (char**) realloc(pDstrb->pReadGrpNames, pDstrb->capacity * sizeof(char*));
            if (pDstrb->pReadGrpNames == NULL)
                SR_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");
        }

        // get the name of the current read group
        const char* readGrpNamePos = strstr(readGrpPos, "ID:");
        if (readGrpNamePos != NULL)
        {
            readGrpNamePos += 3;
        }
        else
        {
            SR_ErrMsg("ERROR: the \"ID\" field is required under the read group tag.\n");
            return SR_ERR;
        }

        const char* readGrpNameEnd = strpbrk(readGrpPos, " \t\n\0");
        size_t readGrpNameLen = readGrpNameEnd - readGrpNamePos;
        if (readGrpNameLen == 0)
        {
            SR_ErrMsg("ERROR: the length of read group ID is zero.\n");
            return SR_ERR;
        }

        pDstrb->pReadGrpNames[pDstrb->size] = (char*) calloc(readGrpNameLen + 1, sizeof(char));
        if (pDstrb->pReadGrpNames[pDstrb->size] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");

        strncpy(pDstrb->pReadGrpNames[pDstrb->size], readGrpPos, readGrpNameLen);

        ++(pDstrb->size);
        ++readGrpPos;
    }

    if (pDstrb->size != 0) // read group names are found. we insert them into a hash table
    {
        pDstrb->hasRG = TRUE;
        pDstrb->pReadGrpHash = (khash_t(readGrpName)*) kh_init(readGrpName);
        kh_resize(readGrpName, pDstrb->pReadGrpHash, pDstrb->capacity);

        int khRet = 0;
        khiter_t khIter;
        for (unsigned int i = 0; i != pDstrb->size; ++i)
        {
            khIter = kh_put(readGrpName, pDstrb->pReadGrpHash, pDstrb->pReadGrpNames[i], &khRet);
            if (khRet == 0)
            {
                SR_ErrMsg("ERROR: Found a non-unique read group ID.\n");
                return SR_ERR;
            }

            kh_value((khash_t(readGrpName)*) pDstrb->pReadGrpHash, khIter) = i;
        }
    }
    else // did not find any read group names. treat all the alignments as one read group
    {
        SR_ErrMsg("WARNING: No read group is found. The alignments in this bam file will be treated as from one read group.\n");

        pDstrb->size = 1;
        pDstrb->hasRG = FALSE;
    }

    // initialize the fragment length range array according to the number of read groups
    pDstrb->pFragLenRange = (double (*)[2]) calloc(2 * pDstrb->size, sizeof(double));
    if (pDstrb->pFragLenRange == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of fragment lenght range array in the fragment lenght distribution object.\n");

    return SR_OK;
}

void SR_FragLenDstrbInitHists(SR_FragLenDstrb* pDstrb)
{
    pDstrb->pFragLenHists = (gsl_histogram**) malloc(pDstrb->size * sizeof(gsl_histogram*));
    if (pDstrb->pFragLenHists == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storageof the histogram in the fragment length distribution object.\n");

    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        double min = pDstrb->pFragLenRange[i][0];
        double max = pDstrb->pFragLenRange[i][1];
        pDstrb->pFragLenHists[i] = gsl_histogram_alloc((size_t) max - min);
        if (pDstrb->pFragLenHists[i] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storageof the histogram in the fragment length distribution object.\n");

        gsl_histogram_set_ranges_uniform(pDstrb->pFragLenHists[i], min, max);
    }
}

