/*
 * =====================================================================================
 *
 *       Filename:  SR_LibInfo.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/04/2012 08:37:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <ctype.h>

#include "khash.h"
#include "SR_Error.h"
#include "SR_LibInfo.h"
#include "SR_Utilities.h"

// index of the mode count for invalid pair mode 
#define INVALID_PAIR_MODE_SET_INDEX 1

#define NUM_ZA_FIELD 7

// sample name hash
KHASH_MAP_INIT_STR(name, uint32_t);

int SR_GetNumMismatchFromZA(const char* cigarStr, unsigned int cigarLen, const char* mdStr, unsigned int mdLen)
{
    if (cigarStr == NULL)
        return -1;

    int numMM = 0;
    for (const char* cigarEnd = strpbrk(cigarStr, "DI"); cigarEnd != NULL && cigarEnd - cigarStr < cigarLen; cigarEnd = strpbrk(cigarEnd + 1, "DI"))
    {
        const char* currPos = cigarEnd - 1;
        while (isdigit(*currPos))
            --currPos;

        numMM += atoi(currPos + 1);
    }

    const char* mdFieldPos = mdStr;
    while (mdFieldPos != NULL && mdFieldPos - mdStr < mdLen)
    {
        if (isdigit(*mdFieldPos))
        {
            ++mdFieldPos;
            continue;
        }

        const char* mdFieldEnd = mdFieldPos + 1;
        while (!isdigit(*mdFieldEnd) && *mdFieldEnd != '\0')
            ++mdFieldEnd;

        if (*mdFieldPos != '^')
            numMM += mdFieldEnd - mdFieldPos;

        mdFieldPos = mdFieldEnd;
    }

    return numMM;
}

int SR_GetNumMismatchFromBam(const bam1_t* pAlgn)
{
    int numMM = 0;
    uint32_t* cigar = bam1_cigar(pAlgn);
    for (unsigned i = 0; i != pAlgn->core.n_cigar; ++i)
    {
        int type = (cigar[i] & BAM_CIGAR_MASK);
        if (type == BAM_CINS || type == BAM_CDEL)
            numMM += (cigar[i] >> BAM_CIGAR_SHIFT);
    }

    uint8_t* mdPos = bam_aux_get(pAlgn, "MD");
    if (mdPos != NULL)
    {
        const char* mdStr = bam_aux2Z(mdPos);
        const char* mdFieldPos = mdStr;
        while (mdFieldPos != NULL && *mdFieldPos != '\0')
        {
            if (isdigit(*mdFieldPos))
            {
                ++mdFieldPos;
                continue;
            }

            const char* mdFieldEnd = mdFieldPos + 1;
            while (!isdigit(*mdFieldEnd) && *mdFieldEnd != '\0')
                ++mdFieldEnd;

            if (*mdFieldPos != '^')
                numMM += mdFieldEnd - mdFieldPos;

            mdFieldPos = mdFieldEnd;
        }
    }

    return numMM;
}

static SR_PairMode SR_GetPairMode(const bam1_t* pAlignment)
{
    int8_t upMode = 0;
    int8_t downMode = 0;

    if ((pAlignment->core.flag & BAM_FREVERSE) != 0)
    {
        upMode |= 1;
    }

    if ((pAlignment->core.flag & BAM_FMREVERSE) != 0)
    {
        downMode |= 1;
    }

    if ((pAlignment->core.flag & BAM_FREAD1) != 0 && (pAlignment->core.flag & BAM_FREAD2) != 0)
        return SR_BAD_PAIR_MODE;
    else if ((pAlignment->core.flag & BAM_FREAD2) != 0)
        upMode |= (1 << 1);
    else if ((pAlignment->core.flag & BAM_FREAD1) != 0)
        downMode |= (1 << 1);
    else
        return SR_BAD_PAIR_MODE;

    if (pAlignment->core.pos > pAlignment->core.mpos)
        SR_SWAP(upMode, downMode, unsigned int);

    return (SR_PairMode) SR_PairModeMap[((upMode << 2) | downMode)];
}

/*
static SR_Status SR_LibInfoTableAddAnchor(SR_LibInfoTable* pTable, const char* tagPos)
{
    const SR_Bool isLoaded = pTable->pAnchorInfo->size == 0 ? FALSE : TRUE;
    SR_AnchorInfo* pAnchorInfo = pTable->pAnchorInfo;

    char buff[255];
    while ((tagPos = strstr(tagPos, "@SQ")) != NULL)
    {
        const char* lineEnd = strpbrk(tagPos, "\n");
        const char* refNamePos = strstr(tagPos, "SN:");

        if (refNamePos == NULL || refNamePos > lineEnd)
        {
            SR_ErrMsg("ERROR: Reference name in the bam header is not specified\n");
            return SR_ERR;
        }

        refNamePos += 3;
        const char* refNameEnd = strpbrk(refNamePos, " \t\n");
        if (refNameEnd == NULL || refNameEnd > lineEnd)
        {
            SR_ErrMsg("ERROR: Reference name in the bam header is not specified\n");
            return SR_ERR;
        }

        int refNameLen = refNameEnd - refNamePos;
        if (refNameLen <= 0)
        {
            SR_ErrMsg("ERROR: Reference name in the bam header is not specified\n");
            return SR_ERR;
        }

        int ret = 0;
        khiter_t khIter = 0;
        unsigned int refIndex = 0;
        SR_Bool isUsedRef = TRUE;

        if (isLoaded)
        {
            buff[refNameLen] = '\0';
            strncpy(buff, refNamePos, refNameLen);

            khIter = kh_put(name, pAnchorInfo->pAnchorHash, buff, &ret);
            if (ret != 0)
            {
                SR_ErrMsg("ERROR: Found a reference name that is not recorded in the previous bam files.\n");
                return SR_ERR;
            }

            refIndex = kh_value((khash_t(name)*) pAnchorInfo->pAnchorHash, khIter);
            if (pAnchorInfo->pLength[refIndex] < 0)
                isUsedRef = FALSE;
        }
        else
        {
            if (pAnchorInfo->size == pAnchorInfo->capacity)
            {
                pAnchorInfo->capacity *= 2;
                pAnchorInfo->pAnchors = (char**) realloc(pAnchorInfo->pAnchors, sizeof(char*) * pAnchorInfo->capacity);
                if (pAnchorInfo->pAnchors == NULL)
                    SR_ErrQuit("ERROR: Not enough memory for the storage of the anchor names.\n");

                pAnchorInfo->pLength = (int32_t*) realloc(pAnchorInfo->pLength, sizeof(int32_t) * pAnchorInfo->capacity);
                if (pAnchorInfo->pLength == NULL)
                    SR_ErrQuit("ERROR: Not enough memory for the storage of the anchor length.\n");

                pAnchorInfo->pMd5s = (char*) realloc(pAnchorInfo->pMd5s, sizeof(char) * MD5_STR_LEN * pAnchorInfo->capacity);
                if (pAnchorInfo->pMd5s == NULL)
                    SR_ErrQuit("ERROR: Not enough memory for the storage of the md5 strings.\n");
            }

            pAnchorInfo->pAnchors[pAnchorInfo->size] = (char*) malloc((refNameLen + 1) * sizeof(char));
            if (pAnchorInfo->pAnchors[pAnchorInfo->size] == NULL)
                SR_ErrQuit("ERROR: Not enough memory for the storage of reference name.\n");

            pAnchorInfo->pAnchors[pAnchorInfo->size][refNameLen] = '\0';
            strncpy(pAnchorInfo->pAnchors[pAnchorInfo->size], refNamePos, refNameLen);

            khIter = kh_put(name, pAnchorInfo->pAnchorHash, pAnchorInfo->pAnchors[pAnchorInfo->size], &ret);
            if (ret == 0)
            {
                SR_ErrMsg("ERROR: Found a duplicated reference name in a bam file.\n");
                return SR_ERR;
            }

            refIndex = pAnchorInfo->size;
            kh_value((khash_t(name)*) pAnchorInfo->pAnchorHash, khIter) = refIndex;

            if (strncmp("GL0", pAnchorInfo->pAnchors[refIndex], 3) == 0
                || strncmp("NC_", pAnchorInfo->pAnchors[refIndex], 3) == 0
                || strncmp("NT_", pAnchorInfo->pAnchors[refIndex], 3) == 0)
            {
                isUsedRef = FALSE;
                pAnchorInfo->pLength[refIndex] = -1;
            }
        }

        if (isUsedRef)
        {
            const char* refLenPos = strstr(tagPos, "LN:");
            if (refLenPos == NULL || refLenPos > lineEnd)
            {
                SR_ErrMsg("ERROR: Cannot find the length of the reference.\n");
                return SR_ERR;
            }

            refLenPos += 3;
            const char* refLenEnd = strpbrk(refLenPos, " \t\n");
            int refLenLen = refLenEnd - refLenPos;
            if (refLenEnd == NULL || refLenLen <= 0 || refLenEnd > lineEnd)
            {
                SR_ErrMsg("ERROR: Cannot find the length of the reference.\n");
                return SR_ERR;
            }

            buff[refLenLen] = '\0';
            strncpy(buff, refLenPos, refLenLen);

            int refLen = atoi(buff);
            if (refLen == 0)
            {
                SR_ErrMsg("ERROR: Cannot find the length of the reference.\n");
                return SR_ERR;
            }

            if (isLoaded)
            {
                if (pAnchorInfo->pLength[refIndex] > 0 && refLen != pAnchorInfo->pLength[refIndex])
                {
                    SR_ErrMsg("ERROR: The length of the reference found in this bam file is inconsistent with that in the previous bam file.\n");
                    return SR_ERR;
                }
            }
            else
                pAnchorInfo->pLength[refIndex] = refLen;

            const char* md5Pos = strstr(tagPos, "M5:");
            if (md5Pos == NULL || md5Pos > lineEnd)
            {
                SR_ErrMsg("ERROR: Cannot find the md5 string of the reference.\n");
                return SR_ERR;
            }

            md5Pos += 3;
            if (isLoaded)
            {
                strncpy(buff, md5Pos, MD5_STR_LEN);
                if (strncmp(buff, pAnchorInfo->pMd5s + refIndex * MD5_STR_LEN, MD5_STR_LEN) != 0)
                {
                    SR_ErrMsg("ERROR: The md5 string of the reference found in this bam file is inconsistent with that in the previous bam file.\n");
                    return SR_ERR;
                }
            }
            else
                strncpy(pAnchorInfo->pMd5s + refIndex * MD5_STR_LEN, md5Pos, MD5_STR_LEN);

        }
        else if (!isLoaded)
            memset(pAnchorInfo->pMd5s + refIndex * MD5_STR_LEN, '-', MD5_STR_LEN);

        if (!isLoaded)
            ++(pAnchorInfo->size);

        ++tagPos;
    }

    return SR_OK;
}

*/

static SR_Status SR_LibInfoTableAddAnchor(SR_LibInfoTable* pTable, const SR_BamHeader* pBamHeader)
{
    const SR_Bool isLoaded = pTable->pAnchorInfo->size == 0 ? FALSE : TRUE;

    int ret = 0;
    khiter_t khIter = 0;

    if (isLoaded)
    {
        if (pBamHeader->pOrigHeader->n_targets != pTable->pAnchorInfo->size)
        {
            SR_ErrMsg("ERROR: The number of reference sequences in this bam file is inconsistent with that in previous bam files.\n");
            return SR_ERR;
        }

        unsigned int refIndex = 0;
        for (unsigned int i = 0; i != pBamHeader->pOrigHeader->n_targets; ++i)
        {
            khIter = kh_put(name, pTable->pAnchorInfo->pAnchorHash, pBamHeader->pOrigHeader->target_name[i], &ret);
            if (ret == 0)
            {
                khash_t(name)* pAnchorHash = pTable->pAnchorInfo->pAnchorHash;
                refIndex = kh_value(pAnchorHash, khIter);
                if (refIndex != i)
                {
                    SR_ErrMsg("ERROR: Reference ID in this bam file is inconsistent with that in the previous bam files.\n");
                    return SR_ERR;
                }

                if (pTable->pAnchorInfo->pLength[i] > 0 && pTable->pAnchorInfo->pLength[i] != pBamHeader->pOrigHeader->target_len[i])
                {
                    SR_ErrMsg("ERROR: The length of the reference sequence in this bam file is inconsistent with that in previous bam files.\n");
                    return SR_ERR;
                }
            }
            else
            {
                SR_ErrMsg("ERROR: Found a reference sequence that is not in the previous bam headers.\n");
                return SR_ERR;
            }

            if (strncmp(pTable->pAnchorInfo->pMd5s + MD5_STR_LEN * i, pBamHeader->pMD5s[i], MD5_STR_LEN) != 0)
            {
                SR_ErrMsg("ERROR: MD5 string in this bam file is inconsistent with that in previous bam files.\n");
                return SR_ERR;
            }
        }
    }
    else
    {
        if (pBamHeader->pOrigHeader->n_targets > pTable->pAnchorInfo->capacity)
        {
            SR_AnchorInfoFree(pTable->pAnchorInfo);
            pTable->pAnchorInfo = SR_AnchorInfoAlloc(pBamHeader->pOrigHeader->n_targets);
        }

        pTable->pAnchorInfo->size = pBamHeader->pOrigHeader->n_targets;

        for (unsigned int i = 0; i != pBamHeader->pOrigHeader->n_targets; ++i)
        {
            unsigned int nameLen = strlen(pBamHeader->pOrigHeader->target_name[i]);

            pTable->pAnchorInfo->pAnchors[i] = (char*) calloc(nameLen + 1, sizeof(char));
            if (pTable->pAnchorInfo->pAnchors[i] == NULL)
                SR_ErrQuit("ERROR: Not enough memory for the storage of reference name.\n");

            strncpy(pTable->pAnchorInfo->pAnchors[i], pBamHeader->pOrigHeader->target_name[i], nameLen);
            khIter = kh_put(name, pTable->pAnchorInfo->pAnchorHash, pTable->pAnchorInfo->pAnchors[i], &ret);

            khash_t(name)* pAnchorHash = pTable->pAnchorInfo->pAnchorHash;
            kh_value(pAnchorHash, khIter) = i;

            // set the length of those unused references to -1 so that
            // we will ignore any reads aligned to them
            pTable->pAnchorInfo->pLength[i] = pBamHeader->pOrigHeader->target_len[i];
            if (strncmp("GL0", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                || strncmp("NC_", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                || strncmp("NT_", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                || strncmp("hs", pTable->pAnchorInfo->pAnchors[i], 2) == 0)
            {
                pTable->pAnchorInfo->pLength[i] = -1;
            }

            strncpy(pTable->pAnchorInfo->pMd5s + MD5_STR_LEN * i, pBamHeader->pMD5s[i], MD5_STR_LEN);
        }
    }

    return SR_OK;
}

static SR_Status SR_LibInfoTableAddSample(int* pSampleID, SR_LibInfoTable* pTable, const char* tagPos, const char* lineEnd)
{
    const char* sampleNamePos = strstr(tagPos, "SM:");
    if (sampleNamePos != NULL && sampleNamePos < lineEnd)
        sampleNamePos += 3;

    const char* sampleNameEnd = strpbrk(sampleNamePos, " \t\n");
    int sampleNameLen = sampleNameEnd - sampleNamePos;

    SR_SampleInfo* pSampleInfo = pTable->pSampleInfo;

    if (sampleNameLen > 0)
    {
        char* buff = (char*) malloc((sampleNameLen + 1) * sizeof(char));
        if (buff == NULL)
            SR_ErrQuit("ERROR: Not enought memory for the read group ID in the fragment length distribution object.\n");

        buff[sampleNameLen] = '\0';
        memcpy(buff, sampleNamePos, sampleNameLen);

        int ret = 0;

        khash_t(name)* pSampleHash = pSampleInfo->pSampleHash;
        khiter_t khIter = kh_put(name, pSampleHash, buff, &ret);

        if (ret == 0)
        {
            free(buff);
            *pSampleID = kh_value(pSampleHash, khIter);
        }
        else
        {
            if (pSampleInfo->size == pSampleInfo->capacity)
            {
                pSampleInfo->capacity *= 2;
                pSampleInfo->pSamples = (char**) realloc(pSampleInfo->pSamples, pSampleInfo->capacity * sizeof(char*));
                if (pSampleInfo->pSamples == NULL)
                    SR_ErrQuit("ERROR: Not enought memory for the sample names in the fragment length distribution object.\n");
            }

            *pSampleID = pSampleInfo->size;
            kh_value(pSampleHash, khIter) = pSampleInfo->size;
            pSampleInfo->pSamples[pSampleInfo->size] = buff;
            ++(pSampleInfo->size);
        }

        return SR_OK;
    }

    return SR_ERR;
}

static SR_Status SR_LibInfoTableAddSeqTech(SR_LibInfoTable* pTable, const char* platformPos)
{
    char buff[30];

    const char* platformEnd = strpbrk(platformPos, " \t\n");
    unsigned int platformLen = platformEnd - platformPos;
    if (platformLen > 29)
        return SR_ERR;

    strncpy(buff, platformPos, platformLen);
    buff[platformLen] = '\0';

    StrToUpper(buff);
    unsigned int lastIndex = pTable->size;
    if (strstr(buff, "LONG") != NULL)
        pTable->pSeqTech[lastIndex] = ST_ILLUMINA_LONG;
    else if (strcmp(buff, "ILLUMINA") == 0)
        pTable->pSeqTech[lastIndex] = ST_ILLUMINA;
    else if (strcmp(buff, "LS454") == 0)
        pTable->pSeqTech[lastIndex] = ST_454;
    else if (strcmp(buff, "SOLID") == 0)
        pTable->pSeqTech[lastIndex] = ST_SOLID;
    else
        return SR_ERR;

    return SR_OK;
}

static SR_Status SR_LibInfoTableAddReadGrp(SR_LibInfoTable* pTable, const char* tagPos, const char* lineEnd, int sampleID)
{
    // get the name of the current read group
    const char* readGrpNamePos = strstr(tagPos, "ID:");
    if (readGrpNamePos != NULL && readGrpNamePos < lineEnd)
        readGrpNamePos += 3;
    else
        return SR_ERR;

    const char* platformPos = strstr(tagPos, "PL:");
    if (platformPos != NULL && platformPos < lineEnd)
        platformPos += 3;
    else
        return SR_ERR;

    // expand the array if necessary
    if (pTable->size == pTable->capacity)
    {
        pTable->capacity *= 2;
        pTable->pReadGrps = (char**) realloc(pTable->pReadGrps, pTable->capacity * sizeof(char*));
        if (pTable->pReadGrps == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of read group names in the library information object.\n");

        pTable->pSampleMap = (int32_t*) realloc(pTable->pSampleMap, pTable->capacity * sizeof(int32_t));
        if (pTable->pSampleMap == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the sample ID map in the library information object.\n");

        pTable->pSeqTech = (int8_t*) realloc(pTable->pSeqTech, pTable->capacity * sizeof(int8_t));
        if (pTable->pSeqTech == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the sequencing technology in the library information  object.\n");
    }

    SR_Status status = SR_LibInfoTableAddSeqTech(pTable, platformPos);
    if (status != SR_OK)
        return SR_ERR;

    const char* readGrpNameEnd = strpbrk(readGrpNamePos, " \t\n\0");
    size_t readGrpNameLen = readGrpNameEnd - readGrpNamePos;
    if (readGrpNameLen > 0)
    {
        pTable->pReadGrps[pTable->size] = (char*) calloc(readGrpNameLen + 1, sizeof(char));
        if (pTable->pReadGrps[pTable->size] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");

        memcpy(pTable->pReadGrps[pTable->size], readGrpNamePos, readGrpNameLen);

        int ret = 0;

        khash_t(name)* pReadGrpHash = pTable->pReadGrpHash;
        khiter_t khIter = kh_put(name, pReadGrpHash, pTable->pReadGrps[pTable->size], &ret);

        if (ret != 0)
        {
            pTable->pSampleMap[pTable->size] = sampleID;
            kh_value(pReadGrpHash, khIter) = pTable->size;
            ++(pTable->size);

            return SR_OK;
        }
        else
        {
            free(pTable->pReadGrps[pTable->size]);
            pTable->pReadGrps[pTable->size] = NULL;
            SR_ErrMsg("ERROR: Found a duplicated read group ID.\n");
        }
    }

    return SR_ERR;
}

static void SR_AnchorInfoWrite(const SR_AnchorInfo* pInfo, FILE* libFile)
{
    unsigned int writeSize = 0;

    writeSize = fwrite(pInfo->pLength, sizeof(int32_t), pInfo->size, libFile);
    if (writeSize != pInfo->size)
        SR_ErrQuit("ERROR: Cannot write the length of anchors into the information file.\n");

    writeSize = fwrite(pInfo->pMd5s, sizeof(char), pInfo->size * MD5_STR_LEN, libFile);
    if (writeSize != pInfo->size * MD5_STR_LEN)
        SR_ErrQuit("ERROR: Cannot write the md5 string into the information file.\n");

    for (unsigned int i = 0; i != pInfo->size; ++i)
    {
        uint32_t anchorNameLen = strlen(pInfo->pAnchors[i]);
        writeSize = fwrite(&(anchorNameLen), sizeof(uint32_t), 1, libFile);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Cannot write the anchor name length into the information file.\n");

        writeSize = fwrite(pInfo->pAnchors[i], sizeof(char), anchorNameLen, libFile);
        if (writeSize != anchorNameLen)
            SR_ErrQuit("ERROR: Cannot write the sample name into the information file.\n");
    }
}

static void SR_SampleInfoWrite(const SR_SampleInfo* pInfo, FILE* libFile)
{
    unsigned int writeSize = 0;

    for (unsigned int i = 0; i != pInfo->size; ++i)
    {
        uint32_t sampleNameLen = strlen(pInfo->pSamples[i]);
        writeSize = fwrite(&(sampleNameLen), sizeof(uint32_t), 1, libFile);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Cannot write the sample name length into the information file.\n");

        writeSize = fwrite(pInfo->pSamples[i], sizeof(char), sampleNameLen, libFile);
        if (writeSize != sampleNameLen)
            SR_ErrQuit("ERROR: Cannot write the sample name into the information file.\n");
    }
}

static void SR_AnchorInfoRead(SR_AnchorInfo* pInfo, FILE* libFile)
{
    unsigned int readSize = 0;
    
    readSize = fread(pInfo->pLength, sizeof(int32_t), pInfo->size, libFile);
    if (readSize != pInfo->size)
        SR_ErrQuit("ERROR: Cannot read the length of anchors from the library file.\n");

    readSize = fread(pInfo->pMd5s, sizeof(char), pInfo->size * MD5_STR_LEN, libFile);
    if (readSize != pInfo->size * MD5_STR_LEN)
        SR_ErrQuit("ERROR: Cannot read the md5 strings from the library file.\n");

    uint32_t anchorNameLen = 0;
    for (unsigned int i = 0; i != pInfo->size; ++i)
    {
        readSize = fread(&anchorNameLen, sizeof(uint32_t), 1, libFile);
        if (readSize != 1)
            SR_ErrQuit("ERROR: Cannot read the length of anchor name from the library file.\n");

        pInfo->pAnchors[i] = (char*) malloc((anchorNameLen + 1) * sizeof(char));
        if (pInfo->pAnchors[i] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the anchor name.\n");

        pInfo->pAnchors[i][anchorNameLen] = '\0';
        
        readSize = fread(pInfo->pAnchors[i], sizeof(char), anchorNameLen, libFile);
        if (readSize != anchorNameLen)
            SR_ErrQuit("ERROR: Cannot read the anchor name from the library file.\n");

        int ret = 0;
        khiter_t khIter = kh_put(name, pInfo->pAnchorHash, pInfo->pAnchors[i], &ret);
        kh_value((khash_t(name)*) pInfo->pAnchorHash, khIter) = i;
    }
}

static void SR_SampleInfoRead(SR_SampleInfo* pInfo, FILE* libFile)
{
    unsigned int readSize = 0;
    uint32_t sampleNameLen = 0;

    for (unsigned int i = 0; i != pInfo->size; ++i)
    {
        readSize = fread(&sampleNameLen, sizeof(uint32_t), 1, libFile);
        if (readSize != 1)
            SR_ErrQuit("ERROR: Cannot read the length of sample name from the library file.\n");

        pInfo->pSamples[i] = (char*) malloc((sampleNameLen + 1) * sizeof(char));
        if (pInfo->pSamples[i] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the sample name.\n");

        pInfo->pSamples[i][sampleNameLen] = '\0';
        
        readSize = fread(pInfo->pSamples[i], sizeof(char), sampleNameLen, libFile);
        if (readSize != sampleNameLen)
            SR_ErrQuit("ERROR: Cannot read the sample name from the library file.\n");

        int ret = 0;
        khiter_t khIter = kh_put(name, pInfo->pSampleHash, pInfo->pSamples[i], &ret);
        kh_value((khash_t(name)*) pInfo->pSampleHash, khIter) = i;
    }
}


SR_AnchorInfo* SR_AnchorInfoAlloc(uint32_t capacity)
{
    SR_AnchorInfo* pNewInfo = (SR_AnchorInfo*) malloc(sizeof(SR_AnchorInfo));
    if (pNewInfo == NULL)
        SR_ErrQuit("ERROR: Not enough memory for an anchor information object.\n");

    pNewInfo->pAnchors = (char**) malloc(capacity * sizeof(char*));
    if (pNewInfo->pAnchors == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of anchor names in the anchor information object.\n");

    pNewInfo->pLength = (int32_t*) malloc(capacity * sizeof(int32_t));
    if (pNewInfo->pLength == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of anchor length in the anchor information object.\n");

    pNewInfo->pMd5s = (char*) malloc(capacity * MD5_STR_LEN * sizeof(char));
    if (pNewInfo->pMd5s == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of md5 string in the anchor information object.\n");

    pNewInfo->size = 0;
    pNewInfo->capacity = capacity;

    pNewInfo->pAnchorHash = kh_init(name);
    kh_resize(name, pNewInfo->pAnchorHash, 2 * capacity);

    return pNewInfo;
}

void SR_AnchorInfoFree(SR_AnchorInfo* pInfo)
{
    if (pInfo != NULL)
    {
        for (unsigned int i = 0; i != pInfo->size; ++i)
            free(pInfo->pAnchors[i]);

        kh_destroy(name, pInfo->pAnchorHash);

        free(pInfo->pAnchors);
        free(pInfo->pLength);
        free(pInfo->pMd5s);
        free(pInfo);
    }
}

SR_SampleInfo* SR_SampleInfoAlloc(uint32_t capacity)
{
    SR_SampleInfo* pNewInfo = (SR_SampleInfo*) malloc(sizeof(SR_SampleInfo));
    if (pNewInfo == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a sample information object.\n");

    pNewInfo->pSamples = (char**) malloc(capacity * sizeof(char*));
    if (pNewInfo->pSamples == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of sample names in the sample information object.\n");

    pNewInfo->pReadFraction = (double*) malloc(capacity * sizeof(double));
    if (pNewInfo->pReadFraction == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of read fraction in the sample information object.\n");

    pNewInfo->pSampleHash = kh_init(name);
    kh_resize(name, pNewInfo->pSampleHash, 2 * capacity);

    pNewInfo->size = 0;
    pNewInfo->capacity = capacity;

    return pNewInfo;
}

void SR_SampleInfoFree(SR_SampleInfo* pInfo)
{
    if (pInfo != NULL)
    {
        if (pInfo->pSamples != NULL)
        {
            for (unsigned int i = 0; i != pInfo->size; ++i)
                free(pInfo->pSamples[i]);

            free(pInfo->pSamples);
        }

        kh_destroy(name, pInfo->pSampleHash);
        free(pInfo->pReadFraction);
        free(pInfo);
    }
}

SR_LibInfoTable* SR_LibInfoTableAlloc(uint32_t capAnchor, uint32_t capSample, uint32_t capReadGrp)
{
    SR_LibInfoTable* pNewTable = (SR_LibInfoTable*) malloc(sizeof(SR_LibInfoTable));
    if (pNewTable == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a library information table object.\n");

    pNewTable->pSampleInfo = SR_SampleInfoAlloc(capSample);

    pNewTable->pLibInfo = (SR_LibInfo*) malloc(capReadGrp * sizeof(SR_LibInfo));
    if (pNewTable->pLibInfo == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of library information in an library table object.\n");

    pNewTable->pReadGrps = (char**) malloc(capReadGrp * sizeof(char*));
    if (pNewTable->pLibInfo == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of read group names in an library table object.\n");

    pNewTable->pSampleMap = (int32_t*) malloc(capReadGrp * sizeof(int32_t));
    if (pNewTable->pSampleMap == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of read-group-to-sample map in an library table object.\n");

    pNewTable->pSeqTech = (int8_t*) malloc(capReadGrp * sizeof(int8_t));
    if (pNewTable->pSeqTech == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of sequencing technologies in an library table object.\n");

    pNewTable->pAnchorInfo = SR_AnchorInfoAlloc(capAnchor);

    pNewTable->pReadGrpHash = kh_init(name);
    kh_resize(name, pNewTable->pReadGrpHash, 2 * capReadGrp);

    pNewTable->size = 0;
    pNewTable->capacity = capReadGrp;
    pNewTable->fragLenMax = 0;
    pNewTable->cutoff = 0.0;
    pNewTable->trimRate = 0.0;

    return pNewTable;
}

void SR_LibInfoTableFree(SR_LibInfoTable* pTable)
{
    if (pTable != NULL)
    {
        if (pTable->pReadGrps != NULL)
        {
            for (unsigned int i = 0; i != pTable->size; ++i)
                free(pTable->pReadGrps[i]);

            free(pTable->pReadGrps);
        }

        SR_SampleInfoFree(pTable->pSampleInfo);
        SR_AnchorInfoFree(pTable->pAnchorInfo);

        kh_destroy(name, pTable->pReadGrpHash);

        free(pTable->pLibInfo);
        free(pTable->pSampleMap);
        free(pTable->pSeqTech);

        free(pTable);
    }
}

SR_Status SR_LoadZAtag(SR_ZAtag* pZAtag, const bam1_t* pUpAlgn)
{
    uint8_t* zaPos = bam_aux_get(pUpAlgn, "ZA");
    const char* zaStr = NULL;
    if (zaPos != NULL)
        zaStr = bam_aux2Z(zaPos);
    else
        return SR_NOT_FOUND;

    // set the first character of the special reference name to be ' '
    // we can later check if the pair hit any special references
    pZAtag->spRef[0][0] = ' ';
    pZAtag->spRef[1][0] = ' ';

    pZAtag->spRef[0][2] = '\0';
    pZAtag->spRef[1][2] = '\0';

    unsigned int whichMate = 0;
    const char* currFieldPos = zaStr + 1;
    if (*currFieldPos == '&')
        whichMate = 1;

    currFieldPos += 2;
    unsigned int currFieldNum = 1;
    const char* currFieldEnd = NULL;

    const char* cigarStr = NULL;
    unsigned int cigarLen = 0;

    const char* mdStr = NULL;
    unsigned int mdLen = 0;

    do
    {
        ++currFieldNum;
        currFieldEnd = strpbrk(currFieldPos, ";>");

        if (currFieldEnd - currFieldPos == 0)
        {
            ++currFieldPos;
            continue;
        }

        switch (currFieldNum)
        {
            case 2:
                pZAtag->bestMQ[whichMate] = atoi(currFieldPos);
                break;
            case 3:
                pZAtag->secMQ[whichMate] = atoi(currFieldPos);
                break;
            case 4:
                pZAtag->spRef[whichMate][0] = *currFieldPos;
                pZAtag->spRef[whichMate][1] = *(currFieldPos + 1);
                break;
            case 5:
                pZAtag->numMappings[whichMate] = atoi(currFieldPos);
                break;
            case 6:
                cigarStr = currFieldPos;
                cigarLen = currFieldEnd - currFieldPos;
                break;
            case 7:
                mdStr = currFieldPos;
                mdLen = currFieldEnd - currFieldPos;
                break;
            default:
                break;
        }

        currFieldPos = currFieldEnd + 1;

    }while(*currFieldEnd != '>');

    if (cigarStr != NULL && mdStr != NULL)
        pZAtag->numMM[whichMate] = SR_GetNumMismatchFromZA(cigarStr, cigarLen, mdStr, mdLen);
    else
        pZAtag->numMM[whichMate] = SR_GetNumMismatchFromBam(pUpAlgn);


    // handle the sencond ZA tag
    currFieldPos += 3;
    currFieldNum = 1;
    currFieldEnd = NULL;

    cigarStr = NULL;
    cigarLen = 0;

    mdStr = NULL;
    mdLen = 0;

    whichMate ^= 1;

    do
    {
        ++currFieldNum;
        currFieldEnd = strpbrk(currFieldPos, ";>");

        if (currFieldEnd - currFieldPos == 0)
        {
            ++currFieldPos;
            continue;
        }

        switch (currFieldNum)
        {
            case 2:
                pZAtag->bestMQ[whichMate] = atoi(currFieldPos);
                break;
            case 3:
                pZAtag->secMQ[whichMate] = atoi(currFieldPos);
                break;
            case 4:
                pZAtag->spRef[whichMate][0] = *currFieldPos;
                pZAtag->spRef[whichMate][1] = *(currFieldPos + 1);
                break;
            case 5:
                pZAtag->numMappings[whichMate] = atoi(currFieldPos);
                break;
            case 6:
                cigarStr = currFieldPos;
                cigarLen = currFieldEnd - currFieldPos;
                break;
            case 7:
                mdStr = currFieldPos;
                mdLen = currFieldEnd - currFieldPos;
            default:
                break;
        }

        currFieldPos = currFieldEnd + 1;

    }while(*currFieldEnd != '>');

    if (cigarStr != NULL && mdStr != NULL)
        pZAtag->numMM[whichMate] = SR_GetNumMismatchFromZA(cigarStr, cigarLen, mdStr, mdLen);
    else
        pZAtag->numMM[whichMate] = SR_GetNumMismatchFromBam(pUpAlgn);

    return SR_OK;
}

SR_Status SR_LoadPairStats(SR_PairStats* pPairStats, const bam1_t* pAlignment, const SR_LibInfoTable* pTable)
{
    pPairStats->pairMode = SR_GetPairMode(pAlignment);
    if (pPairStats->pairMode == SR_BAD_PAIR_MODE)
        return SR_ERR;

    pPairStats->fragLen = abs(pAlignment->core.isize);

    if (pAlignment->core.tid != pAlignment->core.mtid)
        pPairStats->fragLen = -1;

    static const char tagRG[2] = {'R', 'G'};
    uint8_t* rgPos = bam_aux_get(pAlignment, tagRG);
    if (rgPos != NULL)
    {
        const char* RG = bam_aux2Z(rgPos);
        SR_Status status = SR_LibInfoTableGetRGIndex(&(pPairStats->readGrpID), pTable, RG);
        if (status != SR_OK)
            return SR_ERR;

        SR_SeqTech seqTech = pTable->pSeqTech[pPairStats->readGrpID];
        pPairStats->pairMode = SR_SeqTechMap[seqTech][pPairStats->pairMode];
    }
    else
        return SR_ERR;

    return SR_OK;
}

SR_Bool SR_IsNormalPair(SR_PairStats* pPairStats, unsigned int* pBackHistIndex, 
                        const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const SR_LibInfoTable* pTable, unsigned short minMQ)
{
    if (pUpAlgn->core.qual < minMQ || pDownAlgn->core.qual < minMQ)
        return FALSE;

    SR_Status status = SR_LoadPairStats(pPairStats, pUpAlgn, pTable);
    if (status != SR_OK)
        return FALSE;

    *pBackHistIndex = pTable->size - pPairStats->readGrpID;

    if (SV_ReadPairTypeMap[0][pPairStats->pairMode] != SV_NORMAL && SV_ReadPairTypeMap[1][pPairStats->pairMode] != SV_NORMAL)
        return FALSE;

    return TRUE;
}

SR_Status SR_LibInfoTableSetRG(SR_LibInfoTable* pTable, unsigned int* oldSize, const SR_BamHeader* pBamHeader)
{
    SR_Status status = SR_OK;
    *oldSize = pTable->size;
    unsigned int oldCapacity = pTable->capacity;

    const char* tagPos = pBamHeader->pOrigHeader->text;
    status = SR_LibInfoTableAddAnchor(pTable, pBamHeader);
    if (status != SR_OK)
        return SR_ERR;

    while ((tagPos = strstr(tagPos, "@RG")) != NULL)
    {
        const char* lineEnd = strpbrk(tagPos, "\n");
        int32_t sampleID = 0;

        status = SR_LibInfoTableAddSample(&sampleID, pTable, tagPos, lineEnd);
        if (status != SR_OK)
        {
            SR_ErrMsg("ERROR: the \"SM\" field is not found under the read group tag in the bam header.\n");

            ++tagPos;
            continue;
        }

        status = SR_LibInfoTableAddReadGrp(pTable, tagPos, lineEnd, sampleID);
        if (status != SR_OK)
            SR_ErrMsg("ERROR: the \"ID\" field is required under the read group tag.\n");

        ++tagPos;
    }

    if (pTable->capacity > oldCapacity)
    {
        pTable->pLibInfo = (SR_LibInfo*) realloc(pTable->pLibInfo, sizeof(SR_LibInfo) * pTable->capacity);
        if (pTable->pLibInfo == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of library summary in the library information object.\n");
    }

    return SR_OK;
}

SR_Status SR_LibInfoTableGetRGIndex(int32_t* pReadGrpIndex, const SR_LibInfoTable* pTable, const char* pReadGrpName)
{
    *pReadGrpIndex = 0;

    khash_t(name)* pRgHash = pTable->pReadGrpHash;
    khiter_t khIter = kh_get(name, pRgHash, pReadGrpName);
    if (khIter != kh_end(pRgHash))
    {
        *pReadGrpIndex = kh_value(pRgHash, khIter);
    }
    else
    {
        SR_ErrMsg("ERROR: Found a read group name that is not recorded in the library information table.\n");
        return SR_ERR;
    }

    return SR_OK;
}


void SR_LibInfoTableUpdate(SR_LibInfoTable* pTable, const SR_FragLenHistArray* pHistArray, unsigned int oldSize)
{
    for (unsigned int i = 0; i != pHistArray->size; ++i)
    {
        const SR_FragLenHist* pHist = pHistArray->data + i;

        pTable->pLibInfo[oldSize].fragLenMedian = pHist->median;

        double oneSideFlow = pTable->trimRate / 2.0 * pHist->modeCount[0];
        double total = pHist->modeCount[0] * (1.0 - pTable->trimRate);

        double cumFreq = 0.0;
        double oneSideCutoff = pTable->cutoff / 2;
        for (unsigned int i = 0; i != pHist->size; ++i)
        {
            cumFreq += pHist->freq[i];
            if (((cumFreq - oneSideFlow) / total) > oneSideCutoff)
            {
                pTable->pLibInfo[oldSize].fragLenLow = pHist->fragLen[i];
                break;
            }
        }

        cumFreq = 0.0;
        for (int i = pHist->size - 1; i != -1; --i)
        {
            cumFreq += pHist->freq[i];
            if (((cumFreq - oneSideFlow) / total) > oneSideCutoff)
            {
                pTable->pLibInfo[oldSize].fragLenHigh = pHist->fragLen[i];
                break;
            }
        }

        if (pTable->pLibInfo[oldSize].fragLenHigh > pTable->fragLenMax)
            pTable->fragLenMax = pTable->pLibInfo[oldSize].fragLenHigh;

        ++oldSize;
    }
}

void SR_LibInfoTableWrite(const SR_LibInfoTable* pTable, FILE* libFile)
{
    unsigned int writeSize = 0;

    writeSize = fwrite(&(pTable->pAnchorInfo->size), sizeof(uint32_t), 1, libFile);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the number of anchors into the information file.\n");

    writeSize = fwrite(&(pTable->pSampleInfo->size), sizeof(uint32_t), 1, libFile);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the number of samples into the information file.\n");

    writeSize = fwrite(&(pTable->size), sizeof(uint32_t), 1, libFile);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the number of read groups into the information file.\n");

    SR_AnchorInfoWrite(pTable->pAnchorInfo, libFile);
    SR_SampleInfoWrite(pTable->pSampleInfo, libFile);

    writeSize = fwrite(pTable->pSampleMap, sizeof(int32_t), pTable->size, libFile);
    if (writeSize != pTable->size)
        SR_ErrQuit("ERROR: Cannot write the read-group-to-sample map into the information file.\n");

    for (unsigned int i = 0; i != pTable->size; ++i)
    {
        uint32_t readGrpNameLen = strlen(pTable->pReadGrps[i]);
        writeSize = fwrite(&(readGrpNameLen), sizeof(uint32_t), 1, libFile);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Cannot write the readGrp name length into the information file.\n");

        writeSize = fwrite(pTable->pReadGrps[i], sizeof(char), readGrpNameLen, libFile);
        if (writeSize != readGrpNameLen)
            SR_ErrQuit("ERROR: Cannot write the readGrp name into the information file.\n");
    }

    writeSize = fwrite(&(pTable->fragLenMax), sizeof(uint32_t), 1, libFile);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the maximum fragment length into the library file.\n");

    writeSize = fwrite(&(pTable->cutoff), sizeof(double), 1, libFile);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the cutoff into the library information file.\n");

    writeSize = fwrite(&(pTable->trimRate), sizeof(double), 1, libFile);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the trim rate into the library information file.\n");

    writeSize = fwrite(pTable->pLibInfo, sizeof(SR_LibInfo), pTable->size, libFile);
    if (writeSize != pTable->size)
        SR_ErrQuit("ERROR: Cannot write the library information into the output file.\n");

    fflush(libFile);
}


SR_LibInfoTable* SR_LibInfoTableRead(FILE* libFile)
{
    unsigned int readSize = 0;
    uint32_t sizeAC = 0;
    uint32_t sizeSM = 0;
    uint32_t sizeRG = 0;

    readSize = fread(&sizeAC, sizeof(uint32_t), 1, libFile);
    if (readSize != 1)
        SR_ErrQuit("ERROR: Cannot read the number of anchors from the library file.\n");

    readSize = fread(&sizeSM, sizeof(uint32_t), 1, libFile);
    if (readSize != 1)
        SR_ErrQuit("ERROR: Cannot read the number of sample from the library file.\n");

    readSize = fread(&sizeRG, sizeof(uint32_t), 1, libFile);
    if (readSize != 1)
        SR_ErrQuit("ERROR: Cannot read the number of read group from the library file.\n");

    SR_LibInfoTable* pTable = SR_LibInfoTableAlloc(sizeAC, sizeSM, sizeRG);

    pTable->pAnchorInfo->size = sizeAC;
    pTable->pSampleInfo->size = sizeSM;
    pTable->size = sizeRG;

    pTable->pAnchorInfo->capacity = sizeAC;
    pTable->pSampleInfo->capacity = sizeSM;
    pTable->capacity = sizeRG;

    SR_AnchorInfoRead(pTable->pAnchorInfo, libFile);
    SR_SampleInfoRead(pTable->pSampleInfo, libFile);

    readSize = fread(pTable->pSampleMap, sizeof(int32_t), sizeRG, libFile);
    if (readSize != sizeRG)
        SR_ErrQuit("ERROR: Cannot read the read-group-to-sample map from the library file.\n");

    uint32_t readGrpNameLen = 0;

    for (unsigned int i = 0; i != sizeRG; ++i)
    {
        readSize = fread(&readGrpNameLen, sizeof(uint32_t), 1, libFile);
        if (readSize != 1)
            SR_ErrQuit("ERROR: Cannot read the length of read group name from the library file.\n");

        pTable->pReadGrps[i] = (char*) malloc((readGrpNameLen + 1) * sizeof(char));
        if (pTable->pReadGrps[i] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the sample name.\n");

        pTable->pReadGrps[i][readGrpNameLen] = '\0';

        readSize = fread(pTable->pReadGrps[i], sizeof(char), readGrpNameLen, libFile);
        if (readSize != readGrpNameLen)
            SR_ErrQuit("ERROR: Cannot read the read group name from the library file.\n");

        int ret = 0;
        khiter_t khIter = kh_put(name, pTable->pReadGrpHash, pTable->pReadGrps[i], &ret);
        kh_value((khash_t(name)*) pTable->pReadGrpHash, khIter) = i;
    }

    readSize = fread(&(pTable->fragLenMax), sizeof(uint32_t), 1, libFile);
    if (readSize != 1)
        SR_ErrQuit("ERROR: Cannot read the maximum fragment length from the library file.\n");

    readSize = fread(&(pTable->cutoff), sizeof(double), 1, libFile);
    if (readSize != 1)
        SR_ErrQuit("ERROR: Cannot read the cutoff from the library file.\n");

    readSize = fread(&(pTable->trimRate), sizeof(double), 1, libFile);
    if (readSize != 1)
        SR_ErrQuit("ERROR: Cannot read the trim rate from the library file.\n");

    readSize = fread(pTable->pLibInfo, sizeof(SR_LibInfo), pTable->size, libFile);
    if (readSize != pTable->size)
        SR_ErrQuit("ERROR: Cannot read the library information from the library file.\n");

    return pTable;
}
