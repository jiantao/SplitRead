/*
 * =====================================================================================
 *
 *       Filename:  SR_ReadPairDetector.h
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

#ifndef  SR_READPAIRDETECTOR_H
#define  SR_READPAIRDETECTOR_H

#include "SR_FragLenDstrb.h"

typedef enum
{
    SV_UNKNOWN            = 0,

    SV_DELETION           = 1,

    SV_INSERTION          = 2,

    SV_MEI_INSERTION      = 4,

    SV_TANDEM_DUP         = 8,

    SV_INVERSION          = 16,

    SV_INTRA_CHR_TRNSLCTN = 32,

    SV_INTER_CHR_TRNSLCTN = 64

}SV_EventType;

typedef struct SR_ReadPairInfo
{
    int32_t readGrpID;

    int32_t upRefID;

    union
    {
        int32_t dowRefID;
        int32_t fragLen;
    };

    int32_t upPos;

    int32_t downPos;

    int32_t probPower:8, probNum:8, eventType:8, pairMode:8;

}SR_ReadPairInfo;

typedef struct SR_ReadPairInfoArray
{
    SR_ReadPairInfo* data;

    uint32_t size;

    uint32_t capacity;

}SR_ReadPairInfoArray;

SR_Bool SR_ReadPairFilter(SR_BamNode* pBamNode, const void* filterData);

static inline SR_Bool SR_IsSVPair(SR_BamNode** ppUpAlgn, SR_BamNode** ppDownAlgn, unsigned short minMQ)
{
    if ((*ppUpAlgn)->alignment.core.qual < minMQ 
        || (*ppDownAlgn)->alignment.core.qual < minMQ)
    {
        return FALSE;
    }

    return TRUE;
}

SR_Status SR_ReadPairInfoLoad(SR_ReadPairInfo* pInfo, const SR_BamNode* pUpAlgn, const SR_BamNode* pDownAlgn, const SR_FragLenDstrb* pDstrb);

#endif  /*SR_READPAIRDETECTOR_H*/
