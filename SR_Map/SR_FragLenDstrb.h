/*
 * =====================================================================================
 *
 *       Filename:  SR_FragLenDstrb.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/07/2011 03:06:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef  SR_FRAGLENDSTRB_H
#define  SR_FRAGLENDSTRB_H

#include <gsl/gsl_histogram.h>

#include "bam.h"
#include "SR_Types.h"
#include "SR_BamHeader.h"

typedef struct SR_FragLenDstrb
{
    char** pReadGrpNames;

    void* pReadGrpHash;

    gsl_histogram** pFragLenHists;

    uint8_t drctField;

    uint32_t size;

    uint32_t capacity;

}SR_FragLenDstrb;


inline SR_SingleOrnt SR_BamGetOrnt(bam1_t* pAlgn)
{
    if ((pAlgn->core.flag & BAM_FREAD1) != 0)
    {
        if ((pAlgn->core.flag & BAM_FREVERSE) == 0)
            return SR_1F;
        else
            return SR_1R;
    }
    else
    {
        if ((pAlgn->core.flag & BAM_FREVERSE) == 0)
            return SR_2F;
        else
            return SR_2R;
    }
}

SR_FragLenDstrb* SR_FragLenDstrbAlloc(uint8_t drctField, uint32_t capacity);

void SR_FragLenDstrbFree(SR_FragLenDstrb* pDstrb);

SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader);

void SR_FragLenDstrbSetHist(SR_FragLenDstrb* pDstrb, size_t numBins);

#endif  /*SR_FRAGLENDSTRB_H*/
