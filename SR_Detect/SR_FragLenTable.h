/*
 * =====================================================================================
 *
 *       Filename:  SR_FragLenTable.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/22/2011 04:33:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef  SR_FRAGLENTABLE_H
#define  SR_FRAGLENTABLE_H

#include <stdint.h>
#include <stdio.h>

#include "SR_Types.h"

typedef struct SR_FragLenTable
{
    char** pReadGrpNames;

    void* pReadGrpHash;

    uint8_t* pairMode;

    uint64_t* totalFreq;

    double* mean;

    double* median;

    double* stdev;

    uint32_t* histIndex;

    uint32_t* fragLen;

    double*   cdf;

    uint32_t size;

    SR_Bool hasRG;

}SR_FragLenTable;


SR_FragLenTable* SR_FragLenTableAlloc(void);

void SR_FragLenTableFree(SR_FragLenTable* pFragLenTable);

void SR_FragLenTableRead(SR_FragLenTable* pFragLenTable, FILE* fragLenInput);



#endif  /*SR_FRAGLENTABLE_H*/
