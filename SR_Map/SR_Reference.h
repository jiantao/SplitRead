/*
 * =====================================================================================
 *
 *       Filename:  SR_Reference.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/03/2011 08:05:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_REFERENCE_H
#define  SR_REFERENCE_H

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>

#include "SR_Types.h"


// reset the reference object for next reading
#define SR_ReferenceReset(reference, nextChr) \
    do                                        \
    {                                         \
        (reference)->chr = (nextChr);         \
        (reference)->length = 0;              \
    }while(0) 


// an object hold the reference sequence of a chromosome
typedef struct SR_Reference
{
    char* name;                   // chromsome name

    char* sequence;               // reference sequence

    char  md5[MD5_STR_LEN + 1];   // md5 checksum string

    int32_t  id;                  // id of the chromosome

    uint32_t nameCap;             // capacity of reference name

    uint32_t seqLen;              // length of the chromosome

    uint32_t seqCap;              // capacity of reference sequence

    off_t    filePos;             // file position of the chromosome

}SR_Reference;


// create a new reference object
SR_Reference* SR_ReferenceAlloc(void);

// free an existing reference object
void SR_ReferenceFree(SR_Reference* reference);

// read the reference sequence in the fasta file line by line, one chromosome at each time
SR_Bool SR_ReferenceLoad(SR_Reference* reference, unsigned char* nextChr, FILE* faInput);

// skip the reference sequence with unknown chromosome ID
SR_Bool SR_ReferenceSkip(unsigned char* nextChr, FILE* faInput);

// write the reference sequence into a output file in the binary format
off_t SR_ReferenceWrite(FILE* refOutput, const SR_Reference* reference);

// read the reference sequence from an input file in the binary format
SR_Bool SR_ReferenceRead(SR_Reference* reference, FILE* refInput);


#endif  /*SR_REFERENCE_H*/
