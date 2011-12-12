/*
 * =====================================================================================
 *
 *       Filename:  SR_Build_Main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/16/2011 16:16:06
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "SR_Error.h"
#include "SR_Types.h"
#include "SR_Build_GetOpt.h"
#include "SR_Reference.h"
#include "SR_OutHashTable.h"


int main(int argc, char *argv[])
{
    // load and check the parameters from the command line arguments
    SR_Build_Pars buildPars;
    SR_Build_SetPars(&buildPars, argc, argv);

    // write the hash size to the beginning of hash position index file and hash position file
    SR_ReferenceLeaveStart(buildPars.refOutput);
    SR_OutHashTableWriteStart(buildPars.hashSize, buildPars.hashTableOutput);

    // create the reference object and the reference hash table object
    SR_Reference* reference = SR_ReferenceAlloc();
    SR_RefHeader* refHeader = SR_RefHeaderAlloc(DEFAULT_NUM_CHR, DEFAULT_NUM_CHR);
    SR_OutHashTable* refHashTable = SR_OutHashTableAlloc(buildPars.hashSize);

    // a indicator of the end of the input reference file
    SR_Status status = SR_EOF;

    // read the reference sequence from the fasta file chromosome by chromosome and store it in the reference object
    // index the referen sequence with the user-specified hash size and store the hash positions in the hash position file
    // for each different hash its starting position in the hash position array will be stored in the hash position index file

    do
    {
        // this function will read the fasta file line by line until it hits a header line start with '>' or the end of file
        // when it hits the '>' character at the beginning of a line it will set the nextChr variable and return TRUE
        // when it hist the eof it will return FALSE

        if (status == SR_OK && refHeader->names[refHeader->numRefs] == NULL)
        {
            SR_ReferenceReset(reference);
            status = SR_ReferenceSkip(refHeader, buildPars.faInput);
            continue;
        }

        status = SR_ReferenceLoad(reference, refHeader, buildPars.faInput);


        // we won't get any sequence in the first round or any chromosome with an unknown ID
        // so here we skip the following steps
        if (reference->seqLen == 0)
            continue;

        // index every possible hash position in the current chromosome
        // and write the results into hash position index file and hash position file
        SR_OutHashTableLoad(refHashTable, reference->sequence, reference->seqLen, reference->id);
        int64_t htFileOffset = SR_OutHashTableWrite(refHashTable, buildPars.hashTableOutput);
        // reset the reference hash table object for next loading
        SR_OutHashTableReset(refHashTable);

        // write the reference sequence of current chromosome into the reference output file
        int64_t refFileOffset = SR_ReferenceWrite(reference, buildPars.refOutput);
        // reset the reference object for next reading
        SR_ReferenceReset(reference);

        refHeader->refFilePos[refHeader->numRefs - 1] = refFileOffset;
        refHeader->htFilePos[refHeader->numRefs - 1] = htFileOffset;

    }while(status == SR_OK);

    // handle the special references
    if (buildPars.specialRefInput != NULL)
    {
        refHeader->pSpecialRefInfo = SR_SpecialRefInfoAlloc(DEFAULT_NUM_SPECIAL_REF);

        // load the special references
        SR_SpecialRefLoad(reference, refHeader, buildPars.specialRefInput);

        // index every possible hash position in the current chromosome
        // and write the results into hash position index file and hash position file
        SR_OutHashTableLoad(refHashTable, reference->sequence, reference->seqLen, reference->id);
        int64_t htFileOffset = SR_OutHashTableWrite(refHashTable, buildPars.hashTableOutput);

        // write the reference sequence of current chromosome into the reference output file
        int64_t refFileOffset = SR_ReferenceWrite(reference, buildPars.refOutput);

        refHeader->refFilePos[refHeader->numSeqs - 1] = refFileOffset;
        refHeader->htFilePos[refHeader->numSeqs - 1] = htFileOffset;
    }
    
    int64_t refHeaderPos = SR_RefHeaderWrite(refHeader, buildPars.refOutput);
    SR_ReferenceSetStart(refHeaderPos, buildPars.refOutput);
    SR_OutHashTableSetStart(refHeaderPos, buildPars.hashTableOutput);

    // close the open files and free the allocated memory for objects
    SR_Build_Clean(reference, refHeader, refHashTable, &buildPars);

    return EXIT_SUCCESS;
}
