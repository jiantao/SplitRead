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
    WriteHashSize(buildPars.hashTableOutput, buildPars.hashSize);

    // create the reference object and the reference hash table object
    SR_Reference* reference = SR_ReferenceAlloc(DEFAULT_REF_CAPACITY);
    SR_OutHashTable* refHashTable = SR_OutHashTableAlloc(buildPars.hashSize);

    // next chromosome ID in the reference input file
    unsigned char nextChr = DEFAULT_START_CHR;

    // a indicator of the end of the input reference file
    Bool keepLoading = TRUE;

    // read the reference sequence from the fasta file chromosome by chromosome and store it in the reference object
    // index the referen sequence with the user-specified hash size and store the hash positions in the hash position file
    // for each different hash its starting position in the hash position array will be stored in the hash position index file

    do
    {
        // this function will read the fasta file line by line until it hits a header line start with '>' or the end of file
        // when it hits the '>' character at the beginning of a line it will set the nextChr variable and return TRUE
        // when it hist the eof it will return FALSE

        if (nextChr != INVALID_CHR_ID)
            keepLoading = SR_ReferenceLoad(reference, &nextChr, buildPars.faInput);
        else
        {
            SR_ErrMsg("WARNING: Found unrecognized chromosome ID. This chromosome will be skipped.\n");
            keepLoading = SR_ReferenceSkip(&nextChr, buildPars.faInput);
        }

        // we won't get any sequence in the first round or any chromosome with an unknown ID
        // so here we skip the following steps
        if (reference->length == 0)
        {
            reference->chr = nextChr;
            continue;
        }

        // index every possible hash position in the current chromosome
        // and write the results into hash position index file and hash position file
        SR_OutHashTableLoad(refHashTable, reference->md5, reference->sequence, reference->length, reference->chr);
        off_t hsFileOffset = SR_OutHashTableWrite(buildPars.hashTableOutput, refHashTable);
        // reset the reference hash table object for next loading
        SR_OutHashTableReset(refHashTable);

        // write the reference sequence of current chromosome into the reference output file
        off_t refFileOffset = SR_ReferenceWrite(buildPars.refOutput, reference);
        // reset the reference object for next reading

        // if offset file of reference and hash table is specified 
        // we write the offset position of each chromosome into the offset file.
        if (buildPars.offsetOutput != NULL)
            fprintf(buildPars.offsetOutput, "%u\t%lu\t%lu\n", reference->chr, refFileOffset, hsFileOffset);

        SR_ReferenceReset(reference, nextChr);


    }while(keepLoading);

    // close the open files and free the allocated memory for objects
    SR_Build_Clean(reference, refHashTable, &buildPars);

    return EXIT_SUCCESS;
}
