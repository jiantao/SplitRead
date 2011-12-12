/*
 * =====================================================================================
 *
 *       Filename:  SR_Build_GetOpt.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/29/2011 06:42:36 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <string.h>

#include "SR_Error.h"
#include "SR_Build_GetOpt.h"

// total number of arguments we should expect for the split-read build program
#define OPT_BUILD_TOTAL_NUM 6

// total number of required arguments we should expect for the split-read build program
#define OPT_BUILD_REQUIRED_NUM 4

// the index of show help in the option object array
#define OPT_HELP            0

// the index of the fasta input file in the option object array
#define OPT_FA_INPUT_FILE   1

// the index of the reference output file in the option object array
#define OPT_REF_OUTPUT_FILE 2

// the index of the output hash table file in the option object array
#define OPT_HASH_TABLE_FILE 3

// the index of the hash size in the option object array
#define OPT_HASH_SIZE       4

#define OPT_SPECIAL_REF_INPUT 5


// get the options from command line arguemnts
int SR_GetOpt(SR_Option opts[], int argc, char* argv[])
{
    int optNum = 0;
    SR_Bool hasParsed = FALSE;

    for (unsigned int i = 1; i < argc; ++i)
    {
        if (hasParsed)
        {
            hasParsed = FALSE;
            continue;
        }

        if (argv[i][0] != '-')
            SR_ErrQuit("ERROR: Invalid argument %s.\n", argv[i]);

        const char* currOpt = argv[i] + 1;

        for (unsigned int j = 0; ; ++j)
        {
            if (opts[j].name == NULL)
                SR_ErrQuit("ERROR: Ivalid argument \"%s\".\n", argv[i]);

            if (strcmp(currOpt, opts[j].name) == 0)
            {
                opts[j].isFound = TRUE;
                if (i + 1 != argc && argv[i+1][0] != '-')
                {
                    hasParsed = TRUE;
                    opts[j].value = argv[i + 1];
                }
                else
                    opts[j].value = NULL;

                ++optNum;
                break;
            }
        }
    }

    return optNum;
}


// set the parameters for the split-read build program from the parsed command line arguments 
void SR_Build_SetPars(SR_Build_Pars* pars, int argc, char* argv[])
{
    SR_Option opts[] = 
    {
        {"help", NULL, FALSE},
        {"fi",   NULL, FALSE},
        {"ro",   NULL, FALSE},
        {"hto",  NULL, FALSE},
        {"hs",   NULL, FALSE},
        {"sfi",  NULL, FALSE},
        {NULL,   NULL, FALSE}
    };

    int optNum = SR_GetOpt(opts, argc, argv);

    for (unsigned int i = 0; i != OPT_BUILD_TOTAL_NUM; ++i)
    {
        switch (i)
        {
            case OPT_HELP:
                if (opts[i].isFound)
                    SR_Build_ShowHelp();
                break;
            case OPT_FA_INPUT_FILE:
                if (opts[i].value == NULL)
                    SR_ErrQuit("ERROR: The input fasta file is not specified.\n");

                pars->faInput = fopen(opts[i].value, "r");
                if (pars->faInput == NULL)
                    SR_ErrSys("ERROR: Cannot open the fasta file \"%s\" for reading.\n", opts[i].value);

                break;
            case OPT_REF_OUTPUT_FILE:
                if (opts[i].value == NULL)
                    SR_ErrQuit("ERROR: The output reference file is not specified.\n");

                pars->refOutput = fopen(opts[i].value, "wb");
                if (pars->refOutput == NULL)
                    SR_ErrSys("ERROR: Cannot open reference file \"%s\" for writing.\n", opts[i].value);

                break;
            case OPT_HASH_TABLE_FILE:
                if (opts[i].value == NULL)
                    SR_ErrQuit("ERROR: The output hash table file is not specified.\n");

                pars->hashTableOutput = fopen(opts[i].value, "wb");
                if (pars->hashTableOutput == NULL)
                    SR_ErrSys("ERROR: Cannot open hash table file \"%s\" for writing.\n", opts[i].value);

                break;
            case OPT_HASH_SIZE:
                if (opts[i].value == NULL)
                    SR_ErrQuit("ERROR: Hash size is not specified.\n");

                pars->hashSize = atoi(opts[i].value);
                if (pars->hashSize <= 0 || pars->hashSize >= MAX_HASH_SIZE)
                    SR_ErrQuit("ERROR: Invalid hash size. Hash size should be greater than zero and less than %d.\n", MAX_HASH_SIZE);

                break;
            case OPT_SPECIAL_REF_INPUT:
                if (opts[i].isFound)
                {
                    pars->specialRefInput = fopen(opts[i].value, "r");
                    if (pars->specialRefInput == NULL)
                        SR_ErrSys("ERROR: Cannot open special reference fasta file \"%s\" for reading.\n", opts[i].value);
                }
                else
                    pars->specialRefInput = NULL;

                break;
            default:
                SR_ErrQuit("ERROR: Unrecognized argument.\n");
                break;
        }
    }

    if (optNum < OPT_BUILD_REQUIRED_NUM)
        SR_ErrQuit("ERROR: Incorrect number of arguments.\n");
}

// show the help message and quit
void SR_Build_ShowHelp(void)
{
    printf("Usage: SR_Build -fi <input_fasta_file> -ro <reference_output_file> -hto <hash_table_output_file> -hs <hash_size> [-sfi special_fasta_file]\n");
    printf("Read in the reference file in fasta file and ouput the SR format reference file and hash table file.\n\n");

    printf("-fi       input reference file in fasta format\n");
    printf("-ro       output reference file in \"SR\" format\n");
    printf("-hto      output hash table file.\n");
    printf("-hs       hash size parameter(1 - %d)\n", MAX_HASH_SIZE);
    printf("-sfi      input special reference file in fast format (optional)\n");
    printf("-help     display help message and exit\n\n");

    exit(EXIT_SUCCESS);
}

// clean up the resouses used in the split-read build program
void SR_Build_Clean(SR_Reference* reference, SR_RefHeader* refHeader, SR_OutHashTable* refHashTable, SR_Build_Pars* buildPars)
{
    SR_ReferenceFree(reference);
    SR_RefHeaderFree(refHeader);
    SR_OutHashTableFree(refHashTable);

    fclose(buildPars->faInput);
    fclose(buildPars->refOutput);
    fclose(buildPars->hashTableOutput);

    if (buildPars->specialRefInput != NULL)
        fclose(buildPars->specialRefInput);
}
