/*
 * =====================================================================================
 *
 *       Filename:  SR_Build_GetOpt.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/28/2011 05:07:49 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_BUILD_GETOPT_H
#define  SR_BUILD_GETOPT_H

#include <stdio.h>
#include "SR_Types.h"
#include "SR_Reference.h"
#include "SR_OutHashTable.h"


// an object hold the options parsed from command line arguments
typedef struct SR_Option
{
    char* name;          // name of this option, following a '-' character in the command line

    const char* value;   // a pointer points to the value of the option in the argv array

    SR_Bool isFound;        // boolean variable to indicate that if we fount this option or not

}SR_Option;

// an object hold the parameters used in the split-read build program
typedef struct SR_Build_Pars
{
    FILE* faInput;            // input stream of the fasta file

    FILE* refOutput;          // output stream of the reference file

    FILE* hashTableOutput;    // output stream of the hash table file

    unsigned char hashSize;   // hash size used to index the reference

}SR_Build_Pars;

// get the options from command line arguemnts
int SR_GetOpt(SR_Option opts[], int argc, char* argv[]);

// set the parameters for the split-read build program from the parsed command line arguments 
void SR_Build_SetPars(SR_Build_Pars* pars, int argc, char* argv[]);

// show the help message and quit
void SR_Build_ShowHelp(void);

// clean up the resouses used in the split-read build program
void SR_Build_Clean(SR_Reference* reference, SR_RefHeader* refHeader, SR_OutHashTable* refHashTable, SR_Build_Pars* buildPars);

#endif  /*SR_BUILD_GETOPT_H*/
