/*
 * =====================================================================================
 *
 *       Filename:  SR_ReadPairBuild.h
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

#ifndef  SR_READPAIRBUILD_H
#define  SR_READPAIRBUILD_H

#include "SR_LibInfo.h"


//===============================
// Type and constant definition
//===============================

// parameters used for read pair build
typedef struct
{
    double cutoff;           // fragment length cutoff p-value

    double trimRate;         // trim rate from fragment length distribution

    unsigned int binLen;     // bin length for searching mate

    uint32_t detectSet;      // bit set used to store the SV types that should be detected 

    unsigned char minMQ;     // minimum mapping quality for a read pair

    FILE* fileListInput;     // input stream of a file list containing all the bam file names

    char* workingDir;        // working directory for the detector

}SR_ReadPairBuildPars;

// local pair structure(for deletion, tademn duplication and inversion)
typedef struct SR_LocalPairInfo
{
    int32_t readGrpID;                                        // read group ID

    int32_t refID:16, upMapQ:8, downMapQ:8;                   // reference ID, up mate mapping quality, down mate mapping quality

    int32_t upPos;                                            // alignment position of the up mate

    int32_t upEnd;                                            // alignment end of the up mate

    int32_t downPos;                                          // alignment position of the down mate

    int32_t fragLen;                                          // fragment length of the read pair

    int32_t upNumMM:16, downNumMM:16;                         // number of mismatches of the up mate, number of mismatches of the down mate

    int32_t fragLenQual:16, readPairType:8, pairMode:8;       // fragment length quality(-10log(p-value)), read pairt type, pair mode

}SR_LocalPair;

// cross pair structure(two mate hit different reference uniquely)
typedef struct SR_CrossPair
{
    int32_t readGrpID;                                            // read group ID
                                                                                                                                          
    int32_t upRefID:16, downRefID:16;                             // up reference ID, down reference ID
                                                                                                                                          
    int32_t upMapQ:16, downMapQ:16;                               // up mapping quality, down mapping quality
                                                                  
    int32_t upPos;                                                // alignment position of the up mate
                                                                                                                                          
    int32_t upEnd;                                                // alignment end of the up mate
                                                                                                                                          
    int32_t downPos;                                              // alignment position of the down mate
                                                                                                                                          
    int32_t downEnd;                                              // alignment end of the down mate
                                                                                                                                          
    int32_t upNumMM:16, downNumMM:16;                             // number of mismatches of the up mate, number of mismatches of the down
                                                                                                                                          
    int32_t fragLenQual:16, readPairType:8, pairMode:8;           // fragment length quality(-10log(p-value)), read pairt type, pair mode

}SR_CrossPair;

// special pair structure (one mate is unique the other mate hits the special reference)
typedef struct SR_SpecialPair
{
    int32_t readGrpID;                                           // read group ID
                                                                                                                                         
    int16_t refID[2];                                            // up reference ID, down reference ID
                                                                                                                                         
    int32_t pos[2];                                              // alignment position of the up and down mate
                                                                 
    int32_t end[2];                                              // alignment end of the up and down mate

    int32_t numSpeicalHits;                                      // number of hits on the special reference
                                                                                                                                         
    int16_t numMM[2];                                            // number of mismatches of the up mate, number of mismatches of the down
                                                                                                                                         
    uint8_t bestMQ[2];                                           // best mapping quality of the up and down mate
                                                                                                                                         
    uint8_t secMQ[2];                                            // second best mapping quality of the up and down mate

    int32_t fragLenQual:16, readPairType:8, pairMode:8;          // fragment length quality(-10log(p-value)), read pairt type, pair mode

    uint32_t specialID;                                          // special reference ID

}SR_SpecialPair;

// local pair array
typedef struct SR_LocalPairArray
{
    SR_LocalPair* data;     // local pairs

    uint64_t* chrCount;     // number of pairs in each chromosome

    uint64_t size;          // numer of local pairs

    uint64_t capacity;      // capacity of the array

}SR_LocalPairArray;

// cross pair array
typedef struct SR_CrossPairArray
{
    SR_CrossPair* data;     // cross pairs

    uint64_t* chrCount;     // number of pairs in each chromosomes(up mate)

    uint64_t size;          // number of cross pairs

    uint64_t capacity;      // capacity of the array

}SR_CrossPairArray;

// special pair array
typedef struct SR_SpecialPairArray
{
    SR_SpecialPair* data;    // special pairs

    uint64_t* chrCount;      // number of pairs in each chromosomes(anchor mate)

    uint64_t size;           // number of special pairs

    uint64_t capacity;       // capacity of the array

}SR_SpecialPairArray;

// special pair table
typedef struct SR_SpecialPairTable
{
    char (*names)[3];                    // array of special reference names(from ZA tags)

    void* nameHash;                      // special reference name hash

    uint32_t size;                       // number of the special reference names

    uint32_t capacity;                   // capacity of the name array

    SR_SpecialPairArray array;           // array for local special pairs

    SR_SpecialPairArray crossArray;      // array for cross special pairs

}SR_SpecialPairTable;

// read pair table
typedef struct SR_ReadPairTable
{
    SR_LocalPairArray* pLongPairArray;          // long pair array(deletion)

    SR_LocalPairArray* pShortPairArray;         // short pair array(tadem dup)

    SR_LocalPairArray* pReversedPairArray;      // reveresed pair array (tadem dup)

    SR_LocalPairArray* pInvertedPairArray;      // inverted pair array (inversion)

    SR_CrossPairArray* pCrossPairArray;         // cross pair array (inter-chromosome translocation)

    SR_SpecialPairTable* pSpecialPairTable;     // special pair table (retro insertion)

    uint32_t detectSet;                         // SV detect set

    uint32_t numChr;                            // number of references

    uint64_t numPairs;                          // total number of pairs in the read pair table

}SR_ReadPairTable;


//===============================
// Constructors and Destructors
//===============================

SR_SpecialPairTable* SR_SpecialPairTableAlloc(unsigned int numChr);

void SR_SpecialPairTableFree(SR_SpecialPairTable* pSpecialPairTable);

SR_ReadPairTable* SR_ReadPairTableAlloc(uint32_t numChr, uint32_t detectSet);

void SR_ReadPairTableFree(SR_ReadPairTable* pInfoTable);


//======================
// Interface functions
//======================

//===============================================================
// function:
//      read the bam alignments from the primary bam files, 
//      extract the SV-related information out and wirte it into 
//      read pair files according to the read pair type and 
//      location
//
// args:
//      1. pBuildPars: a pointer to the struct of parameters 
//      for read pair build
// 
//=============================================================== 
void SR_ReadPairBuild(const SR_ReadPairBuildPars* pBuildPars);

//===============================================================
// function:
//      clear the read pair table
//
// args:
//      1. pReadPairTable: a pointer to the read pair table
//                         structure
// 
//=============================================================== 
void SR_ReadPairTableClear(SR_ReadPairTable* pReadPairTable);

//======================================================================
// function:
//      update the read pair table with the incoming read pairs
//
// args:
//      1. pReadPairTable: a pointer to the read pair table structure
//      2. pUpAlgn: a pointor to the bam alignment of the up mate
//      3. pDownAlgn: a pointor to the bam alignment of the down mate
//      4. pZAtag: a pointer to the ZA tag structure
//      5. ppairStats: a pointer to the pair stats structure
//      6. pLibTable: a pointer to the library information table
//      7. pHistArray: a pointer to the fragment length histogram array
//      8. minMQ: minimum mapping quality for a read pair
//======================================================================
void SR_ReadPairTableUpdate(SR_ReadPairTable* pReadPairTable, const bam1_t* pUpAlgn, const bam1_t* pDownAlgn, const SR_ZAtag* pZAtag,
                            const SR_PairStats* pPairStats, const SR_LibInfoTable* pLibTable, const SR_FragLenHistArray* pHistArray, uint8_t minMQ);

//================================================================
// function:
//      open a seriers read pair files for output. A read pair
//      file will be open for each read pair type on each 
//      chromosome
//
// args:
//      1. pLibTable: a pointer to the library information table
//      2. detectSet: a bit set indicating which SV event shoud
//                    be detected
//      3. workingDir: the working directory for the detector
//
// return:
//      a pointer ot a hash table of read pair files
//================================================================
void* SR_ReadPairFilesOpen(const SR_LibInfoTable* pLibTable, uint32_t detectSet, const char* workingDir);

//==================================================================
// function:
//      close all the read pair files
//
// args:
//      1. pFileHash: a pointer to a hash table of read pair files
//==================================================================
void SR_ReadPairFilesClose(void* pFileHash);

//=================================================================
// function:
//      write the read pair table into the read pair files
//
// args:
//      1. pReadPairTable: a pointer to a read pair table
//      2. pFileHash: a pointer to a hash table of read pair files
//=================================================================
void SR_ReadPairTableWrite(const SR_ReadPairTable* pReadPairTable, void* pFileHash);

//=================================================================
// function:
//      write the special reference name into the end of the
//      library information file
//
// args:
//      1. pSpecialTable: a pointer to a special pair table
//      2. libOutput: a file pointer to a library information file
//=================================================================
void SR_SpecialPairTableWirteID(const SR_SpecialPairTable* pSpecialTable, FILE* libOutput);

//=================================================================
// function:
//      write the read pair table into the read pair files
//
// args:
//      1. pReadPairTable: a pointer to a read pair table
//      2. pFileHash: a pointer to a hash table of read pair files
//=================================================================
SR_Status SR_CheckWorkingDir(const char* workingDir);

#endif  /*SR_READPAIRBUILD_H*/
