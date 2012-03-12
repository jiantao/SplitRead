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

#include "bam.h"
#include "SR_Types.h"
#include "SR_BamHeader.h"
// #include "SR_BamPairAux.h"

//===============================
// Type and constant definition
//===============================

#define SR_IsValidPairMode(pDstrb, pairMode) ((pDstrb)->validModeMap & (1 << (pairMode)))

#define SR_GetHistMedian(pDstrb, readGrpID) ((pDstrb)->pHists[(readGrpID)].median)

// a map used to map the pair mode into its corresponding number
// negative value means invalid mode
static const int SR_PairModeMap[16] = { 
                                          -1, -1, 0, 1,

                                          -1, -1, 2, 3,

                                          4, 5, -1, -1,

                                          6, 7, -1, -1
                                      };


static const int SR_PairModeSetMap[64] =  { 
                                               0, 0, 0, 1, 0, 0, 0, 1,
                                               0, 0, 1, 0, 0, 1, 0, 0,
                                               0, 1, 0, 0, 0, 0, 1, 0,
                                               1, 0, 0, 0, 1, 0, 0, 0,
                                               0, 0, 0, 1, 0, 0, 0, 1,
                                               0, 1, 0, 0, 0, 0, 1, 0,
                                               0, 0, 1, 0, 0, 1, 0, 0,
                                               1, 0, 0, 0, 1, 0, 0, 0
                                          };
typedef enum
{
    ST_ILLUMINA = 0,

    ST_454 = 1,

    ST_SOLID = 2,

    ST_ILLUMINA_LONG = 3

}SR_SeqTech;


// the object used to hold the basic statistics of a pair of alignments
typedef struct SR_BamPairStats
{
    const char* RG;           // name of the read group

    unsigned int fragLen;     // fragment length of the pair

    int16_t pairMode;        // orientation mode of the pair

}SR_BamPairStats;

// the object used to hold the fragment length histogram of a given read group
typedef struct SR_FragLenHist
{
    void* rawHist;                  // raw fragment length histogram. this is a hash table for histogram building

    uint32_t* fragLen;              // array of the fragment length

    uint32_t* freq;                 // frequency of each fragment length

    double mean;                    // mean of the histogram

    double median;                  // median of the histogram

    double stdev;                   // standard deviation of the histogram

    uint32_t size;                  // number of unique fragment length

    uint64_t modeCount[2];          // total counts of a histogram

}SR_FragLenHist;

typedef struct SR_FragLenHistArray
{
    SR_FragLenHist* data;

    uint32_t size;

    uint32_t capacity;

}SR_FragLenHistArray;

// the object used to hold the fragment length distribution of all the read groups in a bam file
typedef struct SR_FragLenDstrb
{
    char** pSamples;

    void* pSampleHash;

    char** pReadGrps;

    void* pReadGrpHash;

    int8_t* pSeqTech;

    int32_t* pSampleMap;

    SR_FragLenHist* pHists;                        // array of fragment length histograms

    uint32_t sizeSM;

    uint32_t capacitySM;
    
    uint32_t sizeRG;                                 // number of read groups found in the bam file

    uint32_t capacityRG;                             // capacity of the read group name array

    uint32_t sizeHist;

    uint32_t capacityHist;

    int8_t validModeMap;                           // map the pair mode to its corresponding histogram, invalid pair mode will get a negative value
    
    int8_t validMode[NUM_ALLOWED_PAIR_MODE];       // the valid pair modes

}SR_FragLenDstrb;



//===============================
// Constructors and Destructors
//===============================

SR_FragLenDstrb* SR_FragLenDstrbAlloc(void);

void SR_FragLenDstrbFree(SR_FragLenDstrb* pDstrb);


//======================
// Inline functions
//======================

//==================================================================
// function:
//      get the pair orientation mode from a bam alignment
//
// args:
//      1. pBamNode: a pointer to a bam node structure
//
// return:
//      the mode of the pair
//==================================================================
SR_PairMode SR_GetPairMode(const bam1_t* pAlignment);

//=======================================================================
// function:
//      get the basic statistics from a pair of alignments
//
// args:
//      1. pPairStats: a poiter to the pair statistics object
//      2. ppUpAlgn: a pointer to a bam node object for the alignment 
//                   with smaller coordinate
//      3. ppDownAlgn: a pointer to a bam node object for the alignment 
//                     with greater coordinate
//
// return:
//      if an error occurs, return SR_ERR; else return SR_OK
//========================================================================
SR_Status SR_LoadPairStats(SR_BamPairStats* pPairStats, const bam1_t* pAlignment);


void SR_FragLenDstrbSetPairMode(SR_FragLenDstrb* pDstrb, const int8_t* pValidPairMode);

//=======================================================================
// function:
//      retrieve the read group information from the bam header and
//      initialize the fragment length distribution with it
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. pBamHeader: a pointer to a bam header object
//
// return:
//      if an error occurs, return SR_ERR; else return SR_OK
//========================================================================
SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader);

//=======================================================================
// function:
//      given a read group name, find its corresponding read group ID
//
// args:
//      1. pReadGrpIndex: a pointer to the read group index variable
//                        this is the return value
//      2. pDstrb: a pointer to a fragment length distribution object
//      3. pReadGrpName: a pointer to the read group name
//
// return:
//      if an error occurs, return SR_ERR; else return SR_OK
//========================================================================
SR_Status SR_FragLenDstrbGetRGIndex(int32_t* pReadGrpIndex, const SR_FragLenDstrb* pDstrb, const char* pRreadGrpName);

//=======================================================================
// function:
//      update the fragment length distribution with the information
//      from a new read pair
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. pPairStats: a pointer to a read pair statistic object
//
// return:
//      if an error occurs, return SR_ERR; else return SR_OK
//========================================================================
SR_Status SR_FragLenDstrbUpdate(SR_FragLenDstrb* pDstrb, const SR_BamPairStats* pPairStats);

//=======================================================================
// function:
//      after all the read pairs are processed, raw histogram will be
//      transferred into mature one. its mean, median and standard
//      deviation will be calculated.
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//========================================================================
void SR_FragLenDstrbFinalize(SR_FragLenDstrb* pDstrb);

//=======================================================================
// function:
//      set the probability cutoff for all the histograms
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. cutoff: the fragment length probability cutoff
//========================================================================

void SR_FragLenDstrbInitHistHeader(const SR_FragLenDstrb* pDstrb, FILE* histOutput);

void SR_FragLenDstrbWriteHistHeader(uint32_t numHists, FILE* histOutput);
//=======================================================================
// function:
//      write the fragment length distribution into a file
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. dstrbOutput: output stream
//========================================================================
void SR_FragLenDstrbWriteHist(const SR_FragLenDstrb* pDstrb, FILE* histOutput);


void SR_FragLenDstrbWriteInfo(const SR_FragLenDstrb* pDstrb, FILE* infoOutput);

//=======================================================================
// function:
//      read the fragment length distribution from a file
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. dstrbOutput: input stream
//========================================================================
void SR_FragLenDstrbRead(SR_FragLenDstrb* pDstrb, FILE* dstrbInput);

#endif  /*SR_FRAGLENDSTRB_H*/
