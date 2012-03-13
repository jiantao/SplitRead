/*
 * =====================================================================================
 *
 *       Filename:  SR_FragLenHist.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/08/2012 02:24:19 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef  SR_FRAGLENHIST_H
#define  SR_FRAGLENHIST_H

#include <stdint.h>
#include <stdio.h>

#include "SR_Types.h"

#define INVALID_FRAG_LEN_QUAL 255

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

    uint32_t capacity;              // capacity of the histogram

    uint64_t modeCount[2];          // total counts of a histogram

}SR_FragLenHist;

typedef struct SR_FragLenHistArray
{
    SR_FragLenHist* data;

    uint32_t size;

    uint32_t capacity;

}SR_FragLenHistArray;

SR_FragLenHistArray* SR_FragLenHistArrayAlloc(unsigned int capacity);

void SR_FragLenHistArrayFree(SR_FragLenHistArray* pHistArray);

void SR_FragLenHistArrayClear(SR_FragLenHistArray* pHistArray);

void SR_FragLenHistArrayInit(SR_FragLenHistArray* pHistArray, unsigned int size);

SR_Status SR_FragLenHistArrayUpdate(SR_FragLenHistArray* pHistArray, unsigned int backHistIndex, uint32_t fragLen);

int SR_FragLenHistArrayGetFragLenQual(const SR_FragLenHistArray* pHistArray, unsigned int backHistIndex, uint32_t fragLen);

void SR_FragLenHistArrayFinalize(SR_FragLenHistArray* pHistArray);

void SR_FragLenHistArrayWrite(const SR_FragLenHistArray* pHistArray, FILE* output);

#endif  /*SR_FRAGLENHIST_H*/
