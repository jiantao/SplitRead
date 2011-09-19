/*
 * =====================================================================================
 *
 *       Filename:  SR_Types.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/08/2011 21:43:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_TYPES_H
#define  SR_TYPES_H

typedef enum
{
    FALSE = 0,
    TRUE  = 1

}SR_Bool;

typedef enum
{
    SR_OK  =           0,
    SR_EOF =          -1,
    SR_ERR =          -2,
    SR_FULL      =    -98,
    SR_NOT_FOUND =    -99,
    SR_OUT_OF_RANGE = -100

}SR_Status;

typedef enum
{
    SR_A = 1,
    SR_C = 2,
    SR_G = 4,
    SR_T = 8,
    SR_N = 15

}SR_Base;

// search direction for query region relative to the positon of anchor mate
typedef enum
{
    SR_UPSTREAM,     // search the query region uptream to the anchor mate (smaller coordinate)

    SR_DOWNSTREAM    // search the query region downstream to the anchor mate (larger coordinate)

}SR_Direction;

// action applied on a certain DNA sequence
typedef enum
{
    SR_FORWARD,

    SR_REVERSE_COMP,    // change the sequence to its reverse complement format

    SR_INVERSE,         // change the sequence to its inverse format

    SR_COMP

}SR_SeqAction;

typedef enum
{
    SR_1F = 0,
    SR_1R = 1,
    SR_2F = 2,
    SR_2R = 3

}SR_SingleMode;

typedef enum
{
    SR_1F2F = ((SR_1F << 2) | SR_2F),
    SR_1F2R = ((SR_1F << 2) | SR_2R),
    SR_1R2F = ((SR_1R << 2) | SR_2F),
    SR_1R2R = ((SR_1R << 2) | SR_2R),
    SR_2F1F = ((SR_2F << 2) | SR_1F),
    SR_2F1R = ((SR_2F << 2) | SR_1R),
    SR_2R1F = ((SR_2R << 2) | SR_1F),
    SR_2R1R = ((SR_2R << 2) | SR_1R)

}SR_PairMode;

typedef SR_SeqAction SR_Strand;

#define MD5_CHECKSUM_LEN 16

#define MD5_STR_LEN 32

#define MAX_HASH_SIZE 12

#define SR_EMPTY 0

#endif  /*SR_TYPES_H*/
