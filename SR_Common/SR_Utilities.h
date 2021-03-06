/*
 * =====================================================================================
 *
 *       Filename:  SR_Utilities->h
 *
 *    Description:  
 *
 *        Version:  1->0
 *        Created:  07/10/2011 03:01:14 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_UTILITIES_H
#define  SR_UTILITIES_H

#include <math.h>
#include <ctype.h>
#include <stdio.h>

#include "SR_Types.h"


#define SR_TO_STRING_PRE(obj) #obj

#define SR_TO_STRING(obj) SR_TO_STRING_PRE(obj)

#define SR_SWAP(objectOne, objectTwo, type)      \
    do                                           \
    {                                            \
        type temp;                               \
        temp = (objectOne);                      \
        (objectOne) = (objectTwo);               \
        (objectTwo) = temp;                      \
                                                 \
    }while(0)


//=======================
// Array utilities
//=======================

#define SR_ARRAY_GET(pArray, i) ((pArray)->data[(i)])

#define SR_ARRAY_GET_PT(pArray, i) ((pArray)->data + (i))

#define SR_ARRAY_GET_SIZE(pArray) ((pArray)->size)

#define SR_ARRAY_GET_FIRST(pArray) ((pArray)->data[0])

#define SR_ARRAY_GET_FIRST_PT(pArray) ((pArray)->data)

#define SR_ARRAY_GET_LAST(pArray) ((pArray)->data[(pArray)->size - 1])

#define SR_ARRAY_GET_LAST_PT(pArray) ((pArray)->data + (pArray)->size - 1)

#define SR_ARRAY_ALLOC(pArray, pArrayCap, objType, dataType)              \
    do                                                                    \
    {                                                                     \
        (pArray) = (objType*) malloc(sizeof(objType));                    \
        if ((pArray) == NULL)                                             \
            SR_ErrSys("ERROR: not enough memory for an object->\n");      \
                                                                          \
        SR_ARRAY_INIT((pArray), pArrayCap, dataType);                     \
                                                                          \
    }while(0)

#define SR_ARRAY_INIT(pArray, pArrayCap, dataType)                           \
    do                                                                       \
    {                                                                        \
        (pArray)->size = SR_EMPTY;                                           \
        (pArray)->capacity = pArrayCap;                                      \
                                                                             \
        (pArray)->data = (dataType*) calloc(pArrayCap, sizeof(dataType));    \
        if ((pArray)->data == NULL)                                          \
            SR_ErrSys("ERROR: not enough memory for the data storage->\n");  \
                                                                             \
    }while(0)


#define SR_ARRAY_PUSH(pArray, pNewElt, dataType)                                                              \
    do                                                                                                        \
    {                                                                                                         \
        if ((pArray)->size == (pArray)->capacity)                                                             \
        {                                                                                                     \
            (pArray)->capacity *= 2;                                                                          \
            (pArray)->data = (dataType*) realloc((pArray)->data, sizeof(dataType) * (pArray)->capacity);      \
            if ((pArray)->data == NULL)                                                                       \
                SR_ErrSys("ERROR: not enough memory to expand the storage of data in an array.\n");           \
        }                                                                                                     \
                                                                                                              \
        SR_ARRAY_GET(pArray, (pArray)->size) = *(pNewElt);                                                    \
        ++((pArray)->size);                                                                                   \
                                                                                                              \
    }while(0)

#define SR_ARRAY_RESIZE(pArray, newCapacity, dataType)                                                        \
    do                                                                                                        \
    {                                                                                                         \
        if ((newCapacity) > (pArray)->capacity)                                                               \
        {                                                                                                     \
            (pArray)->capacity = (newCapacity);                                                               \
            (pArray)->data = (dataType*) realloc((pArray)->data, sizeof(dataType) * (pArray)->capacity);      \
            if ((pArray)->data == NULL)                                                                       \
                SR_ErrSys("ERROR: not enough memory to expand the storage of data in an array.\n");           \
                                                                                                              \
        }                                                                                                     \
                                                                                                              \
    }while(0)

#define SR_ARRAY_RESIZE_NO_COPY(pArray, newCapacity, dataType)                                                \
    do                                                                                                        \
    {                                                                                                         \
        if ((newCapacity) > (pArray)->capacity)                                                               \
        {                                                                                                     \
            free(pArray->data);                                                                               \
            (pArray)->capacity = newCapacity;                                                                 \
            (pArray)->data = (dataType*) calloc(sizeof(dataType), (pArray)->capacity);                        \
            if ((pArray)->data == NULL)                                                                       \
                SR_ErrSys("ERROR: not enough memory to expand the storage of data in an array.\n");           \
                                                                                                              \
        }                                                                                                     \
                                                                                                              \
        (pArray)->size = 0;                                                                                   \
                                                                                                              \
    }while(0)

#define SR_ARRAY_POP(pArray)   \
    do                         \
    {                          \
        (pArray)->size -= 1;   \
                               \
    }while(0)


#define SR_ARRAY_RESET(pArray) (pArray)->size = SR_EMPTY


#define SR_ARRAY_FREE(pArray, freeObj)   \
    do                                   \
    {                                    \
        if ((pArray) != NULL)            \
        {                                \
            free((pArray)->data);        \
            if (freeObj)                 \
                free(pArray);            \
            else                         \
                (pArray)->data = NULL;   \
        }                                \
                                         \
    }while(0)


#define SR_ARRAY_IS_FULL(pArray) ((pArray)->size == (pArray)->capacity)

int FindKthSmallestInt(int array[], int size, int k);

int FindMedianInt(int array[], int size);

unsigned int FindKthSmallestUint(unsigned int array[], unsigned int size, unsigned int k);

unsigned int FindMedianUint(unsigned int array[], unsigned int size);

static inline int DoubleRoundToInt(double number)
{
    return ((int) floor(number + 0.5));
}

static inline void StrToUpper(char* str)
{
    for (unsigned int i = 0; str[i] != '\0'; ++i)
        str[i] = toupper(str[i]);
}

char* SR_CreateFileName(const char* workingDir, const char* fileName);

SR_Status SR_GetNextLine(char* buff, unsigned int buffSize, FILE* input);


#endif  /*SR_UTILITIES_H*/
