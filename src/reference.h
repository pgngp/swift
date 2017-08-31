/*
 * Copyright 2012, St. Jude Children's Research Hospital.
 * Written by Pankaj Gupta, pankaj.gupta@stjude.org.
 *
 * This file is part of Swift.  Swift is free software:  you can redistribute
 * it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 2 of the License,
 * or (at your option) any later version.
 *
 * Swift is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Swift.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include "refPosList.h"
#include "hitList.h"
#include "query.h"


RefPosList **refPreprocess(const char *refFile, uint seedLength);

void refPreprocess2(const char *refFile, uint seedLength);

HitList *refSearchQuery(RefPosList **hashTable, const Query *query,
		uint seqLength, uint isReverseComplement, uint seedLength,
		uint maxNumHits);

void refSearchQuery_paired(RefPosList **hashTable, const Query *qry1,
		const Query *qry2, uint seqLen, uint seedLen, uint maxNumHits,
		uint minFragSize, uint maxFragSize, HitList **result);

__global__ void refSearchQuery_kernel(long *keys, long *vals, char *qrySeq,
		long *results, uint qryLen, uint seedLength, uint maxNumHits);

__global__ void refSearchQuery2_kernel(uint *keys, uint *vals, int numKeys,
		int numVals, char *qrySeq, char *refIdxArr2, int *shiftArr2,
		int *posArr2, int qryLen, int seedLength, int maxNumHits, int randNum);

__global__ void bubbleSort_wrap(long *arr, uint size);

__device__ void bubbleSort(long *arr, uint size);

__global__ void bubbleSort2_wrap(char *refIdxArr, int *shiftArr, int *posArr,
		uint size);

__device__ void bubbleSort2(char *refIdxArr, int *shiftArr, int *posArr,
		uint size);

__global__ void getBestMatches_wrap(long *arr, uint size, uint qryLen,
		long *matches, int maxNumMatches);

__device__ void getBestMatches(long *arr, uint size, uint qryLen,
		long *matches, int maxNumMatches);

__global__ void getBestMatches2_wrap(char *refIdxArr, int *shiftArr,
		int *posArr, int size, int maxNumMatches, char *refIdxArr2,
		int *shiftArr2, int *posArr2);

__device__ void getBestMatches2(char *refIdxArr, int *shiftArr, int *posArr,
		int size, int maxNumMatches, char *refIdxArr2, int *shiftArr2,
		int *posArr2);

static HitList *getBestScoringPosition(HitList *list, uint maxNumHits);

static HitList *getBestScoringPosition2(HitList *list, int maxNumHits);

static int compare(const void *a, const void *b);

static void getHits(RefPosList **hashTable, const Query *qry1, int isQry1RevComp,
		const Query *qry2, int isQry2RevComp, uint seqLen, uint seedLen,
		int maxNumHits, uint minFragSize, uint maxFragSize, HitList **result1,
		HitList **result2);

static void serializeHashTable(RefPosList **hashTable, int numKeys,
		uint numPos, long *keys, long *vals);

void serializeHashTable_wrap(RefPosList **hashTable, int numKeys,
		uint numPos, long *keys, long *vals);

void serializeHashTable2_wrap(RefPosList **hashTable, uint numKeys,
		uint numPos, uint *keys, uint *vals, uint seedLen,
		uint tupleIgnoreThres);

static void serializeHashTable2(RefPosList **hashTable, uint numKeys,
		uint numPos, uint *keys, uint *vals, uint seedLen,
		uint tupleIgnoreThres);

__global__ void getHash_gpu2_wrap(char *s, int length, int *hash);

__device__ static int getHash_gpu2(char *s, int length);

__device__ void arrGetRandomNums_gpu(uint n, uint lowerLimit, uint upperLimit,
		uint *arr);

__global__ void arrGetRandomNums_gpu_wrap(uint n, uint lowerLimit,
		uint upperLimit, uint *arr);

__global__ void arrSearch_gpu_wrap(uint *arr, int arrSize, uint num,
		int *index);

__device__ int arrSearch_gpu(uint *arr, int arrSize, uint num);

__device__ uint getRandNum();

__global__ void getRandNum_wrap(uint *arr, int size);

void refDeleteHashTable(RefPosList **hashTable, uint numKeys);

__global__ void strncpy_gpu_wrap(char *s, char *t, int n);

__device__ void strncpy_gpu(char *s, char *t, int n);

__global__ void pow_gpu_wrap(int base, int n, int *p);

__device__ int pow_gpu2(int base, int n);

void refPrintSerialHashTable(uint *keys, uint numKeys, uint *vals, uint numVals,
		int seedLen);

void refPrintHashTable(RefPosList **hashTable, int numKeys, int seedLen);

__device__ void swap_int(int *arr, int i, int j);

__device__ void swap_char(char *arr, int i, int j);

__global__ void quickSort2_gpu_wrap(char *refIdxArr, int *shiftArr, int *posArr,
		int arrSize);

__device__ int quickSort2_gpu(char *refIdxArr, int *shiftArr, int *posArr,
		int arrSize);

__global__ void insertionSort3_gpu_wrap(char *refIdxArr, int *shiftArr,
		int *posArr, int arrSize);

__device__ void insertionSort3_gpu(char *refIdxArr, int *shiftArr, int *posArr,
		int arrSize);

__device__ void insertionSort4_gpu(char *refIdxArr, int *shiftArr, int *posArr,
		int *idxArr, int arrSize);

#endif /* REFERENCE_H_ */
