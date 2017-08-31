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

#ifndef LOOKUPTABLE5_H_
#define LOOKUPTABLE5_H_

#include "query.h"

void lookupTable5Create_wrap(const char *refFile, int seedLen,
		int maxHitsPerQry, int **lookupTable, char **refIdx, int **refPos,
		int *numDistinctTuples, int **numRepeatsPerTuple,
		int tupleIgnoreThreshold);

void lookupTable5Create(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold, int *totalTuples);

void lookupTable5Create2(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold, int *totalTuples);

void lookupTable5Create(const char *keysFile, const char *valsFile, int maxHits,
		int seedLen, int tupleIgnoreThres);

void lookupTable5CpyHashTableToGPU(int **keys, int *numKeys, int **values,
		int *numValues, int **numRepeatsPerTuple);

void lookupTable5Delete();

int lookupTable5MapQry(const Query *query, uint qryLen, uint isRevComp,
		char *refIdx, int *shift, int *refPos);

int lookupTable5MapQry2(char *qrySeq, int qryLen, char *refIdx_bestHits,
		int *shift_bestHits, int *refPos_bestHits);

__global__ void lookupTable5MapQry_gpu(int *keys, int *values,
		int *numTuplesPerHash, char *qrs, uchar *qryLen, char *refIdx,
		int *refPos, int maxHits, int seedLen, int randNum, int arrSize);

__global__ void lookupTable5MapQry2_gpu(int *keys, int *values,
		int *numRptsPerTuple, char *qrs, uchar *qryLen, short maxQrySeqLen,
		char *refIdx, int *refPos, int maxHits, int seedLen, int randNum,
		int arrSize);

__device__ int arrSearch_gpu(int *arr, int arrSize, int num);

__device__ void arrGetRandomNums_gpu(int n, int lLimit, int uLimit, int *arr);

__device__ int getRandNum_gpu();

__device__ int pow_gpu(int base, int n);

__global__ void getHash_gpu_wrap(char *str, int len, int *hash);

__device__ int getHash_gpu(char *str, int len);

__global__ void sort_gpu_wrap(char *refIdx, int *shift, int *refPos, int arrSize);

__device__ void sort_gpu(char *refIdx, int *shift, int *refPos, int arrSize,
		short threadId);

__global__ void createClusters_gpu_wrap(char *refIdx, int *shift,
		int *clusterSize, int arrSize, short *bgstClustSize, short *numClusters);

__device__ short createClusters_gpu(char *refIdx, int *shift, int *clusterSize,
		int arrSize,  short *bgstClustSize);

__global__ void findBiggestClusters_gpu_wrap(short numClust,
		short bgstClustSize, int *clusterSize, int *bgstClust,
		char *numBgstClust);

__device__ char findBiggestClusters_gpu(short numClust, short bgstClustSize,
		int *clusterSize, int *bgstClust);

__global__ void assignResults_gpu_wrap(char maxHits, char numBgstHits,
		char *refIdx, char *refIdx_global, int *refPos, int *refPos_global,
		int *bgstClust, int randNum);

__device__ void assignResults_gpu(int blockId, char maxHits, char numBgstHits,
		char *refIdx, char *refIdx_global, int *refPos, int *refPos_global,
		int *bgstClust, int *randNums);

__global__ void intializeShrMem_gpu_wrap(char *refIdx, int *shift, int *refPos,
		int *clusterSize, int arrSize);

__device__ void initializeShrMem_gpu(char *refIdx, int *shift, int *refPos,
		int *clusterSize, int arrSize, short threadId);

__global__ void cpyHitsFromGlobalToShr_gpu_wrap(char *refIdx, int *shift,
		int *refPos, int *hashes, int arrSize, int *keys, int *values,
		int *numRptsPerTuple, short numQryTuples, int seedLen);

__device__ void cpyHitsFromGlobalToShr_gpu(char *refIdx, int *shift,
		int *refPos, int *hashes, int arrSize, int *keys, int *values,
		int *numRptsPerTuple, short numQryTuples, int seedLen);

void lookupTable5CpyConstMemToGPU();

#endif /* LOOKUPTABLE5_H_ */
