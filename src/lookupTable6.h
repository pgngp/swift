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

#ifndef LOOKUPTABLE6_H_
#define LOOKUPTABLE6_H_

void lookupTable6Create2(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold, int *totalTuples);

void lookupTable6Reset();

void lookupTable6Delete();

void lookupTable6CpyHashTableToGPU(int **keys, int *numKeys, int **values,
		int *numValues, int **numRepeatsPerTuple);

int lookupTable6MapQry2(char *qrySeq, int qryLen, char *refIdx_bestHits,
		int *shift_bestHits, int *refPos_bestHits);

void lookupTable6CpyConstMemToGPU();

__global__ void lookupTable6MapQry2_gpu(int *keys, int *values,
		int *numRptsPerTuple, char *qrs, uchar *qryLen, short maxQrySeqLen,
		char *refIdx, int *refPos, int maxHits, int seedLen, int randNum,
		int arrSize);

__global__ void cpyHitsFromGlobalToShr_gpu_wrap(long *refPos, int *hashes,
		int arrSize, int *keys, int *values, int *numRptsPerTuple,
		short numQryTuples, int seedLen, short maxQryLen);

__device__ void cpyHitsFromGlobalToShr_gpu(long *refPos, int *hashes,
		int arrSize, int *keys, int *values, int *numRptsPerTuple,
		short numQryTuples, int seedLen, short maxQryLen);

__global__ void intializeShrMem_gpu_wrap(long *refPos, int *clusterSize,
		int arrSize);

__device__ void initializeShrMem_gpu(long *refPos, int *clusterSize, int arrSize,
		short threadId);

__global__ void assignResults_gpu_wrap(char maxHits, char numBgstHits,
		long *refPos, char *refIdx_global, int *refPos_global, short *bgstClust,
		int randNum);

__device__ void assignResults_gpu(int blockId, char maxHits, char numBgstHits,
		long *refPos, char *refIdx_global, int *refPos_global, short *bgstClust,
		int *randNums);

__global__ void findBiggestClusters_gpu_wrap(short numClust,
		short bgstClustSize, int *clusterSize, short *bgstClust,
		char *numBgstClust);

__device__ char findBiggestClusters_gpu(short numClust, short bgstClustSize,
		int *clusterSize, short *bgstClust);

__global__ void createClusters_gpu_wrap(long *refPos, int *clusterSize,
		int arrSize, short *bgstClustSize, short *numClusters);

__device__ short createClusters_gpu(long *refPos, int *clusterSize, int arrSize,
		short *bgstClustSize);

__global__ void sort_gpu_wrap6(long *refPos, int arrSize);

__device__ void sort_gpu6(long *refPos, int arrSize, short threadId);

__global__ void getHash_gpu_wrap6(char *str, int len, int *hash);

__device__ int getHash_gpu6(char *str, int len);

__device__ int arrSearch_gpu6(int *arr, int arrSize, int num);

__device__ void arrGetRandomNums_gpu6(int n, int lLimit, int uLimit, int *arr);

__device__ int getRandNum_gpu6();

__device__ int pow_gpu6(int base, int n);

#endif /* LOOKUPTABLE6_H_ */
