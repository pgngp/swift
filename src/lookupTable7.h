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

#ifndef LOOKUPTABLE7_H_
#define LOOKUPTABLE7_H_

void lookupTable7Create2(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold, int *totalTuples);

void lookupTable7Create3(const char *keysFile, const char *valsFile,
		int maxHits, int seedLen, int tupleIgnoreThres);

void lookupTable7Reset();

void lookupTable7Delete();

void lookupTable7CpyHashTableToGPU(int **keys, int *numKeys, int **values,
		int *numValues, int **numRepeatsPerTuple);

void lookupTable7FetchHashTable(int **keys, int *numKeys, int **values,
		int *numValues, int **numRepeatsPerTuple);

int lookupTable7MapQry2(char *qrySeq, int qryLen, char *refIdx_bestHits,
		int *shift_bestHits, int *refPos_bestHits);

void lookupTable7CpyConstMemToGPU();

void lookupTable7CreateConstMem();

__global__ void lookupTable7MapQry2_gpu2(int *keys, int *values,
		int *numRptsPerTuple, char *qrs, uchar *qryLen, short maxQrySeqLen,
		char *refIdx, int *refPos, int maxHits, int seedLen, int randNum,
		int arrSize, int maxNumChr);

void lookupTable7MapQry2_cpu2(int *keys, int *values, int *numRptsPerTuple,
		char *qrs, uchar *qryLen, short maxQrySeqLen, char *refIdx,
		int *refPos, int maxHits, int seedLen, int randNum, int arrSize,
		int maxNumChr, int blockId);

__global__ void cpyHitsFromGlobalToShr_gpu_wrap7(long *refPos, int *hashes,
		int arrSize, int *keys, int *values, int *numRptsPerTuple,
		short numQryTuples, int seedLen, short maxQryLen);

__device__ void cpyHitsFromGlobalToShr_gpu7(long *refPos, int *hashes,
		int arrSize, int *keys, int *values, int *numRptsPerTuple,
		short numQryTuples, int seedLen, short maxQryLen);

__global__ void intializeShrMem_gpu_wrap7(long *refPos, int *clusterSize,
		int arrSize);

__device__ void initializeShrMem_gpu7(long *refPos, int *clusterSize, int arrSize,
		short threadId);

__global__ void assignResults_gpu_wrap7(char maxHits, char numBgstHits,
		long *refPos, char *refIdx_global, int *refPos_global, short *bgstClust,
		int randNum);

__device__ void assignResults_gpu7(int blockId, char maxHits, char numBgstHits,
		long *refPos, char *refIdx_global, int *refPos_global, short *bgstClust,
		int *randNums);

__global__ void findBiggestClusters_gpu_wrap7(short numClust,
		short bgstClustSize, int *clusterSize, short *bgstClust,
		char *numBgstClust);

__device__ char findBiggestClusters_gpu7(short numClust, short bgstClustSize,
		int *clusterSize, short *bgstClust);

__global__ void createClusters_gpu_wrap7(long *refPos, int *clusterSize,
		int arrSize, short *bgstClustSize, short *numClusters);

__device__ short createClusters_gpu7(long *refPos, int *clusterSize, int arrSize,
		short *bgstClustSize);

__global__ void createClusters_gpu7_2_wrap(int *shift, int *refPos, short start,
		short end, short *clusterSize, int *distinctShift, int arrSize);

__device__ void createClusters_gpu7_2(int *shift, int *refPos, short start,
		short end, short *clusterSize, int *distinctShift);

void createClusters_cpu7_2(int *shift, int *refPos, short start, short end,
		short *clusterSize, int *distinctShift);

__global__ void sort_gpu_wrap7(long *refPos, int arrSize);

__device__ void sort_gpu7(long *refPos, int arrSize, short threadId);

__device__ int quickSort_gpu7(char *refIdxArr, int *shiftArr, int *posArr,
		int arrSize, short *begin, short *end);

__global__ void getHash_gpu_wrap7(char *str, int len, int *hash);

__device__ int getHash_gpu7(char *str, int len);

int getHash_cpu7(char *str, int len);

__device__ int arrSearch_gpu7(int *arr, int arrSize, int num);

int arrSearch_cpu7(int *arr, int arrSize, int num);

__device__ void arrGetRandomNums_gpu7(int n, int lLimit, int uLimit, int *arr);

void arrGetRandomNums_cpu7(int n, int lLimit, int uLimit, int *arr,
		uint randNumSeed);

__device__ int getRandNum_gpu7();

int getRandNum_cpu7(uint randNumSeed);

__device__ int pow_gpu7(int base, int n);

__device__ short getBiggestClustSize(short *clustArr, short arrSize);

__device__ short getBiggestHits(short *clustArr, short arrSize,
		short bgstClustSize, int *idxArr, short *bstMtchsArr);

__device__ void assignRsltsToGlobMem(char *refIdx_glb, int *refPos_glb,
		char *refIdx_shr, int *refPos_shr, short numBgstHits,
		short *bestMatches_shr, int blockId, int maxHits);

#endif /* LOOKUPTABLE7_H_ */
