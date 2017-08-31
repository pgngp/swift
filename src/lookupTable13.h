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

#ifndef LOOKUPTABLE13_H_
#define LOOKUPTABLE13_H_

void lookupTable13Create(const char *refFile, int seedLen, int tupleIgnoreThres,
		int *totalTuples);

void lookupTable13Create2(const char *keysFile, const char *valsFile,
		const char *refPosHashMapFile, int maxHits, int seedLen,
		int tupleIgnoreThres, int *numTotalTuples);

void lookupTable13CpyHashTableToGPU(int **keys, int *numKeys, int **values,
		int *numValues, int **numRepeatsPerTuple, int **refPosHashMap,
		int *maxRefPosComposite);

__global__ void lookupTable13MapQry_gpu(int *keys, int *vals,
		int *numTuplesPerHash, int *refPosHashMap, char *qrs, uchar *qryLen,
		int *qryDistance, int maxQryLen, char *refIdx, int *refPos,
		int maxNumHits, int seedLen, int randNum, int idealArrSize);

__global__ void lookupTable13MapQry2_gpu(int *keys, int *vals,
		int *numTuplesPerHash, int *refPosHashMap, char *qrs, uchar *qryLen,
		int *qryDistance, int maxQryLen, char *refIdx, int *refPos,
		int maxNumHits, int seedLen, int randNum, int idealArrSize);

void lookupTable13CpyConstMemToGPU();

__device__ int getHash13_gpu(char *str, int len);

__device__ void getUnhash13_gpu(int hash, char *s, int seedLen);

__device__ void getRefIdxAndPos13_gpu(int compositeVal, int *refIdx, int *refPos,
		int seedLen);

__device__ int arrSearch13_gpu(int *arr, int arrSize, int num);

__device__ int getRandNum13_gpu();

__device__ void arrGetRandomNums13_gpu(int n, int lLimit, int uLimit, int *arr);

#endif /* LOOKUPTABLE13_H_ */
