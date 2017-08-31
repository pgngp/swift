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

#include <check.h>
#include <stdio.h>
#include "../src/common.h"
#include "testLookupTable6.h"
#include "../src/lookupTable6.h"
#include "../src/refMap.h"
#include "../src/common.h"
#include "../src/preprocess.h"
#include <limits.h>


/**
 * Tests @a sort_gpu6 function.
 */
START_TEST(sort_gpu)
{
	/* Array size is 2. */
	{
		long refPos[] = {10, 5};
		int arrSize = 2;

		long *refPos_d;
		cudaMalloc(&refPos_d, arrSize * sizeof(long));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos_d, refPos, arrSize * sizeof(long),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Using 1 thread. */
		sort_gpu_wrap6<<<1, 1>>>(refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(long),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refPos[0] != 5 || refPos[1] != 10)
			fail("Incorrect behavior when the array size is 2 and number of "
					"threads launched is 1.\n");

		/* Using 2 threads. */
		sort_gpu_wrap6<<<1, 2>>>(refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(long),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refPos[0] != 5 || refPos[1] != 10)
			fail("Incorrect behavior when the array size is 2 and number of "
					"threads launched is 2.\n");

		/* Using 3 threads. */
		sort_gpu_wrap6<<<1, 3>>>(refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(long),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refPos[0] != 5 || refPos[1] != 10)
			fail("Incorrect behavior when the array size is 2 and number of "
					"threads launched is 3.\n");

		cudaFree(refPos_d);
	}

	/* Array size is 5. */
	{
		long refPos[] = {18, 13, 10, 5, 0};
		int arrSize = 5;

		long *refPos_d;
		cudaMalloc(&refPos_d, arrSize * sizeof(long));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos_d, refPos, arrSize * sizeof(long),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Using 1 thread. */
		sort_gpu_wrap6<<<1, 1>>>(refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(long),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refPos[0] != 0 || refPos[1] != 5 || refPos[2] != 10
				|| refPos[3] != 13 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 1.\n");

		/* Using 2 threads. */
		sort_gpu_wrap6<<<1, 2>>>(refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(long),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refPos[0] != 0 || refPos[1] != 5 || refPos[2] != 10
				|| refPos[3] != 13 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 2.\n");

		/* Using 3 threads. */
		sort_gpu_wrap6<<<1, 3>>>(refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(long),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refPos[0] != 0 || refPos[1] != 5 || refPos[2] != 10
				|| refPos[3] != 13 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 3.\n");

		/* Using 4 threads. */
		sort_gpu_wrap6<<<1, 4>>>(refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(long),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refPos[0] != 0 || refPos[1] != 5 || refPos[2] != 10
				|| refPos[3] != 13 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 4.\n");

		/* Using 5 threads. */
		sort_gpu_wrap6<<<1, 5>>>(refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(long),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refPos[0] != 0 || refPos[1] != 5 || refPos[2] != 10
				|| refPos[3] != 13 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 5.\n");

		/* Using 6 threads. */
		sort_gpu_wrap6<<<1, 6>>>(refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(long),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refPos[0] != 0 || refPos[1] != 5 || refPos[2] != 10
				|| refPos[3] != 13 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 6.\n");

		cudaFree(refPos_d);
	}
}
END_TEST


/**
 * Tests the @a getHash_gpu6 function.
 */
START_TEST(getHash_gpu)
{
	char str[] = "ACGT";
	int len = strlen(str);

	char *str_d;
	cudaMalloc((void **) &str_d, len * sizeof(char));
	PRINT_CUDA_ERROR()
	cudaMemcpy(str_d, str, len * sizeof(char), cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	int *hash_d;
	cudaMalloc((void **) &hash_d, sizeof(int));
	PRINT_CUDA_ERROR()

	lookupTable6CpyConstMemToGPU();

	getHash_gpu_wrap6<<<1, 1>>>(str_d, len, hash_d);
	PRINT_CUDA_ERROR()

	int hash;
	cudaMemcpy(&hash, hash_d, sizeof(int), cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	if (hash != 228)
		fail("Incorrect behavior.\n");

	cudaFree(str_d);
	cudaFree(hash_d);
}
END_TEST


/**
 * Tests @a createClusters_gpu function.
 */
START_TEST(createClusters_gpu)
{
	/* Array size is 5, there are 4 clusters, and 1 biggest cluster. */
	{
		long refPos[] = {53418655744, 72057649067196416, 216172834458697728,
				216172834458697728, 648518401907490816};
		int arrSize = 5;

		long *refPos_d;
		cudaMalloc(&refPos_d, arrSize * sizeof(long));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos_d, refPos, arrSize * sizeof(long),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *clusterSize_d;
		cudaMalloc(&clusterSize_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		short *biggestClusterSize_d;
		cudaMalloc(&biggestClusterSize_d, sizeof(short));
		PRINT_CUDA_ERROR()

		short *numClusters_d;
		cudaMalloc(&numClusters_d, sizeof(short));
		PRINT_CUDA_ERROR()

		createClusters_gpu_wrap<<<1, 1>>>(refPos_d, clusterSize_d,
				arrSize, biggestClusterSize_d, numClusters_d);

		int clusterSize[arrSize];
		cudaMemcpy(clusterSize, clusterSize_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short biggestClusterSize;
		cudaMemcpy(&biggestClusterSize, biggestClusterSize_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short numClusters;
		cudaMemcpy(&numClusters, numClusters_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (numClusters != 4 || biggestClusterSize != 2 || clusterSize[0] != 1
				|| clusterSize[1] != 1 || clusterSize[2] != 2
				|| clusterSize[3] != 1)
			fail("Incorrect behavior when array size is 5, there are 4 clusters, "
					" and 1 biggest cluster.\n");

		cudaFree(refPos_d);
		cudaFree(clusterSize_d);
		cudaFree(biggestClusterSize_d);
		cudaFree(numClusters_d);
	}

	/* Array size is 5, there is 1 cluster and 1 biggest cluster. */
	{
		long refPos[] = {72057646382841856, 72057646382841856, 72057646382841856,
				72057646382841856, 72057646382841856};
		int arrSize = 5;

		long *refPos_d;
		cudaMalloc(&refPos_d, arrSize * sizeof(long));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos_d, refPos, arrSize * sizeof(long),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *clusterSize_d;
		cudaMalloc(&clusterSize_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		short *biggestClusterSize_d;
		cudaMalloc(&biggestClusterSize_d, sizeof(short));
		PRINT_CUDA_ERROR()

		short *numClusters_d;
		cudaMalloc(&numClusters_d, sizeof(short));
		PRINT_CUDA_ERROR()

		createClusters_gpu_wrap<<<1, 1>>>(refPos_d, clusterSize_d,
				arrSize, biggestClusterSize_d, numClusters_d);

		int clusterSize[arrSize];
		cudaMemcpy(clusterSize, clusterSize_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short biggestClusterSize;
		cudaMemcpy(&biggestClusterSize, biggestClusterSize_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short numClusters;
		cudaMemcpy(&numClusters, numClusters_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (numClusters != 1 || biggestClusterSize != 5 || clusterSize[0] != 5)
			fail("Incorrect behavior when array size is 5, there is 1 cluster "
					"and 1 biggest cluster.\n");

		cudaFree(refPos_d);
		cudaFree(clusterSize_d);
		cudaFree(biggestClusterSize_d);
		cudaFree(numClusters_d);
	}
}
END_TEST


/**
 * Tests @a findBiggestClusters_gpu function.
 */
START_TEST(findBiggestClusters_gpu)
{
	/* Array size is 5, there are 4 clusters, and 1 biggest cluster. */
	{
		long refPos[] = {53418655744, 72057649067196416, 216172834458697728,
				216172834458697728, 648518401907490816};
		int arrSize = 5;

		long *refPos_d;
		cudaMalloc(&refPos_d, arrSize * sizeof(long));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos_d, refPos, arrSize * sizeof(long),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *clusterSize_d;
		cudaMalloc(&clusterSize_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		short *biggestClusterSize_d;
		cudaMalloc(&biggestClusterSize_d, sizeof(short));
		PRINT_CUDA_ERROR()

		short *numClusters_d;
		cudaMalloc(&numClusters_d, sizeof(short));
		PRINT_CUDA_ERROR()

		createClusters_gpu_wrap<<<1, 1>>>(refPos_d, clusterSize_d,
				arrSize, biggestClusterSize_d, numClusters_d);

		int clusterSize[arrSize];
		cudaMemcpy(clusterSize, clusterSize_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short biggestClusterSize;
		cudaMemcpy(&biggestClusterSize, biggestClusterSize_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short numClusters;
		cudaMemcpy(&numClusters, numClusters_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short *biggestClusters_d;
		cudaMalloc(&biggestClusters_d, arrSize * sizeof(short));
		PRINT_CUDA_ERROR()

		char *numBiggestClusters_d;
		cudaMalloc(&numBiggestClusters_d, sizeof(char));
		PRINT_CUDA_ERROR()

		findBiggestClusters_gpu_wrap<<<1, 1>>>(numClusters, biggestClusterSize,
				clusterSize_d, biggestClusters_d, numBiggestClusters_d);
		PRINT_CUDA_ERROR()

		short biggestClusters[arrSize];
		cudaMemcpy(biggestClusters, biggestClusters_d, arrSize * sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		char numBiggestClusters;
		cudaMemcpy(&numBiggestClusters, numBiggestClusters_d,
				sizeof(char), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (numBiggestClusters != 1 || biggestClusters[0] != 2)
			fail("Incorrect behavior when array size is 5, there are 4 clusters, "
					" and 1 biggest cluster.\n");

		cudaFree(refPos_d);
		cudaFree(clusterSize_d);
		cudaFree(biggestClusterSize_d);
		cudaFree(numClusters_d);
		cudaFree(biggestClusters_d);
		cudaFree(numBiggestClusters_d);
	}

	/* Array size is 5, there are 3 clusters, and 2 biggest clusters. */
	{
		long refPos[] = {1342177280, 1342177280, 216172834458697728,
				216172834458697728, 648518401907490816};
		int arrSize = 5;

		long *refPos_d;
		cudaMalloc(&refPos_d, arrSize * sizeof(long));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos_d, refPos, arrSize * sizeof(long),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *clusterSize_d;
		cudaMalloc(&clusterSize_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		short *biggestClusterSize_d;
		cudaMalloc(&biggestClusterSize_d, sizeof(short));
		PRINT_CUDA_ERROR()

		short *numClusters_d;
		cudaMalloc(&numClusters_d, sizeof(short));
		PRINT_CUDA_ERROR()

		createClusters_gpu_wrap<<<1, 1>>>(refPos_d, clusterSize_d, arrSize,
				biggestClusterSize_d, numClusters_d);

		int clusterSize[arrSize];
		cudaMemcpy(clusterSize, clusterSize_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short biggestClusterSize;
		cudaMemcpy(&biggestClusterSize, biggestClusterSize_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short numClusters;
		cudaMemcpy(&numClusters, numClusters_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short *biggestClusters_d;
		cudaMalloc(&biggestClusters_d, arrSize * sizeof(short));
		PRINT_CUDA_ERROR()

		char *numBiggestClusters_d;
		cudaMalloc(&numBiggestClusters_d, sizeof(char));
		PRINT_CUDA_ERROR()

		findBiggestClusters_gpu_wrap<<<1, 1>>>(numClusters, biggestClusterSize,
				clusterSize_d, biggestClusters_d, numBiggestClusters_d);
		PRINT_CUDA_ERROR()

		short biggestClusters[arrSize];
		cudaMemcpy(biggestClusters, biggestClusters_d, arrSize * sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		char numBiggestClusters;
		cudaMemcpy(&numBiggestClusters, numBiggestClusters_d,
				sizeof(char), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (numBiggestClusters != 2 || biggestClusters[0] != 0
				|| biggestClusters[1] != 2)
			fail("Incorrect behavior when array size is 5, there are 3 clusters, "
					" and 2 biggest clusters.\n");

		cudaFree(refPos_d);
		cudaFree(clusterSize_d);
		cudaFree(biggestClusterSize_d);
		cudaFree(numClusters_d);
		cudaFree(biggestClusters_d);
		cudaFree(numBiggestClusters_d);
	}
}
END_TEST


/**
 * Tests @a assignResults_gpu function.
 */
START_TEST(assignResults_gpu)
{
	/* Num hits is 1 and max allowed hits is 1. */
	{
		char maxHits = 1;
		char numBgstHits = 1;
		long refPos[] = {1, 2, 5, 72057595380105219, 72057595380105219};
		int arrSize = 5;

		char *refIdx2_d;
		cudaMalloc(&refIdx2_d, maxHits * sizeof(char));
		PRINT_CUDA_ERROR()

		long *refPos1_d;
		cudaMalloc(&refPos1_d, arrSize * sizeof(long));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos1_d, refPos, arrSize * sizeof(long),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *refPos2_d;
		cudaMalloc(&refPos2_d, sizeof(int));
		PRINT_CUDA_ERROR()

		short *bgstClust_d;
		cudaMalloc(&bgstClust_d, sizeof(short));
		PRINT_CUDA_ERROR()

		short bgstClust = 3;
		cudaMemcpy(bgstClust_d, &bgstClust, numBgstHits * sizeof(short),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		srand(time(NULL));
		int randNum = rand();
		assignResults_gpu_wrap<<<1, 1>>>(maxHits, numBgstHits, refPos1_d,
				refIdx2_d, refPos2_d, bgstClust_d, randNum);
		PRINT_CUDA_ERROR()

		char refIdx2;
		cudaMemcpy(&refIdx2, refIdx2_d, numBgstHits * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int refPos2;
		cudaMemcpy(&refPos2, refPos2_d, numBgstHits * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx2 != 1 || refPos2 != 3)
			fail("Incorrect behavior when number of hits is 1 and max allowed "
					"hits is 1.\n");

		cudaFree(refIdx2_d);
		cudaFree(refPos1_d);
		cudaFree(refPos2_d);
		cudaFree(bgstClust_d);
	}

	/* Num hits is 2 and max allowed hits is 1. */
	{
		char maxHits = 1;
		char numBgstHits = 2;
		long refPos[] = {2, 2, 5, 72057595380105219, 72057595380105219};
		int arrSize = 5;

		char *refIdx2_d;
		cudaMalloc(&refIdx2_d, maxHits * sizeof(char));
		PRINT_CUDA_ERROR()

		long *refPos1_d;
		cudaMalloc(&refPos1_d, arrSize * sizeof(long));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos1_d, refPos, arrSize * sizeof(long),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *refPos2_d;
		cudaMalloc(&refPos2_d, maxHits * sizeof(int));
		PRINT_CUDA_ERROR()

		short *bgstClust_d;
		cudaMalloc(&bgstClust_d, arrSize * sizeof(short));
		PRINT_CUDA_ERROR()

		short bgstClust[] = {0, 3};
		cudaMemcpy(bgstClust_d, &bgstClust, arrSize * sizeof(short),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		srand(time(NULL));
		int randNum = rand();
		assignResults_gpu_wrap<<<1, 1>>>(maxHits, numBgstHits, refPos1_d,
				refIdx2_d, refPos2_d, bgstClust_d, randNum);
		PRINT_CUDA_ERROR()

		char refIdx2;
		cudaMemcpy(&refIdx2, refIdx2_d, maxHits * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int refPos2;
		cudaMemcpy(&refPos2, refPos2_d, maxHits * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if ((refIdx2 != 0 && refIdx2 != 1) || (refPos2 != 2 && refPos2 != 3))
			fail("Incorrect behavior when number of hits is 2 and max allowed "
					"hits is 1.\n");

		cudaFree(refIdx2_d);
		cudaFree(refPos1_d);
		cudaFree(refPos2_d);
		cudaFree(bgstClust_d);
	}

	/* Num hits is 1 and max allowed hits is 2. */
	{
		char maxHits = 2;
		char numBgstHits = 1;

		long refPos[] = {1, 2, 5, 72057595380105219, 72057595380105219};
		int arrSize = 5;

		char *refIdx2_d;
		cudaMalloc(&refIdx2_d, maxHits * sizeof(char));
		PRINT_CUDA_ERROR()

		long *refPos1_d;
		cudaMalloc(&refPos1_d, arrSize * sizeof(long));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos1_d, refPos, arrSize * sizeof(long),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *refPos2_d;
		cudaMalloc(&refPos2_d, maxHits * sizeof(int));
		PRINT_CUDA_ERROR()

		short *bgstClust_d;
		cudaMalloc(&bgstClust_d, arrSize * sizeof(short));
		PRINT_CUDA_ERROR()

		short bgstClust[] = {3};
		cudaMemcpy(bgstClust_d, &bgstClust, arrSize * sizeof(short),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		srand(time(NULL));
		int randNum = rand();
		assignResults_gpu_wrap<<<1, 1>>>(maxHits, numBgstHits, refPos1_d,
				refIdx2_d, refPos2_d, bgstClust_d, randNum);
		PRINT_CUDA_ERROR()

		char refIdx2;
		cudaMemcpy(&refIdx2, refIdx2_d, maxHits * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int refPos2;
		cudaMemcpy(&refPos2, refPos2_d, maxHits * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx2 != 1 || refPos2 != 3)
			fail("Incorrect behavior when number of hits is 1 and max allowed "
					"hits is 2.\n");

		cudaFree(refIdx2_d);
		cudaFree(refPos1_d);
		cudaFree(refPos2_d);
		cudaFree(bgstClust_d);
	}
}
END_TEST


/**
 * Tests @a cpyHitsFromGlobalToShr_gpu function.
 */
START_TEST(cpyHitsFromGlobalToShr_gpu)
{
	int seedLen = 5;
	int numDistinctTuples = (int) pow((float) DNA_ALPHABET_SIZE,
			(int) seedLen);
	int *keys = (int *) calloc(numDistinctTuples, sizeof(int));
	keys[228] = 1;
	keys[377] = 3;
	keys[916] = 4;

	int numTotalTuples = 4;
	int *vals = (int *) calloc((numTotalTuples + 1), sizeof(int));
	vals[0] = -1;
	vals[1] = 1;
	vals[2] = 33554435;
	vals[3] = 2;
	vals[4] = 0;

	int *numRptsPerTuple = (int *) calloc(numDistinctTuples, sizeof(int));
	numRptsPerTuple[228] = 2;
	numRptsPerTuple[377] = 1;
	numRptsPerTuple[916] = 1;

	int *keys_d;
	cudaMalloc(&keys_d, numDistinctTuples * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(keys_d, keys, numDistinctTuples * sizeof(int),
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	int *vals_d;
	cudaMalloc(&vals_d, (numTotalTuples + 1) * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(vals_d, vals, (numTotalTuples + 1) * sizeof(int),
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	int *numRptsPerTuple_d;
	cudaMalloc(&numRptsPerTuple_d, numDistinctTuples * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(numRptsPerTuple_d, numRptsPerTuple,
			numDistinctTuples * sizeof(int), cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	int arrSize = 100;
	long *refPos_d;
	cudaMalloc(&refPos_d, arrSize * sizeof(long));
	PRINT_CUDA_ERROR()

	int numQryTuples = 6;
	int *hashes = (int *) calloc(numQryTuples, sizeof(int));
	hashes[0] = 228;
	hashes[1] = 313;
	hashes[2] = 590;
	hashes[3] = 915;
	hashes[4] = 484;
	hashes[5] = 377;

	int *hashes_d;
	cudaMalloc(&hashes_d, numQryTuples * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(hashes_d, hashes, numQryTuples * sizeof(int),
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	short maxQryLen = 200;
	cpyHitsFromGlobalToShr_gpu_wrap<<<1, 1>>>(refPos_d, hashes_d, arrSize,
			keys_d, vals_d, numRptsPerTuple_d, numQryTuples, seedLen, maxQryLen);
	PRINT_CUDA_ERROR()

	long *refPos = (long *) calloc(arrSize, sizeof(long));
	cudaMemcpy(refPos, refPos_d, arrSize * sizeof(long), cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	if (refPos[0] != 55029268485 || refPos[1] != 72057651751550991
			|| refPos[2] != 55029268490)
		fail("Incorrect behavior.\n");

	cudaFree(keys_d);
	cudaFree(vals_d);
	cudaFree(numRptsPerTuple_d);
	cudaFree(refPos_d);
	cudaFree(hashes_d);
	free(keys);
	free(vals);
	free(numRptsPerTuple);
	free(refPos);
	free(hashes);
}
END_TEST


/**
 * Tests @a lookupTable5MapQry_gpu function.
 */
START_TEST(lookupTable6MapQry_gpu)
{
	/* Case 1. */
	{
		int seedLen = 5;
		int numDistinctTuples = (int) pow((float) DNA_ALPHABET_SIZE,
				(int) seedLen);
		int *keys = (int *) calloc(numDistinctTuples, sizeof(int));
		keys[228] = 1;
		keys[377] = 3;
		keys[916] = 4;

		int numTotalTuples = 4;
		int *vals = (int *) calloc((numTotalTuples + 1), sizeof(int));
		vals[0] = -1;
		vals[1] = 1;
		vals[2] = 3;
		vals[3] = 2;
		vals[4] = 0;

		int *numRptsPerTuple = (int *) calloc(numDistinctTuples, sizeof(int));
		numRptsPerTuple[228] = 2;
		numRptsPerTuple[377] = 1;
		numRptsPerTuple[916] = 1;

		int *keys_d;
		cudaMalloc(&keys_d, numDistinctTuples * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(keys_d, keys, numDistinctTuples * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *vals_d;
		cudaMalloc(&vals_d, (numTotalTuples + 1) * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(vals_d, vals, (numTotalTuples + 1) * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *numRptsPerTuple_d;
		cudaMalloc(&numRptsPerTuple_d, numDistinctTuples * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(numRptsPerTuple_d, numRptsPerTuple,
				numDistinctTuples * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int maxQrySeqLen = 12;
		int numQrs = 2;
		char *qrs = (char *) calloc(maxQrySeqLen * numQrs, sizeof(char));
		strcpy(qrs, "ACGTACGTCC");
		strcpy(qrs + maxQrySeqLen, "GGACGTACGT");

		uchar *qryLen = (uchar *) calloc(numQrs, sizeof(uchar));
		qryLen[0] = 10;
		qryLen[1] = 10;

		char *qrs_d;
		cudaMalloc(&qrs_d, maxQrySeqLen * numQrs * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(qrs_d, qrs, maxQrySeqLen * numQrs * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		uchar *qryLen_d;
		cudaMalloc(&qryLen_d, numQrs * sizeof(uchar));
		PRINT_CUDA_ERROR()
		cudaMemcpy(qryLen_d, qryLen, numQrs * sizeof(uchar),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		char *refIdx_d;
		cudaMalloc(&refIdx_d, numQrs * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemset(refIdx_d, CHAR_MAX, numQrs * sizeof(char));
		PRINT_CUDA_ERROR()

		int *refPos_d;
		cudaMalloc(&refPos_d, numQrs * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemset(refPos_d, -1, numQrs * sizeof(int));
		PRINT_CUDA_ERROR()

		lookupTable6CpyConstMemToGPU();

		int maxHits = 1;
		srand(time(NULL));
		int randNum = rand();
		int arrSize = 100;
		int shrMem = (4 * arrSize * sizeof(int)) + (maxHits * sizeof(int));
		lookupTable6MapQry2_gpu<<<2, 1, shrMem>>>(keys_d, vals_d,
				numRptsPerTuple_d, qrs_d, qryLen_d, maxQrySeqLen, refIdx_d,
				refPos_d, maxHits, seedLen, randNum, arrSize);
		PRINT_CUDA_ERROR()

		char *refIdx = (char *) malloc(numQrs * sizeof(char));
		cudaMemcpy(refIdx, refIdx_d, numQrs * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *refPos = (int *) malloc(numQrs * sizeof(int));
		cudaMemcpy(refPos, refPos_d, numQrs * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx[0] != 0 || refPos[0] != 5 || refIdx[1] != 0
				|| (refPos[1] != 5 && refPos[1] != 15))
			fail("Incorrect behavior.\n");

		cudaFree(keys_d);
		cudaFree(vals_d);
		cudaFree(numRptsPerTuple_d);
		cudaFree(qrs_d);
		cudaFree(qryLen_d);
		cudaFree(refIdx_d);
		cudaFree(refPos_d);
		free(keys);
		free(vals);
		free(numRptsPerTuple);
		free(qrs);
		free(qryLen);
		free(refIdx);
		free(refPos);
	}
}
END_TEST


/**
 * Tests @a initializeShrMem_gpu function.
 */
START_TEST(initializeShrMem_gpu)
{
#define	REF_POS_MASK3	9151314442816847872 /* Binary: 0111 1111 0000 0000 0000
0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 */
	int arrSize = 5;
	long *refPos = (long *) calloc(arrSize, sizeof(long));
	long *refPos_d;
	cudaMalloc(&refPos_d, arrSize * sizeof(long));
	PRINT_CUDA_ERROR()

	int *clusterSize = (int *) calloc(arrSize, sizeof(int));
	int *clusterSize_d;
	cudaMalloc(&clusterSize_d, arrSize * sizeof(int));
	PRINT_CUDA_ERROR()

	intializeShrMem_gpu_wrap<<<1, 5>>>(refPos_d, clusterSize_d, arrSize);

	cudaMemcpy(refPos, refPos_d, arrSize * sizeof(long),
			cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	cudaMemcpy(clusterSize, clusterSize_d, arrSize * sizeof(int),
			cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	int i;
	for (i = 0; i < arrSize; ++i)
	{
		if (refPos[i] != REF_POS_MASK3)
			fail("Incorrect behavior: Reference position array is not "
					"initialized properly.\n");
		if (clusterSize[i] != -1)
			fail("Incorrect behavior: Cluster size array is not "
					"initialized properly.\n");
	}

	cudaFree(refPos_d);
	cudaFree(clusterSize_d);
	free(refPos);
	free(clusterSize);
}
END_TEST


/**
 * Creates test suite.
 */
Suite *lookupTable6Suite(void)
{
	Suite *s = suite_create("lookupTable6");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, getHash_gpu);
	tcase_add_test(testCaseCore, sort_gpu);
	tcase_add_test(testCaseCore, createClusters_gpu);
	tcase_add_test(testCaseCore, findBiggestClusters_gpu);
	tcase_add_test(testCaseCore, assignResults_gpu);
	tcase_add_test(testCaseCore, initializeShrMem_gpu);
	tcase_add_test(testCaseCore, cpyHitsFromGlobalToShr_gpu);
//	tcase_add_test(testCaseCore, lookupTable6MapQry_gpu);
	suite_add_tcase (s, testCaseCore);

	return s;
}
