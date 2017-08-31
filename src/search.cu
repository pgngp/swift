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

#include "search.h"
#include "common.h"
#include "input.h"
#include "preprocess.h"
#include <string.h>
#include "reference.h"
#include "refPosMap.h"
#include "refNameMap.h"
#include <limits.h>
#include <time.h>
#include "lookupTable.h"
#include "lookupTable2.h"
#include "lookupTable3.h"
#include "lookupTable4.h"
#include "lookupTable5.h"
#include "lookupTable6.h"
#include "lookupTable7.h"
#include "lookupTable8.h"
#include "lookupTable9.h"
#include "lookupTable10.h"
#include "lookupTable11.h"
#include <pthread.h>

#define NUM_RESULTS		8
#define	FORWARD_STRAND	0
#define REV_STRAND		1
#define	MIN_QRS_FOR_GPU	1000
#define	SHR_MEM_MARGIN	256


pthread_mutex_t mutex2;


typedef struct data
{
	Query *qry;
	int qryLen;
	char *refIdx1;
	int *shift1;
	int *refPos1;
	char *refIdx2;
	int *shift2;
	int *refPos2;
	int numHits;
	int numHits1;
	int numHits2;
} Data;


typedef struct data2
{
	char *qrySeq;
	int qryLen;
	char *refIdx;
	int *shift;
	int *refPos;
} Data2;



/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param qryFileName Query file name.
 * @param numQrs Number of queries in the query file.
 * @param qryLength Length of each query in the query file.
 * @param refFileName Reference file name.
 * @param seedLen Seed length.
 * @param matchFile File in which filtered reference sequence fragment
 * positions are stored.
 * @param maxNumHits Max number of hits per query.
 * @param numMatches Number of matches obtained.
 */
void searchQueries(char *qryFileName, int numQrs, int qryLen,
		char *refFileName, int seedLen, char *matchFile, uint maxNumHits,
		int *numMatches)
{
	if (qryFileName == NULL || qryLen == 0 || refFileName == NULL
			|| seedLen <= 0 || matchFile == NULL || numQrs <= 0)
		return;

	/* Step 1: Preprocess reference */
	fprintf(stderr, "Preprocessing reference sequences...\n");
	lookupTable3Create(refFileName, seedLen, maxNumHits, TUPLE_IGNORE_THRES);
	fprintf(stderr, "done.\n");

	/* Step 2: Preprocess query */
	fprintf(stderr, "Preprocessing query sequences...");
	Query *qry = qryCreate2();
	qryInitialize(qryFileName, qryLen);
	fprintf(stderr, "done.\n");

	/* Step 3: Find matching reference positions for each query and its
	 * reverse complement */
	fprintf(stderr, "Mapping queries to references...\n");
	FILE *filePtr = fopen(matchFile, "w");
	HitList *bestMatches = NULL;
	int isRevComp = 0, i;
	time_t startTime, endTime;
	double diffTime = 0.0;
	*numMatches = 0;
	for (i = 0; i < numQrs; ++i)
	{
		time(&startTime);

		qry = qryGetNext();

		/* Find matching reference positions for query */
		isRevComp = 0;
		bestMatches = lookupTable3MapQry(qry, qryLen, isRevComp);
		printHeuristicMatches(qry, isRevComp, bestMatches, filePtr, seedLen,
				numMatches);
		hitListDelete(bestMatches);

		/* Find matching reference position for reverse complement of query */
		isRevComp = 1;
		bestMatches = lookupTable3MapQry(qry, qryLen, isRevComp);
		printHeuristicMatches(qry, isRevComp, bestMatches, filePtr, seedLen,
				numMatches);
		hitListDelete(bestMatches);

		time(&endTime);
		diffTime += difftime(endTime, startTime);
		if ((i + 1) % 100000 == 0)
		{
			fprintf(stderr, "   Searched %d queries...\n", i + 1);
			fprintf(stderr, "      Time = %.2lf secs\n", diffTime);
			diffTime = 0.0;
//			break;
		}
	}

	/* Clean up resources */
	fclose(filePtr);
	qryClean();
//	qryDelete(qry);

	/* Free memory occupied by lookup table. */
	fprintf(stderr, "Deleting lookup table...");
	lookupTable3Delete();
	fprintf(stderr, "done.\n");
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param qryFileName Query file name.
 * @param numQrs Number of queries in the query file.
 * @param qryLength Length of each query in the query file.
 * @param refFileName Reference file name.
 * @param seedLen Seed length.
 * @param matchFile File in which filtered reference sequence fragment
 * positions are stored.
 * @param maxNumHits Max number of hits per query.
 * @param numMatches Number of matches obtained.
 */
void searchQueries2(char *qryFileName, int numQrs, int qryLen,
		char *refFileName, int seedLen, char *matchFile, uint maxNumHits,
		int *numMatches)
{
	if (qryFileName == NULL || qryLen == 0 || refFileName == NULL
			|| seedLen <= 0 || matchFile == NULL || numQrs <= 0)
		return;

	/* Step 1: Preprocess reference */
	fprintf(stderr, "Preprocessing reference sequences...\n");
	lookupTable4Create(refFileName, seedLen, maxNumHits, TUPLE_IGNORE_THRES);
	fprintf(stderr, "done.\n");

	/* Step 2: Preprocess query */
	fprintf(stderr, "Preprocessing query sequences...");
	Query *qry = qryCreate2();
	qryInitialize(qryFileName, qryLen);
	fprintf(stderr, "done.\n");

	/* Step 3: Find matching reference positions for each query and its
	 * reverse complement */
	fprintf(stderr, "Mapping queries to references...\n");
	FILE *filePtr = fopen(matchFile, "w");
	int isRevComp = 0, i, numHits;
	time_t startTime, endTime;
	double diffTime = 0.0;
	char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
	int *shift = (int *) calloc(maxNumHits, sizeof(int));
	int *refPos = (int *) calloc(maxNumHits, sizeof(int));
	*numMatches = 0;
	numHits = 0;
	for (i = 0; i < numQrs; ++i)
	{
		time(&startTime);

		qry = qryGetNext();

		/* Find matching reference positions for query */
		isRevComp = 0;
		numHits = lookupTable4MapQry(qry, qryLen, isRevComp, refIdx, shift,
				refPos);
		printHeuristicMatches4(qry, isRevComp, refIdx, shift, refPos, numHits,
				filePtr);
		*numMatches += numHits;

		/* Find matching reference position for reverse complement of query */
		isRevComp = 1;
		numHits = lookupTable4MapQry(qry, qryLen, isRevComp, refIdx, shift,
				refPos);
		printHeuristicMatches4(qry, isRevComp, refIdx, shift, refPos, numHits,
				filePtr);
		*numMatches += numHits;

		time(&endTime);
		diffTime += difftime(endTime, startTime);
		if ((i + 1) % 100000 == 0)
		{
			fprintf(stderr, "   Searched %d queries...\n", i + 1);
			fprintf(stderr, "      Time = %.2lf secs\n", diffTime);
			diffTime = 0.0;
			break;
		}
	}

	/* Clean up resources */
	fclose(filePtr);
	qryClean();
//	qryDelete(qry);

	/* Free memory occupied by lookup table. */
	fprintf(stderr, "Deleting lookup table...");
	lookupTable4Delete();
	free(refIdx);
	free(shift);
	free(refPos);
	fprintf(stderr, "done.\n");
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param qryFileName Query file name.
 * @param numQrs Number of queries in the query file.
 * @param qryLength Length of each query in the query file.
 * @param refFileName Reference file name.
 * @param seedLen Seed length.
 * @param matchFile File in which filtered reference sequence fragment
 * positions are stored.
 * @param maxNumHits Max number of hits per query.
 * @param numMatches Number of matches obtained.
 * @param numThreads Number of CPU threads to use.
 * @param tupleIgnoreThres Tuple ignore threshold.
 */
void searchQueries3(char *qryFileName, int numQrs, int qryLen,
		char *refFileName, int seedLen, char *matchFile, uint maxNumHits,
		int *numMatches, int numThreads, int tupleIgnoreThres)
{
	if (qryFileName == NULL || refFileName == NULL || seedLen <= 0
			|| matchFile == NULL || numQrs <= 0)
		return;

	/* Step 1: Preprocess reference */
	fprintf(stderr, "Preprocessing reference sequences...\n");
	int totalTuples;
	lookupTable5Create(refFileName, seedLen, maxNumHits, tupleIgnoreThres,
			&totalTuples);
//	lookupTable5Create("/home/pgupta/data/gpusw/keys1000.txt",
//			"/home/pgupta/data/gpusw/vals1000.txt", maxNumHits, seedLen,
//			tupleIgnoreThres);
	fprintf(stderr, "done.\n");

	/* Step 2: Preprocess query */
	fprintf(stderr, "Preprocessing query sequences...");
//	qryInitialize(qryFileName, qryLen);
	qryInitialize(qryFileName, MAX_QRY_SEQ_LENGTH);
	fprintf(stderr, "done.\n");

	/* Step 3: Find matching reference positions for each query and its
	 * reverse complement */
	fprintf(stderr, "Mapping queries to references...\n");
	FILE *filePtr = fopen(matchFile, "w");
	int i;
	time_t startTime, endTime;
	double diffTime = 0.0;

	Data **data = (Data **) malloc(numThreads * sizeof(Data *));
	for (i = 0; i < numThreads; ++i)
	{
		data[i] = (Data *) malloc(sizeof(Data));
		data[i]->refIdx1 = (char *) calloc(maxNumHits, sizeof(char));
		data[i]->shift1 = (int *) calloc(maxNumHits, sizeof(int));
		data[i]->refPos1 = (int *) calloc(maxNumHits, sizeof(int));
		data[i]->refIdx2 = (char *) calloc(maxNumHits, sizeof(char));
		data[i]->shift2 = (int *) calloc(maxNumHits, sizeof(int));
		data[i]->refPos2 = (int *) calloc(maxNumHits, sizeof(int));
		data[i]->qry = qryCreate2();
	}

	FILE *qryFilePtr = fopen(qryFileName, "r");

	*numMatches = 0;
	pthread_t thread[numThreads];
	int j;
	for (i = 0; i < numQrs; i += numThreads)
	{
		time(&startTime);

		/* Create threads. */
		for (j = 0; j < numThreads; ++j)
		{
			qryGetNext2(qryFilePtr, data[j]->qry);
			pthread_create(&(thread[j]), NULL, getHits, data[j]);
		}

		/* Wait for threads to complete, then print the intermediate results. */
		for (j = 0; j < numThreads; ++j)
		{
			pthread_join(thread[j], NULL);
			printHeuristicMatches4(data[j]->qry, FORWARD_STRAND, data[j]->refIdx1,
					data[j]->shift1, data[j]->refPos1, data[j]->numHits1, filePtr);
			printHeuristicMatches4(data[j]->qry, REV_STRAND, data[j]->refIdx2,
					data[j]->shift2, data[j]->refPos2, data[j]->numHits2, filePtr);
			*numMatches += data[j]->numHits;
		}

		/* Show progress. */
		time(&endTime);
		diffTime += difftime(endTime, startTime);
		if ((i + numThreads) % 25000 == 0)
		{
			fprintf(stderr, "   Searched %d queries...", (i + numThreads));
			fprintf(stderr, "Time = %.2lf secs\n", diffTime);
			diffTime = 0.0;
//			break;
		}
	}

	/* Clean up resources */
	fclose(filePtr);
	fclose(qryFilePtr);
	qryClean();
//	qryDelete(qry);

	/* Free memory occupied by lookup table. */
	fprintf(stderr, "Deleting lookup table...");
	lookupTable5Delete();
	for (i = 0; i < numThreads; ++i)
	{
		free(data[i]->qry);
		free(data[i]->refIdx1);
		free(data[i]->shift1);
		free(data[i]->refPos1);
		free(data[i]->refIdx2);
		free(data[i]->shift2);
		free(data[i]->refPos2);
		free(data[i]);
	}
	free(data);

	fprintf(stderr, "done.\n");
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param	qryFileName	Query file name.
 * @param	refFileName	Reference file name.
 * @param	matchFile	File in which filtered reference sequence fragment
 * positions are stored.
 * @param	seedLen		Seed length.
 * @param	maxNumHits	Max number of hits per query.
 * @param 	numMatches	Total number of matches obtained.
 */
void searchQueries4(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches)
{
	if (qryFileName == NULL || refFileName == NULL || matchFile == NULL
			|| seedLen <= 0)
		return;

	cudaSetDevice(SELECTED_DEVICE_INDEX);
	PRINT_CUDA_ERROR()

	/* Preprocess reference */
	fprintf(stderr, "Preprocessing reference sequences...\n");
	int numTotalTuples;
//	lookupTable5Create2(refFileName, seedLen, maxNumHits, TUPLE_IGNORE_THRES,
//			&numTotalTuples);
	lookupTable7Create2(refFileName, seedLen, maxNumHits, TUPLE_IGNORE_THRES,
			&numTotalTuples);
	int *keys, *values, numKeys, numValues, *numTuplesPerHash;
//	lookupTable5CpyHashTableToGPU(&keys, &numKeys, &values, &numValues,
//			&numTuplesPerHash);
	lookupTable7CpyHashTableToGPU(&keys, &numKeys, &values, &numValues,
			&numTuplesPerHash);
	fprintf(stderr, "done.\n");
	fprintf(stderr, "Total number of ref tuples = %d\n", numValues);

	/* Query the GPU for configuration. */
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, SELECTED_DEVICE_INDEX);
	PRINT_CUDA_ERROR()

	/* Find how many queries can fit on the GPU. */
	long globalMem = deviceProp.totalGlobalMem;
	long memHashTableKeys = pow((float) DNA_ALPHABET_SIZE, seedLen) * sizeof(int);
	long memHashTableVals = numTotalTuples * sizeof(int);
	long memHashTableNumTuples = memHashTableKeys;
	long memHashTable = memHashTableKeys + memHashTableVals
			+ memHashTableNumTuples;
	if (memHashTable >= globalMem)
	{
		fprintf(stderr, "Error: In Phase 1, hash table size (%ld B) is "
				"greater than the GPU global memory size (%ld B). Program is "
				"exiting.\n", memHashTable, globalMem);
		exit(1);
	}
	int memPerQryInput1 = (MAX_QRY_SEQ_LENGTH + 1) * sizeof(char);
	int memPerQryInput2 = sizeof(uchar);
	int memPerQryOutput = maxNumHits * (sizeof(char) + sizeof(int));
	/* Multiplying by 2 below to account for the reverse complement of a query. */
	int memPerQry = (2 * memPerQryInput1) + memPerQryInput2 + memPerQryOutput;
	int unusedGlobalMem = globalMem - memHashTable - GLOBAL_MEM_MARGIN;
	int numQrsPerIter = 0;
	/* Note: I'm not using 2 dimensional blocks because from experimentation I
	 * found that they don't seem to improve the speed at all. */
	int totalBlocks = deviceProp.maxGridSize[0];
	while (unusedGlobalMem > 0)
	{
		unusedGlobalMem -= memPerQry;
		++numQrsPerIter;
		if (numQrsPerIter >= totalBlocks)
			break;
	}
	if (numQrsPerIter % 2 == 1)
		--numQrsPerIter;
	if (numQrsPerIter < MIN_QRS_FOR_GPU)
	{
		fprintf(stderr, "Error: In Phase 1, there is not enough memory on the "
				"GPU global memory for %d queries.\n", MIN_QRS_FOR_GPU);
		exit(1);
	}
	int numBlocksXDim = numQrsPerIter;
	if (numBlocksXDim > deviceProp.maxGridSize[0])
		numBlocksXDim = deviceProp.maxGridSize[0];
	int numBlocksYDim = (int) ceil(((float) numQrsPerIter)
			/ deviceProp.maxGridSize[0]);
	fprintf(stderr, "Number of queries per iteration = %d\n", numQrsPerIter);
	fprintf(stderr, "Number of blocks: x = %d, y = %d\n", numBlocksXDim,
			numBlocksYDim);

	/* Allocate memory on the CPU. */
	char *qrs = (char *) malloc(numQrsPerIter * memPerQryInput1);
	uchar *qryLen = (uchar *) malloc(numQrsPerIter * sizeof(uchar));
	char *refIdx = (char *) malloc(maxNumHits * numQrsPerIter * sizeof(char));
	int *shift = (int *) malloc(maxNumHits * numQrsPerIter * sizeof(int));
	int *refPos = (int *) malloc(maxNumHits * numQrsPerIter * sizeof(int));
	int halfNumQrsPerIter = numQrsPerIter / 2;
	Query **qry = (Query **) malloc(halfNumQrsPerIter * sizeof(Query *));
	int i;
	for (i = 0; i < halfNumQrsPerIter; ++i)
		qry[i] = qryCreate2();

	/* Copy memory from host to constant memory on GPU. */
//	lookupTable5CpyConstMemToGPU();
	lookupTable7CpyConstMemToGPU();

	/* Allocate memory on the GPU. */
	char *qrs_d;
	cudaMalloc((void **) &qrs_d, numQrsPerIter * memPerQryInput1);
	PRINT_CUDA_ERROR()
	uchar *qryLen_d;
	cudaMalloc((void **) &qryLen_d, numQrsPerIter * sizeof(uchar));
	PRINT_CUDA_ERROR()
	char *refIdx_d;
	cudaMalloc((void **) &refIdx_d, maxNumHits * numQrsPerIter * sizeof(char));
	PRINT_CUDA_ERROR()
	int *refPos_d;
	cudaMalloc((void **) &refPos_d, maxNumHits * numQrsPerIter * sizeof(int));
	PRINT_CUDA_ERROR()

//	/* Find the ideal array size for the arrays in shared memory. */
//	int availShrMem = deviceProp.sharedMemPerBlock - SHR_MEM_MARGIN;
//	int currentShrMem = maxNumHits * sizeof(int);
//	int idealArrSize = 0;
//	int intSize = sizeof(int), charSize = sizeof(char);
//	int sizePerElement = (3 * intSize) + charSize;
//	do
//	{
//		currentShrMem += sizePerElement;
//		++idealArrSize;
//	} while (currentShrMem <= availShrMem);
//	currentShrMem = (idealArrSize * sizePerElement) + (maxNumHits * intSize);
//	currentShrMem += (intSize - (currentShrMem % intSize));
//	fprintf(stderr, "Array size on shared memory = %d\n", idealArrSize);
//	fprintf(stderr, "Dynamically allocated shared memory = %d\n", currentShrMem);

	/* Find the ideal array size for the arrays in shared memory. */
	int availShrMem = deviceProp.sharedMemPerBlock - SHR_MEM_MARGIN;
	int currentShrMem = maxNumHits * sizeof(int);
	int idealArrSize = 0;
	int intSize = sizeof(int), shortSize = sizeof(short), longSize = sizeof(long);
	int charSize = sizeof(char);
	int sizePerElement = longSize + intSize + shortSize;
	do
	{
		currentShrMem += sizePerElement;
		++idealArrSize;
	} while (currentShrMem <= availShrMem);
	currentShrMem = (idealArrSize * sizePerElement) + (maxNumHits * intSize);
	currentShrMem += (longSize - (currentShrMem % longSize));
	fprintf(stderr, "Array size on shared memory = %d\n", idealArrSize);
	fprintf(stderr, "Dynamically allocated shared memory = %d\n", currentShrMem);

	/* Search for matches. */
	FILE *qryFilePtr = fopen(qryFileName, "r");
	FILE *outputFilePtr = fopen(matchFile, "w");
	int j, k, randNum, isDone = FALSE, totalQrs = 0;
	dim3 dimGrid, dimBlock;
	dimGrid.x = numBlocksXDim;
	dimGrid.y = numBlocksYDim;
	dimBlock.x = deviceProp.maxThreadsDim[0];
	/* Make sure the number of threads launched is less than or equal to
	 * the number of elements in the arrays in the shared memory. */
	if (dimBlock.x > idealArrSize)
		dimBlock.x = idealArrSize;
	fprintf(stderr, "Number of threads per block = %u\n", dimBlock.x);
	dimBlock.y = 1;
	dimBlock.z = 1;
	*numMatches = 0;
	time_t startTime, endTime;
	double diffTime = 0.0;
	int numResultsPerIter, memResultsPerIter_char, memResultsPerIter_int;
	int numQryOpsCPU = 0;
	fprintf(stderr, "Filtering reference hits...\n");
	while (isDone == FALSE)
	{
		time(&startTime);

		/* Get query sequences. */
		fprintf(stderr, "   Fetching query sequences...");
		k = 0;
		for (j = 0; j < halfNumQrsPerIter; ++j)
		{
			if (qryGetNext2(qryFilePtr, qry[j]) == NULL)
			{
				isDone = TRUE;
				break;
			}
			strcpy(qrs + (k * MAX_QRY_SEQ_LENGTH), qry[j]->seq);
			qryLen[k] = qry[j]->length;
			++k;
			strcpy(qrs + (k * MAX_QRY_SEQ_LENGTH), qry[j]->revCompSeq);
			qryLen[k] = qry[j]->length;
			++k;
		}
		fprintf(stderr, "done.\n");
		numQrsPerIter = k;
		numResultsPerIter = maxNumHits * numQrsPerIter;
		memResultsPerIter_char = numResultsPerIter * charSize;
		memResultsPerIter_int = numResultsPerIter * intSize;

		/* Copy queries and their length to the GPU. */
		fprintf(stderr, "   Copying queries from CPU to GPU...");
		cudaMemcpy(qrs_d, qrs, numQrsPerIter * memPerQryInput1,
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qryLen_d, qryLen, numQrsPerIter * memPerQryInput2,
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Initialize output array on the GPU. */
		fprintf(stderr, "   Initializing output array on GPU...");
		cudaMemset(refIdx_d, CHAR_MAX, memResultsPerIter_char);
		PRINT_CUDA_ERROR()
		cudaMemset(refPos_d, -1, memResultsPerIter_int);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Call kernel function. */
		fprintf(stderr, "   Calling kernel function...");
		randNum = rand();
//		lookupTable5MapQry2_gpu<<<dimGrid, dimBlock, currentShrMem>>>(keys,
//				values, numTuplesPerHash, qrs_d, qryLen_d, MAX_QRY_SEQ_LENGTH,
//				refIdx_d, refPos_d, maxNumHits, seedLen, randNum, idealArrSize);
//		lookupTable7MapQry2_gpu<<<dimGrid, dimBlock, currentShrMem>>>(keys,
//				values, numTuplesPerHash, qrs_d, qryLen_d, MAX_QRY_SEQ_LENGTH,
//				refIdx_d, refPos_d, maxNumHits, seedLen, randNum, idealArrSize);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Copy results from GPU to CPU. */
		fprintf(stderr, "   Copying results from GPU to CPU...");
		cudaMemcpy(refIdx, refIdx_d, memResultsPerIter_char,
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos, refPos_d, memResultsPerIter_int,
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* If filtering could not be performed on the GPU due to
		 * memory limitations, do it on the CPU. */
		fprintf(stderr, "   Filtered ");
		numQryOpsCPU = 0;
		for (j = 0; j < numQrsPerIter; j += maxNumHits)
		{
			if (refIdx[j] == CHAR_MAX)
			{
				++numQryOpsCPU;
				lookupTable6MapQry2(qrs + (j * MAX_QRY_SEQ_LENGTH), qryLen[j],
						refIdx + j, shift + j, refPos + j);
			}
		}
		fprintf(stderr, "%d / %d queries on GPU\n",
				(numQrsPerIter - numQryOpsCPU), numQrsPerIter);

		/* Print results to file. */
		printHeuristicMatches5(qry, refIdx, refPos, numResultsPerIter,
				numMatches, outputFilePtr);

		/* Show progress. */
		time(&endTime);
		diffTime = difftime(endTime, startTime);
		totalQrs += numQrsPerIter;
		fprintf(stderr, "   Searched %d queries (total: %d)...%.2lf secs\n\n",
				numQrsPerIter, totalQrs, diffTime);
//		break;
	}

	/* Clean up. */
	free(qryFilePtr);
	free(outputFilePtr);
	free(qrs);
	free(qryLen);
	for (i = 0; i < halfNumQrsPerIter; ++i)
		qryDelete(qry[i]);
	free(qry);
	free(refIdx);
	free(shift);
	free(refPos);
	cudaFree(qrs_d);
	cudaFree(qryLen_d);
	cudaFree(refIdx_d);
	cudaFree(refPos_d);
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param	qryFileName	Query file name.
 * @param	refFileName	Reference file name.
 * @param	matchFile	File in which filtered reference sequence fragment
 * positions are stored.
 * @param	seedLen		Seed length.
 * @param	maxNumHits	Max number of hits per query.
 * @param 	numMatches	Total number of matches obtained.
 */
void searchQueries4_2(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches)
{
	if (qryFileName == NULL || refFileName == NULL || matchFile == NULL
			|| seedLen <= 0)
		return;

	/* Preprocess reference */
	fprintf(stderr, "Preprocessing reference sequences...\n");
	int numTotalTuples;
	lookupTable7Create2(refFileName, seedLen, maxNumHits, TUPLE_IGNORE_THRES,
			&numTotalTuples);
	int *keys_d, *vals_d, numKeys, numVals, *numRptsPerTuple_d;
	lookupTable7CpyHashTableToGPU(&keys_d, &numKeys, &vals_d, &numVals,
			&numRptsPerTuple_d);
	fprintf(stderr, "done.\n");
	fprintf(stderr, "Total number of distinct tuples = %d\n", numKeys);
	fprintf(stderr, "Total number of ref tuples = %d\n", numVals);

	lookupTable7CpyConstMemToGPU();

	/* Allocate memory on the GPU. */
	int numQrsPerIter = NUM_BLOCKS_PHASE1;
	char *qrs_d;
	cudaMalloc((void **) &qrs_d,
			numQrsPerIter * (MAX_QRY_SEQ_LENGTH + 1) * sizeof(char));
	PRINT_CUDA_ERROR()
	uchar *qryLen_d;
	cudaMalloc((void **) &qryLen_d, numQrsPerIter * sizeof(uchar));
	PRINT_CUDA_ERROR()
	char *refIdx_d;
	cudaMalloc((void **) &refIdx_d, maxNumHits * numQrsPerIter * sizeof(char));
	PRINT_CUDA_ERROR()
	int *refPos_d;
	cudaMalloc((void **) &refPos_d, maxNumHits * numQrsPerIter * sizeof(int));
	PRINT_CUDA_ERROR()

	/* Query the GPU for configuration. */
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, SELECTED_DEVICE_INDEX);
	PRINT_CUDA_ERROR()

	/* Find the ideal array size for the arrays in shared memory. */
	int availShrMem = deviceProp.sharedMemPerBlock - SHR_MEM_MARGIN;
	int maxNumChr = MAX_NUM_REFS;
	int extraMem = (maxNumHits * sizeof(int)) + (maxNumChr * sizeof(short))
			+ (maxNumChr * sizeof(short));
	int currentShrMem = extraMem;
	int idealArrSize = 0;
	int intSize = sizeof(int), shortSize = sizeof(short), charSize = sizeof(char);
	int sizePerElement = charSize + (3 * intSize) + shortSize;
	do
	{
		currentShrMem += sizePerElement;
		++idealArrSize;
	} while (currentShrMem <= availShrMem);
	currentShrMem = (idealArrSize * sizePerElement) + extraMem;
	currentShrMem += (intSize - (currentShrMem % intSize));
	fprintf(stderr, "Array size on shared memory = %d\n", idealArrSize);
	fprintf(stderr, "Dynamically allocated shared memory = %d\n", currentShrMem);

	int halfNumQrsPerIter = numQrsPerIter / 2;
	Query **qry = (Query **) malloc(halfNumQrsPerIter * sizeof(Query *));
	int i;
	for (i = 0; i < halfNumQrsPerIter; ++i)
		qry[i] = qryCreate2();

	char *qrs = (char *) malloc(
			numQrsPerIter * (MAX_QRY_SEQ_LENGTH + 1) * sizeof(char));
	uchar *qryLen = (uchar *) malloc(numQrsPerIter * sizeof(uchar));
	char *refIdx = (char *) malloc(maxNumHits * numQrsPerIter * sizeof(char));
	int *shift = (int *) malloc(maxNumHits * numQrsPerIter * sizeof(int));
	int *refPos = (int *) malloc(maxNumHits * numQrsPerIter * sizeof(int));
	int isDone = FALSE;
	int j, k, randNum;
	int numResultsPerIter, memResultsPerIter_char, memResultsPerIter_int;
	FILE *qryFilePtr = fopen(qryFileName, "r");
	dim3 dimGrid = numQrsPerIter, dimBlock = maxNumChr;
	srand(time(NULL));
	FILE *outputFilePtr = fopen(matchFile, "w");
	time_t startTime, endTime;
	double diffTime = 0.0, totalTime = 0.0;
	int totalQrs = 0, numQryOpsCPU = 0;
	fprintf(stderr, "Filtering reference hits...\n");
	while (isDone == FALSE)
	{
		time(&startTime);

		/* Get query sequences. */
		fprintf(stderr, "   Fetching query sequences...");
		k = 0;
		for (j = 0; j < halfNumQrsPerIter; ++j)
		{
			if (qryGetNext2(qryFilePtr, qry[j]) == NULL)
			{
				isDone = TRUE;
				break;
			}
			strcpy(qrs + (k * MAX_QRY_SEQ_LENGTH), qry[j]->seq);
			qryLen[k] = qry[j]->length;
			++k;
			strcpy(qrs + (k * MAX_QRY_SEQ_LENGTH), qry[j]->revCompSeq);
			qryLen[k] = qry[j]->length;
			++k;
		}
		fprintf(stderr, "done.\n");
		numQrsPerIter = k;
		numResultsPerIter = maxNumHits * numQrsPerIter;
		memResultsPerIter_char = numResultsPerIter * sizeof(char);
		memResultsPerIter_int = numResultsPerIter * sizeof(int);

		/* Copy queries and their length to the GPU. */
		fprintf(stderr, "   Copying queries from CPU to GPU...");
		cudaMemcpy(qrs_d, qrs,
				numQrsPerIter * (MAX_QRY_SEQ_LENGTH + 1) * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qryLen_d, qryLen, numQrsPerIter * sizeof(uchar),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Initialize output array on the GPU. */
		fprintf(stderr, "   Initializing output array on GPU...");
		cudaMemset(refIdx_d, 0, memResultsPerIter_char);
		PRINT_CUDA_ERROR()
		cudaMemset(refPos_d, 0, memResultsPerIter_int);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Call kernel function. */
		fprintf(stderr, "   Calling kernel function...");
		randNum = rand();
		lookupTable7MapQry2_gpu2<<<dimGrid, dimBlock, currentShrMem>>>(keys_d,
				vals_d, numRptsPerTuple_d, qrs_d, qryLen_d, MAX_QRY_SEQ_LENGTH,
				refIdx_d, refPos_d, maxNumHits, seedLen, randNum, idealArrSize,
				maxNumChr);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Copy results from GPU to CPU. */
		fprintf(stderr, "   Copying results from GPU to CPU...");
		cudaMemcpy(refIdx, refIdx_d, memResultsPerIter_char,
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos, refPos_d, memResultsPerIter_int,
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* If filtering could not be performed on the GPU due to
		 * memory limitations, do it on the CPU. */
		fprintf(stderr, "   Filtered ");
		numQryOpsCPU = 0;
		for (j = 0; j < numQrsPerIter; j += maxNumHits)
		{
//			if (refIdx[j] == -1)
//			{
				++numQryOpsCPU;
				lookupTable7MapQry2(qrs + (j * MAX_QRY_SEQ_LENGTH), qryLen[j],
						refIdx + j, shift + j, refPos + j);
//			}
		}
		fprintf(stderr, "%d / %d queries on GPU\n",
				(numQrsPerIter - numQryOpsCPU), numQrsPerIter);

		/* Print results to file. */
		printHeuristicMatches5_2(qry, refIdx, refPos, numQrsPerIter, maxNumHits,
				numMatches, outputFilePtr);

		/* Show progress. */
		time(&endTime);
		diffTime = difftime(endTime, startTime);
		totalTime += diffTime;
		totalQrs += numQrsPerIter;
		fprintf(stderr, "   Searched %d queries (total: %d)...%.2lf secs "
				"(total: %.2lf secs)\n\n", numQrsPerIter, totalQrs, diffTime,
				totalTime);

		break;
	}
	fclose(qryFilePtr);
	fclose(outputFilePtr);
	for (i = 0; i < halfNumQrsPerIter; ++i)
		free(qry[i]);
	free(qry);
	free(qrs);
	free(qryLen);
	free(refIdx);
	free(refPos);
	cudaFree(keys_d);
	cudaFree(vals_d);
	cudaFree(numRptsPerTuple_d);
	cudaFree(qrs_d);
	cudaFree(qryLen_d);
	cudaFree(refIdx_d);
	cudaFree(refPos_d);
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param	qryFileName	Query file name.
 * @param	refFileName	Reference file name.
 * @param	matchFile	File in which filtered reference sequence fragment
 * positions are stored.
 * @param	seedLen		Seed length.
 * @param	maxNumHits	Max number of hits per query.
 * @param 	numMatches	Total number of matches obtained.
 * @param	numThreads	Number of threads.
 * @param	tupleIgnoreThres	Threshold value used to ignore reference
 * tuples.
 */
void searchQueries4_3(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches, int numThreads,
		int tupleIgnoreThres)
{
	if (qryFileName == NULL || refFileName == NULL || matchFile == NULL
			|| seedLen <= 0)
		return;

	/* Preprocess reference */
	fprintf(stderr, "Preprocessing reference sequences...\n");
	int numTotalTuples;
	lookupTable7Create2(refFileName, seedLen, maxNumHits, tupleIgnoreThres,
			&numTotalTuples);
//	lookupTable7Create3("/home/pgupta/data/gpusw/keys100.txt",
//			"/home/pgupta/data/gpusw/vals100.txt", maxNumHits, seedLen);
	int *keys_d, *keys, *vals_d, *vals, numKeys, numVals, *numRptsPerTuple_d;
	int *numRptsPerTuple;
	lookupTable7CpyHashTableToGPU(&keys_d, &numKeys, &vals_d, &numVals,
			&numRptsPerTuple_d);
//	lookupTable7FetchHashTable(&keys, &numKeys, &vals, &numVals,
//			&numRptsPerTuple);
	fprintf(stderr, "done.\n");
	fprintf(stderr, "Total number of distinct tuples = %d\n", numKeys);
	fprintf(stderr, "Total number of ref tuples = %d\n", numVals);

	lookupTable7CpyConstMemToGPU();
//	lookupTable7CreateConstMem();

	/* Allocate memory on the GPU. */
	int numQrsPerIter = NUM_BLOCKS_PHASE1;
	char *qrs_d;
	cudaMalloc((void **) &qrs_d,
			numQrsPerIter * (MAX_QRY_SEQ_LENGTH + 1) * sizeof(char));
	PRINT_CUDA_ERROR()
	uchar *qryLen_d;
	cudaMalloc((void **) &qryLen_d, numQrsPerIter * sizeof(uchar));
	PRINT_CUDA_ERROR()
	char *refIdx_d;
	cudaMalloc((void **) &refIdx_d, maxNumHits * numQrsPerIter * sizeof(char));
	PRINT_CUDA_ERROR()
	int *refPos_d;
	cudaMalloc((void **) &refPos_d, maxNumHits * numQrsPerIter * sizeof(int));
	PRINT_CUDA_ERROR()

	/* Query the GPU for configuration. */
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, SELECTED_DEVICE_INDEX);
	PRINT_CUDA_ERROR()

	/* Find the ideal array size for the arrays in shared memory. */
	int availShrMem = deviceProp.sharedMemPerBlock - SHR_MEM_MARGIN;
	int maxNumChr = MAX_NUM_REFS;
	int extraMem = (maxNumHits * sizeof(int)) + (maxNumChr * sizeof(short))
			+ (maxNumChr * sizeof(short));
	int currentShrMem = extraMem;
	int idealArrSize = 0;
	int intSize = sizeof(int), shortSize = sizeof(short), charSize = sizeof(char);
	int sizePerElement = charSize + (3 * intSize) + shortSize;
	do
	{
		currentShrMem += sizePerElement;
		++idealArrSize;
	} while (currentShrMem <= availShrMem);
	currentShrMem = (idealArrSize * sizePerElement) + extraMem;
	currentShrMem += (intSize - (currentShrMem % intSize));
	fprintf(stderr, "Array size on shared memory = %d\n", idealArrSize);
	fprintf(stderr, "Dynamically allocated shared memory = %d\n", currentShrMem);

	int halfNumQrsPerIter = numQrsPerIter / 2;
	Query **qry = (Query **) malloc(halfNumQrsPerIter * sizeof(Query *));
	int i;
	for (i = 0; i < halfNumQrsPerIter; ++i)
		qry[i] = qryCreate2();

	char *qrs = (char *) malloc(
			numQrsPerIter * (MAX_QRY_SEQ_LENGTH + 1) * sizeof(char));
	uchar *qryLen = (uchar *) malloc(numQrsPerIter * sizeof(uchar));
	char *refIdx = (char *) malloc(maxNumHits * numQrsPerIter * sizeof(char));
	int *shift = (int *) malloc(maxNumHits * numQrsPerIter * sizeof(int));
	int *refPos = (int *) malloc(maxNumHits * numQrsPerIter * sizeof(int));

	int isDone = FALSE;
	int j, k, l, m, n, randNum;
	int numResultsPerIter, memResultsPerIter_char, memResultsPerIter_int;
	FILE *qryFilePtr = fopen(qryFileName, "r");
	dim3 dimGrid = numQrsPerIter, dimBlock = maxNumChr;
	srand(time(NULL));
	FILE *outputFilePtr = fopen(matchFile, "w");
	time_t startTime, endTime;
	double diffTime = 0.0, totalTime = 0.0;
	int totalQrs = 0, numQryOpsCPU = 0, numQrsFiltered;
	*numMatches = 0;
	fprintf(stderr, "Filtering reference hits...\n");
	int count = 0;
	while (isDone == FALSE)
	{
		time(&startTime);

		/* Get query sequences. */
		fprintf(stderr, "   Fetching query sequences...");
		k = 0;
		for (j = 0; j < halfNumQrsPerIter; ++j)
		{
			if (qryGetNext2(qryFilePtr, qry[j]) == NULL)
			{
				isDone = TRUE;
				break;
			}
			strcpy(qrs + (k * MAX_QRY_SEQ_LENGTH), qry[j]->seq);
			qryLen[k] = qry[j]->length;
			++k;
			strcpy(qrs + (k * MAX_QRY_SEQ_LENGTH), qry[j]->revCompSeq);
			qryLen[k] = qry[j]->length;
			++k;
		}
		fprintf(stderr, "done.\n");
		numQrsPerIter = k;
		numResultsPerIter = maxNumHits * numQrsPerIter;
		memResultsPerIter_char = numResultsPerIter * sizeof(char);
		memResultsPerIter_int = numResultsPerIter * sizeof(int);

		/* Copy queries and their length to the GPU. */
		fprintf(stderr, "   Copying queries from CPU to GPU...");
		cudaMemcpy(qrs_d, qrs,
				numQrsPerIter * (MAX_QRY_SEQ_LENGTH + 1) * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qryLen_d, qryLen, numQrsPerIter * sizeof(uchar),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Call kernel function. */
		fprintf(stderr, "   Calling kernel function...");
		randNum = rand();
		lookupTable7MapQry2_gpu2<<<dimGrid, dimBlock, currentShrMem>>>(keys_d,
				vals_d, numRptsPerTuple_d, qrs_d, qryLen_d, MAX_QRY_SEQ_LENGTH,
				refIdx_d, refPos_d, maxNumHits, seedLen, randNum, idealArrSize,
				maxNumChr);
		PRINT_CUDA_ERROR()
//		for (j = 0; j < numQrsPerIter; ++j)
//		{
//			lookupTable7MapQry2_cpu2(keys, vals, numRptsPerTuple,
//					qrs, qryLen, MAX_QRY_SEQ_LENGTH, refIdx, refPos,
//					maxNumHits, seedLen, randNum, idealArrSize, maxNumChr, j);
//		}
		fprintf(stderr, "done.\n");

		/* Copy results from GPU to CPU. */
		fprintf(stderr, "   Copying results from GPU to CPU...");
		cudaMemcpy(refIdx, refIdx_d, memResultsPerIter_char,
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos, refPos_d, memResultsPerIter_int,
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* If filtering could not be performed on the GPU due to
		 * memory limitations, do it on the CPU. */
		fprintf(stderr, "   Filtered ");
		numQryOpsCPU = 0;
		for (j = 0; j < numQrsPerIter; ++j)
		{
			k = j * maxNumHits;
			if (refIdx[k] == -1)
			{
				++numQryOpsCPU;
//				lookupTable7MapQry2(qrs + (j * MAX_QRY_SEQ_LENGTH), qryLen[j],
//						refIdx + k, shift + k, refPos + k);
			}
		}
		numQrsFiltered = numQrsPerIter - numQryOpsCPU;
		fprintf(stderr, "%d / %d queries on GPU\n", numQrsFiltered,
				numQrsPerIter);

		/* Print results to file. */
		printHeuristicMatches5_2(qry, refIdx, refPos, numQrsPerIter, maxNumHits,
				numMatches, outputFilePtr);

		/* Show progress. */
		time(&endTime);
		diffTime = difftime(endTime, startTime);
		totalTime += diffTime;
		totalQrs += numQrsPerIter;
		fprintf(stderr, "   Searched %d queries (total: %d)...%.2lf secs "
				"(total: %.2lf secs)\n\n", numQrsFiltered, totalQrs, diffTime,
				totalTime);
		++count;
		if (count == 1)
			break;
	}
	fprintf(stderr, "total numMatches = %d\n", *numMatches);
	fclose(qryFilePtr);
	fclose(outputFilePtr);
	for (i = 0; i < halfNumQrsPerIter; ++i)
		free(qry[i]);
	free(qry);
	free(qrs);
	free(qryLen);
	free(refIdx);
	free(refPos);
	cudaFree(keys_d);
	cudaFree(vals_d);
	cudaFree(numRptsPerTuple_d);
	cudaFree(qrs_d);
	cudaFree(qryLen_d);
	cudaFree(refIdx_d);
	cudaFree(refPos_d);
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param	qryFileName	Query file name.
 * @param	refFileName	Reference file name.
 * @param	matchFile	File in which filtered reference sequence fragment
 * positions are stored.
 * @param	seedLen		Seed length.
 * @param	maxNumHits	Max number of hits per query.
 * @param 	numMatches	Total number of matches obtained.
 * @param	tupleIgnoreThres	Threshold value used to ignore reference
 * tuples.
 */
void searchQueries5(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches, int tupleIgnoreThres)
{
	if (qryFileName == NULL || refFileName == NULL || matchFile == NULL
			|| seedLen <= 0)
		return;

	/* Preprocess reference */
	fprintf(stderr, "Preprocessing reference sequences...\n");
	lookupTable8Create(refFileName, seedLen, maxNumHits, tupleIgnoreThres);
	fprintf(stderr, "done.\n");

	Query *qry = (Query *) malloc(sizeof(Query));
	qry = qryCreate2();
	char *refIdx = (char *) malloc(maxNumHits * sizeof(char));
	int *refPos = (int *) malloc(maxNumHits * sizeof(int));
	int numHits = 0, numQrs = 0, isDone = FALSE, totalQrs = 0;
	FILE *qryFilePtr = fopen(qryFileName, "r");
	FILE *outputFilePtr = fopen(matchFile, "w");
	srand(time(NULL));
	time_t startTime, endTime;
	double diffTime = 0.0, totalTime = 0.0;
	*numMatches = 0;
	fprintf(stderr, "Filtering reference hits...\n");
	while (isDone == FALSE)
	{
		time(&startTime);

		/* Get query sequences. */
		if (qryGetNext2(qryFilePtr, qry) == NULL)
		{
			isDone = TRUE;
			break;
		}

		/* Filter query hits. */
		numHits = lookupTable8MapQry(qry->seq, qry->length, refIdx, refPos);
		*numMatches += numHits;
		printHeuristicMatches6(qry, FORWARD_STRAND, refIdx, refPos, numHits,
				outputFilePtr);
		++numQrs;

		numHits = lookupTable8MapQry(qry->revCompSeq, qry->length, refIdx,
				refPos);
		*numMatches += numHits;
		printHeuristicMatches6(qry, REV_STRAND, refIdx, refPos, numHits,
				outputFilePtr);
		++numQrs;

		/* Show progress. */
		time(&endTime);
		diffTime = difftime(endTime, startTime);
		totalTime += diffTime;
		if (numQrs % 55000 == 0)
		{
			totalQrs += numQrs;
			fprintf(stderr, "   Searched %d queries in %.2lf secs "
					"(total qrs: %d; total time: %.2lf secs)\n\n",
					numQrs, diffTime, totalQrs, totalTime);
			break;
		}
	}
	fclose(qryFilePtr);
	fclose(outputFilePtr);
	free(qry);
	free(refIdx);
	free(refPos);
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param	qryFileName	Query file name.
 * @param	refFileName	Reference file name.
 * @param	matchFile	File in which filtered reference sequence fragment
 * positions are stored.
 * @param	seedLen		Seed length.
 * @param	maxNumHits	Max number of hits per query.
 * @param 	numMatches	Total number of matches obtained.
 * @param	tupleIgnoreThres	Threshold value used to ignore reference
 * tuples.
 */
void searchQueries9(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches, int tupleIgnoreThres)
{
	if (qryFileName == NULL || refFileName == NULL || matchFile == NULL
			|| seedLen <= 0)
		return;

	/* Preprocess reference */
	fprintf(stderr, "Preprocessing reference sequences...\n");
//	lookupTable9Create(refFileName, seedLen, maxNumHits, tupleIgnoreThres);
	lookupTable9Create("/home/pgupta/data/gpusw/keys100.txt",
			"/home/pgupta/data/gpusw/vals100.txt", maxNumHits, seedLen,
			tupleIgnoreThres);
	fprintf(stderr, "done.\n");

	Query *qry = (Query *) malloc(sizeof(Query));
	qry = qryCreate2();
	char *refIdx = (char *) malloc(maxNumHits * sizeof(char));
	int *refPos = (int *) malloc(maxNumHits * sizeof(int));
	int numHits = 0, numQrs = 0, isDone = FALSE, totalQrs = 0;
	FILE *qryFilePtr = fopen(qryFileName, "r");
	FILE *outputFilePtr = fopen(matchFile, "w");
	srand(time(NULL));
	time_t startTime, endTime;
	double diffTime = 0.0, totalTime = 0.0;
	*numMatches = 0;
	fprintf(stderr, "Filtering reference hits...\n");
	while (isDone == FALSE)
	{
		time(&startTime);

		/* Get query sequences. */
		if (qryGetNext2(qryFilePtr, qry) == NULL)
		{
			isDone = TRUE;
			break;
		}

		/* Filter query hits. */
		numHits = lookupTable9MapQry(qry->seq, qry->length, refIdx, refPos);
		*numMatches += numHits;
		printHeuristicMatches6(qry, FORWARD_STRAND, refIdx, refPos, numHits,
				outputFilePtr);
		++numQrs;

//		numHits = lookupTable9MapQry(qry->revCompSeq, qry->length, refIdx,
//				refPos);
//		*numMatches += numHits;
//		printHeuristicMatches6(qry, REV_STRAND, refIdx, refPos, numHits,
//				outputFilePtr);
//		++numQrs;

		/* Show progress. */
		time(&endTime);
		diffTime = difftime(endTime, startTime);
		totalTime += diffTime;
		if (numQrs == 1)
		{
			totalQrs += numQrs;
			fprintf(stderr, "   Searched %d queries in %.2lf secs "
					"(total qrs: %d; total time: %.2lf secs)\n\n",
					numQrs, diffTime, totalQrs, totalTime);
			break;
		}
	}
	fclose(qryFilePtr);
	fclose(outputFilePtr);
	free(qry);
	free(refIdx);
	free(refPos);
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param	qryFileName	Query file name.
 * @param	refFileName	Reference file name.
 * @param	matchFile	File in which filtered reference sequence fragment
 * positions are stored.
 * @param	seedLen		Seed length.
 * @param	maxNumHits	Max number of hits per query.
 * @param[out] 	numMatches	Total number of matches obtained.
 * @param	tupleIgnoreThres	Threshold value used to ignore reference
 * tuples.
 */
void searchQueries10(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches, int tupleIgnoreThres)
{
	if (qryFileName == NULL || refFileName == NULL || matchFile == NULL
			|| seedLen <= 0)
		return;

	/* Preprocess reference */
	fprintf(stderr, "Preprocessing reference sequences...\n");
	lookupTable10Create(refFileName, seedLen, maxNumHits, tupleIgnoreThres);
//	lookupTable10Create("/home/pgupta/data/gpusw/keys100_12M5.txt",
//			"/home/pgupta/data/gpusw/vals100.txt", maxNumHits, seedLen,
//			tupleIgnoreThres);
	fprintf(stderr, "done.\n");

	Query *qry = (Query *) malloc(sizeof(Query));
	qry = qryCreate2();
	char *refIdx = (char *) malloc(maxNumHits * sizeof(char));
	int *refPos = (int *) malloc(maxNumHits * sizeof(int));
	int *refSize = (int *) malloc(maxNumHits * sizeof(int));
	int numHits = 0, numQrs = 0, isDone = FALSE, totalQrs = 0;
	FILE *qryFilePtr = fopen(qryFileName, "r");
	FILE *outputFilePtr = fopen(matchFile, "w");
	srand(time(NULL));
	time_t startTime, endTime;
	double diffTime = 0.0, totalTime = 0.0;
	*numMatches = 0;
	fprintf(stderr, "Filtering reference hits...\n");
	while (isDone == FALSE)
	{
		time(&startTime);

		/* Get query sequences. */
		if (qryGetNext2(qryFilePtr, qry) == NULL)
		{
			isDone = TRUE;
			break;
		}

		/* Filter query hits. */
		numHits = lookupTable10MapQry(qry->seq, qry->length, refIdx, refPos,
				refSize);
		printHeuristicMatches10(qry, FORWARD_STRAND, refIdx, refPos, refSize,
				numHits, numMatches, outputFilePtr);
		++numQrs;

		numHits = lookupTable10MapQry(qry->revCompSeq, qry->length, refIdx,
				refPos, refSize);
		printHeuristicMatches10(qry, REV_STRAND, refIdx, refPos, refSize,
				numHits, numMatches, outputFilePtr);
		++numQrs;

		/* Show progress. */
		time(&endTime);
		diffTime = difftime(endTime, startTime);
		totalTime += diffTime;
		if (numQrs % 100000 == 0)
		{
			totalQrs += numQrs;
			fprintf(stderr, "   Searched %d queries in %.2lf secs "
					"(total qrs: %d; total time: %.2lf secs)\n\n",
					numQrs, diffTime, totalQrs, totalTime);
			numQrs = 0;
//			break;
		}
	}
	fclose(qryFilePtr);
	fclose(outputFilePtr);
	free(qry);
	free(refIdx);
	free(refPos);
	lookupTable10Delete();
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param	qryFileName	Query file name.
 * @param	refFileName	Reference file name.
 * @param	matchFile	File in which filtered reference sequence fragment
 * positions are stored.
 * @param	seedLen		Seed length.
 * @param	maxNumHits	Max number of hits per query.
 * @param[out] 	numMatches	Total number of matches obtained.
 * @param	tupleIgnoreThres	Threshold value used to ignore reference
 * tuples.
 */
void searchQueries11(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches, int tupleIgnoreThres)
{
	if (qryFileName == NULL || refFileName == NULL || matchFile == NULL
			|| seedLen <= 0)
		return;

	/* Preprocess reference */
	fprintf(stderr, "Preprocessing reference sequences...\n");
//	lookupTable11Create(refFileName, seedLen, maxNumHits, tupleIgnoreThres);
//	lookupTable11Create2(refFileName, seedLen, maxNumHits, tupleIgnoreThres);
	lookupTable11Create("/home/pgupta/data/gpusw/keys1000.txt",
			"/home/pgupta/data/gpusw/vals1000.txt",
			"/home/pgupta/data/gpusw/keys1000_2.txt", maxNumHits, seedLen,
			tupleIgnoreThres);
	fprintf(stderr, "done.\n");

	Query *qry = (Query *) malloc(sizeof(Query));
	qry = qryCreate2();
	char *refIdx = (char *) malloc(maxNumHits * sizeof(char));
	int *refPos = (int *) malloc(maxNumHits * sizeof(int));
	int numHits = 0, numQrs = 0, isDone = FALSE, totalQrs = 0;
	FILE *qryFilePtr = fopen(qryFileName, "r");
	FILE *outputFilePtr = fopen(matchFile, "w");
	srand(time(NULL));
	time_t startTime, endTime;
	double diffTime = 0.0, totalTime = 0.0;
	*numMatches = 0;
	int x = 7727;
	fprintf(stderr, "Filtering reference hits...\n");
	while (isDone == FALSE)
	{
		time(&startTime);

		/* Get query sequences. */
		if (qryGetNext2(qryFilePtr, qry) == NULL)
		{
			isDone = TRUE;
			break;
		}
//		while (--x > 0)
//			qryGetNext2(qryFilePtr, qry);

		/* Filter query hits. */
//		numHits = lookupTable11MapQry(qry->seq, qry->length, refIdx, refPos);
		numHits = lookupTable11MapQry2(qry->seq, qry->length, refIdx, refPos);
		printHeuristicMatches11(qry, FORWARD_STRAND, refIdx, refPos, numHits,
				numMatches, outputFilePtr);
		++numQrs;

//		numHits = lookupTable11MapQry(qry->revCompSeq, qry->length, refIdx,
//				refPos);
		numHits = lookupTable11MapQry2(qry->revCompSeq, qry->length, refIdx,
				refPos);
		printHeuristicMatches11(qry, REV_STRAND, refIdx, refPos, numHits,
				numMatches, outputFilePtr);
		++numQrs;

		/* Show progress. */
		time(&endTime);
		diffTime = difftime(endTime, startTime);
		totalTime += diffTime;
		if (numQrs % 50000 == 0)
		{
			totalQrs += numQrs;
			fprintf(stderr, "   Searched %d queries in %.2lf secs "
					"(total qrs: %d; total time: %.2lf secs)\n",
					numQrs, diffTime, totalQrs, totalTime);
			numQrs = 0;
//			break;
		}
	}
	fclose(qryFilePtr);
	fclose(outputFilePtr);
	free(qry);
	free(refIdx);
	free(refPos);
	lookupTable10Delete();
}


/**
 * Get hits for forward strand.
 *
 * @param arg	An argument.
 */
void *getHits(void *arg)
{
	Data *data = (Data *) arg;
	data->numHits1 = lookupTable5MapQry(data->qry, data->qry->length,
			FORWARD_STRAND, data->refIdx1, data->shift1, data->refPos1);
	data->numHits2 = lookupTable5MapQry(data->qry, data->qry->length,
			REV_STRAND, data->refIdx2, data->shift2, data->refPos2);
	data->numHits = data->numHits1 + data->numHits2;
	return NULL;
}


/**
 * Get hits for forward strand.
 *
 * @param arg	An argument.
 */
void *getHits2(void *arg)
{
	Data2 *data = (Data2 *) arg;
	lookupTable7MapQry2(data->qrySeq, data->qryLen, data->refIdx, data->shift,
			data->refPos);
	return NULL;
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param 	qryFile		Query file name.
 * @param 	numQrs 		Number of queries in the query file.
 * @param 	qryLen 		Length of each query in the query file.
 * @param 	refFile		Reference file name.
 * @param 	numRefs 	Number of references in the reference file.
 * @param 	refLen 		Length of each reference sequence in the reference file.
 * @param 	seedLen 	Seed length.
 * @param 	matchFile 	File in which filtered reference sequence fragment
 * positions are stored.
 * @param 	maxNumHits 	Max number of hits per query.
 */
void searchQueries_gpu(char *qryFile, int numQrs, int qryLen, char *refFile,
		int numRefs, int refLen, int seedLen, char *matchFile, uint maxNumHits)
{
	if (qryFile == NULL || numQrs <= 0 || qryLen == 0
			|| refFile == NULL || numRefs <= 0 || refLen <= 0
			|| seedLen <= 0 || matchFile == NULL)
		return;

	/* Set the GPU device. */
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	cudaDeviceProp deviceProp[deviceCount];
	int device;
	for (device = 0; device < deviceCount; ++device)
		cudaGetDeviceProperties(&deviceProp[device], device);
	cudaGetDevice(&device);
	PRINT_CUDA_ERROR()
	if (device != SELECTED_DEVICE_INDEX)
	{
		cudaSetDevice(SELECTED_DEVICE_INDEX);
		PRINT_CUDA_ERROR()
	}

	/* Preprocess reference. */
	fprintf(stderr, "   Creating hash table...");
	RefPosList **hashTable = refPreprocess(refFile, seedLen);
	fprintf(stderr, "done.\n");

	/* Serialize the hash table. */
	fprintf(stderr, "   Serializing the hash table...");
	int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, seedLen) + 1;
	int numKeys2 = numKeys - 1;
	int i;
	uint numVals = 0;
	RefPosList *tmpRefPosList;
	for (i = 0; i < numKeys2; ++i)
	{
		tmpRefPosList = hashTable[i];
		while (tmpRefPosList != NULL)
		{
			++numVals;
			tmpRefPosList = tmpRefPosList->next;
		}
	}
	uint *keys = (uint *) malloc(numKeys * sizeof(uint));
	uint *vals = (uint *) malloc(numVals * sizeof(uint));
	for (i = 0; i < numVals; ++i)
		vals[i] = UINT_MAX;
	uint tupleIgnoreThres = TUPLE_IGNORE_THRES;
	serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals, seedLen,
			tupleIgnoreThres);
	fprintf(stderr, "done.\n");

	/* Map local positions (reference index and tuple index) to global
	 * positions (file offset). */
	fprintf(stderr, "   Creating position map...");
	refMapCreate(hashTable, numKeys2, numVals, seedLen);
	fprintf(stderr, "done.\n");
	refDeleteHashTable(hashTable, numKeys2);

	/* Map sequence offsets to name offsets. */
	fprintf(stderr, "   Creating name map...");
	refNameMapCreate(refFile);
	fprintf(stderr, "done.\n");

	/* Copy serialized hash table to GPU. */
	fprintf(stderr, "   Copying serialized hash table to GPU...");
	uint *keys_d;
	cudaMalloc(&keys_d, numKeys * sizeof(uint));
	PRINT_CUDA_ERROR()
	cudaMemcpy(keys_d, keys, numKeys * sizeof(uint), cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()
	uint *vals_d;
	cudaMalloc(&vals_d, numVals * sizeof(uint));
	PRINT_CUDA_ERROR()
	cudaMemcpy(vals_d, vals, numVals * sizeof(uint), cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()
	fprintf(stderr, "done.\n");

	/* Calculate number of iterations, number of queries per iteration, and
	 * number of queries in the last iteration. */
	fprintf(stderr, "   Calculating num of iters and num queries per iter...");
	int numIter;
	long numQrySeqPerIter, numQrySeqLastIter;
	int keysTotalMem = numKeys * sizeof(uint);
	int valsTotalMem = numVals * sizeof(uint);
	int qrySeqMem = qryLen * sizeof(char);
	int refIdxMem = sizeof(char);
	int shiftMem = sizeof(int);
	int posMem = sizeof(int);
	long gpuTotalMem = deviceProp[SELECTED_DEVICE_INDEX].totalGlobalMem;

	getNumIters(&numIter, &numQrySeqPerIter, &numQrySeqLastIter, numQrs * 2,
			maxNumHits, gpuTotalMem, GLOBAL_MEM_MARGIN, keysTotalMem,
			valsTotalMem, qrySeqMem, refIdxMem, shiftMem, posMem);
	long numQrsPerIter = numQrySeqPerIter / 2;

	fprintf(stderr, "done.\n");

	fprintf(stderr, "      Number of iterations = %d\n", numIter);
	fprintf(stderr, "      Number of query seqs per iter = %ld\n",
			numQrySeqPerIter);
	fprintf(stderr, "      Number of query seqs in last iter = %ld\n",
			numQrySeqLastIter);
	fprintf(stderr, "      Number of queries = %d\n", numQrs);

	/* Allocate CPU and GPU memory for query sequences. */
	long totalQrySeqLen = (numQrySeqPerIter * qryLen * sizeof(char))
			+ sizeof(char);
	char *qrySeq_d;
	cudaMalloc(&qrySeq_d, totalQrySeqLen);
	PRINT_CUDA_ERROR()

	Query **qry = (Query **) malloc(numQrsPerIter * sizeof(Query *));
	for (i = 0; i < numQrsPerIter; ++i)
		qry[i] = qryCreate2();
	FILE *qryFilePtr = fopen(qryFile, "r");
	char *qrySeq = (char *) malloc(totalQrySeqLen);

	/* Allocate GPU memory for results. */
	char *refIdxArr_d;
	cudaMalloc(&refIdxArr_d, numQrySeqPerIter * maxNumHits * sizeof(char));
	PRINT_CUDA_ERROR()

	int *shiftArr_d;
	cudaMalloc(&shiftArr_d, numQrySeqPerIter * maxNumHits * sizeof(int));
	PRINT_CUDA_ERROR()

	int *posArr_d;
	cudaMalloc(&posArr_d, numQrySeqPerIter * maxNumHits * sizeof(int));
	PRINT_CUDA_ERROR()

	/* Allocate CPU memory for results. */
	char *refIdxArr = (char *) malloc(numQrySeqPerIter * maxNumHits
			* sizeof(char));
	int *shiftArr = (int *) malloc(numQrySeqPerIter * maxNumHits * sizeof(int));
	int *posArr = (int *) malloc(numQrySeqPerIter * maxNumHits * sizeof(int));

	/* Search queries. */
	FILE *matchFilePtr = fopen(matchFile, "w");
	int randNum, j, count;
	time_t startTime, endTime;
	for (i = 0; i < numIter; ++i)
	{
		time(&startTime);

		fprintf(stderr, "\nIteration %d of %d:\n", (i + 1), numIter);

		/* Check if this is last iteration. */
		if (i == numIter - 1)
		{
			numQrySeqPerIter = numQrySeqLastIter;
			numQrsPerIter = numQrySeqPerIter / 2;
		}

		/* Copy query sequences from CPU to GPU. */
		fprintf(stderr, "   Copying queries from CPU to GPU...");
		count = 0;
		for (j = 0; j < numQrsPerIter; ++j)
		{
			qry[j] = qryGetNext2(qryFilePtr, qry[j]);
			strcpy(qrySeq + count, qry[j]->seq);
			count += qryLen;
			strcpy(qrySeq + count, qry[j]->revCompSeq);
			count += qryLen;
		}
		cudaMemcpy(qrySeq_d, qrySeq, totalQrySeqLen, cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

//		cudaMemset(refIdxArr_d, 3, numQrySeqPerIter * maxNumHits * sizeof(char));
//		PRINT_CUDA_ERROR()
//		cudaMemset(shiftArr_d, 3, numQrySeqPerIter * maxNumHits * sizeof(int));
//		PRINT_CUDA_ERROR()
//		cudaMemset(posArr_d, 3, numQrySeqPerIter * maxNumHits * sizeof(int));
//		PRINT_CUDA_ERROR()

		/* Call kernel function. */
		fprintf(stderr, "   Calling kernel function...");
		randNum = rand();
		dim3 dimGrid, dimBlock;
		dimBlock.x = 1;
		dimBlock.y = 1;
		dimBlock.z = 1;
		dimGrid.x = numQrySeqPerIter;
		dimGrid.y = 1;
		refSearchQuery2_kernel<<<dimGrid, dimBlock>>>(keys_d, vals_d,
				numKeys, numVals, qrySeq_d, refIdxArr_d, shiftArr_d, posArr_d,
				qryLen, seedLen, maxNumHits, randNum);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Copy results from GPU to CPU. */
		fprintf(stderr, "   Copying results from GPU to CPU...");
		cudaMemcpy(refIdxArr, refIdxArr_d, (numQrySeqPerIter * maxNumHits
				* sizeof(char)), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, (numQrySeqPerIter * maxNumHits
				* sizeof(int)), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, (numQrySeqPerIter * maxNumHits
				* sizeof(int)), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Print matches to a file. */
		fprintf(stderr, "   Printing intermediate results...");
		printHeuristicMatches2(qry, refIdxArr, shiftArr, posArr,
				numQrySeqPerIter, maxNumHits, matchFilePtr, seedLen);
		fprintf(stderr, "done.\n");

		time(&endTime);
		fprintf(stderr, "Time = %.2lf secs\n", difftime(endTime, startTime));

//		break;
//		if (i >= 19)
//			break;
	}

	/* Free resources. */
	fclose(matchFilePtr);
	fclose(qryFilePtr);
	refMapDelete();
	refNameMapDelete();
	free(keys);
	free(vals);
	free(qrySeq);
	free(refIdxArr);
	free(shiftArr);
	free(posArr);
	cudaFree(keys_d);
	cudaFree(vals_d);
	cudaFree(qrySeq_d);
	cudaFree(refIdxArr_d);
	cudaFree(shiftArr_d);
	cudaFree(posArr_d);
	for (i = 0; i < numQrsPerIter; ++i)
		free(qry[i]);
	free(qry);
}


/**
 * Prints matches obtained from heuristics to the given file.
 *
 * @param	qry			Query object array.
 * @param 	refIdxArr 	Reference index array.
 * @param 	shiftArr	Shift array (shift = ref pos - query pos).
 * @param	posArr		Reference position array.
 * @param	numQrsPerIter	Number of queries per iteration.
 * @param	maxNumHits	Maximum number of hits.
 * @param 	filePtr 	Pointer to the file where output would be written.
 * The file must be already open for writing.
 * @param	seedLen		Seed length.
 */
static void printHeuristicMatches2(Query **qry, char *refIdxArr,
		int *shiftArr, int *posArr, long numQrsPerIter, int maxNumHits,
		FILE *filePtr, int seedLen)
{
	static int j, qryIdx, arrIdx, pos;
	static long i;
	static uint localPos, globalPos;
	static int refIdxShiftBits = (sizeof(int) * NUM_BITS_IN_A_BYTE)
			- REF_IDX_BITS;
	for (i = 0; i < numQrsPerIter; i = i + 2)
	{
		qryIdx = i / 2;
		for (j = 0; j < maxNumHits; ++j)
		{
			arrIdx = (i * maxNumHits) + j;
			if (refIdxArr[arrIdx] != -1)
			{
				pos = posArr[arrIdx];
				pos = pos - (pos % seedLen);
				if ((pos - seedLen) >= 0)
					pos -= seedLen;
				localPos = refIdxArr[arrIdx] << refIdxShiftBits;
				localPos = localPos | (pos / seedLen);
				globalPos = refMapGetGlobalPos(localPos);
				fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\t%u\t%u\n",
						qry[qryIdx]->name, 0, refIdxArr[arrIdx],
						shiftArr[arrIdx], pos,
						qry[qryIdx]->nameOffset, qry[qryIdx]->seqOffset,
						refNameMapGetNameOffset(globalPos), globalPos);
			}
		}

		for (j = 0; j < maxNumHits; ++j)
		{
			arrIdx = ((i + 1) * maxNumHits) + j;
			if (refIdxArr[arrIdx] != -1)
			{
				pos = posArr[arrIdx];
				pos = pos - (pos % seedLen);
				if ((pos - seedLen) >= 0)
					pos -= seedLen;
				localPos = refIdxArr[arrIdx] << refIdxShiftBits;
				localPos = localPos | (pos / seedLen);
				globalPos = refMapGetGlobalPos(localPos);
				fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\t%u\t%u\n",
						qry[qryIdx]->name, 1, refIdxArr[arrIdx],
						shiftArr[arrIdx], pos,
						qry[qryIdx]->nameOffset, qry[qryIdx]->seqOffset,
						refNameMapGetNameOffset(globalPos), globalPos);
			}
		}
	}
}


/**
 * Prints matches obtained from heuristics to the given file.
 *
 * @param	id			Query ID.
 * @param	qry			Query for which the matches will be printed.
 * @param	isRevComp	A value of '1' indicates that the matches
 * are for reverse complement; a value of '0' indicates otherwise.
 * @param	mateIdx		This value is 1 if the query is the first
 * query of the pair, and 2 if it is the second query of the pair.
 * @param	matches		Best-matching reference coordinates obtained
 * from heuristic method.
 * @param	filePtr		Pointer to the file where output would be
 * written. The file must be already open for writing.
 * @param	seedLen		Seed length.
 */
static void printHeuristicMatches3(int id, Query *qry, int isRevComp,
		int mateIdx, HitList *matches, FILE *filePtr, int seedLen,
		int *numMatches)
{
	while (matches != NULL)
	{
		++(*numMatches);
		fprintf(filePtr, "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%u\t%u\n", id,
				qry->name, isRevComp, mateIdx, matches->index, matches->shift,
				matches->offset, qry->nameOffset, qry->seqOffset);
		matches = matches->next;
	}
}


/**
 * This is a wrapper function that wraps @a getNumIters function. This
 * function has been added so that @a getNumIters can be unit-tested.
 *
 * @param[out]	numIter			Number of iterations required to match
 * all queries. This value will be set to -1 if there is no memory
 * available for storing query sequences and results on the GPU global memory.
 * @param[out]	numQrsPerIter	Number of queries that the GPU can
 * accomodate per iteration. This value will be set to -1 if there is no memory
 * available for storing query sequences and results on the GPU global memory.
 * @param[out]	numQrsLastIter	Number of queries left for the last iteration.
 * This value will be set to -1 if there is no memory available for storing
 * query sequences and results on the GPU global memory.
 * @param		numQrs			Total number of queries.
 * @param		maxNumHits		Maximum number of hits per query.
 * @param		gpuTotalMem		Total GPU global memory.
 * @param		gpuMemMargin	Amount of memory that should not be included
 * in the calculation of number of iterations, number of queries in each
 * iteration, and number of queries in the last iteration.
 * @param		keysTotalMem	Total memory used by @a keys data structure.
 * @param		valsTotalMem	Total memory used by @a vals data structure.
 * @param		qrySeqMem		Memory required for each query sequence.
 * @param		refIdxMem		Memory required by each reference index.
 * @param		shiftMem		Memory required by each shift value.
 * @param		posMem			Memory required by each position value.
 */
void getNumIters_wrap(int *numIter, long *numQrsPerIter, long *numQrsLastIter,
		int numQrs, int maxNumHits, long gpuTotalMem, int gpuMemMargin,
		int keysTotalMem, int valsTotalMem, int qrySeqMem, int refIdxMem,
		int shiftMem, int posMem)
{
	getNumIters(numIter, numQrsPerIter, numQrsLastIter, numQrs, maxNumHits,
			gpuTotalMem, gpuMemMargin, keysTotalMem, valsTotalMem, qrySeqMem,
			refIdxMem, shiftMem, posMem);
}


/**
 * Calculates the number of queries that the GPU memory can accomodate
 * per iteration. It also calculates the number of queries left for the
 * last iteration. The outbound parameters will be set to -1 if there
 * is no memory available on the GPU global memory to store query sequences
 * and results.
 *
 * @param[out]	numIter			Number of iterations required to match
 * all queries. This value will be set to -1 if there is no memory
 * available for storing query sequences and results on the GPU global memory.
 * @param[out]	numQrsPerIter	Number of queries that the GPU can
 * accomodate per iteration. This value will be set to -1 if there is no memory
 * available for storing query sequences and results on the GPU global memory.
 * @param[out]	numQrsLastIter	Number of queries left for the last iteration.
 * This value will be set to -1 if there is no memory available for storing
 * query sequences and results on the GPU global memory.
 * @param		numQrs			Total number of queries.
 * @param		maxNumHits		Maximum number of hits per query.
 * @param		gpuTotalMem		Total GPU global memory.
 * @param		gpuMemMargin	Amount of memory that should not be included
 * in the calculation of number of iterations, number of queries in each
 * iteration, and number of queries in the last iteration.
 * @param		keysTotalMem	Total memory used by @a keys data structure.
 * @param		valsTotalMem	Total memory used by @a vals data structure.
 * @param		qrySeqMem		Memory required for each query sequence.
 * @param		refIdxMem		Memory required by each reference index.
 * @param		shiftMem		Memory required by each shift value.
 * @param		posMem			Memory required by each position value.
 */
static void getNumIters(int *numIter, long *numQrsPerIter, long *numQrsLastIter,
		int numQrs, int maxNumHits, long gpuTotalMem, int gpuMemMargin,
		int keysTotalMem, int valsTotalMem, int qrySeqMem, int refIdxMem,
		int shiftMem, int posMem)
{
	long remainingMem = gpuTotalMem - keysTotalMem - valsTotalMem - gpuMemMargin;

	/* If the remaining memory is less than or equal to 0, then set the
	 * outbound parameters to -1. */
	if (remainingMem <= 0)
	{
		*numIter = -1;
		*numQrsPerIter = -1;
		*numQrsLastIter = -1;
		return;
	}

	/* Calculate number of queries per iteration. */
	int resultMem = refIdxMem + shiftMem + posMem;
	*numQrsPerIter = remainingMem / (qrySeqMem + (maxNumHits * resultMem));
	if (*numQrsPerIter > numQrs)
		*numQrsPerIter = numQrs;

	/* Number of queries should not be greater than the total number of
	 * blocks. */
	int device;
	cudaGetDevice(&device);
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, device);
	long maxNumBlocksX = deviceProp.maxGridSize[0];
	PRINT_CUDA_ERROR()

	if (*numQrsPerIter > maxNumBlocksX)
		*numQrsPerIter = maxNumBlocksX;

	/* @a numQrsPerIter needs to be an even number because each query has
	 * 2 sequences: forward strand and reverse strand. */
	if ((*numQrsPerIter % 2) != 0)
		--(*numQrsPerIter);

	/* Calculate number of iterations and number of queries in last iteration. */
	*numIter = (int) ceil((float) numQrs / (*numQrsPerIter));
	int numQrsSecondLastIter = (*numIter - 1) * (*numQrsPerIter);
	*numQrsLastIter = numQrs - numQrsSecondLastIter;
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param qryFile1 Query file name.
 * @param qryFile2 Paired query file name.
 * @param numQrs Number of queries in the query file.
 * @param qryLen Length of each query in the query file.
 * @param refFile Reference file name.
 * @param numRefs Number of references in the reference file.
 * @param refLen Length of each reference sequence in the reference file.
 * @param seedLen Seed length.
 * @param matchFile File in which filtered reference sequence fragment
 * positions are stored.
 * @param maxNumHits Max number of hits per query.
 * @param minFragSize Minimum DNA fragment size used in sequencing.
 * @param maxFragSize Maximum DNA fragment size used in sequencing.
 */
void searchQueries_paired(char *qryFile1, char *qryFile2, int numQrs,
		int qryLen, char *refFile, int numRefs, int refLen, int seedLen,
		char *matchFile, uint maxNumHits, uint minFragSize, uint maxFragSize)
{
	if (qryFile1 == NULL || qryFile2 == NULL || numQrs <= 0 || qryLen == 0
			|| refFile == NULL || numRefs <= 0 || refLen <= 0 || seedLen <= 0
			|| matchFile == NULL)
		return;

	fprintf(stderr, "**************Phase I******************\n");
	fprintf(stderr, "Filtering reference sequences...\n");

	/* Preprocess reference. */
	fprintf(stderr, "   Preprocessing reference sequences...\n");
	RefPosList **hashTable = refPreprocess(refFile, seedLen);

	/* Preprocess query. */
	fprintf(stderr, "   Preprocessing query sequences...\n");
	int i;
	Query *qry1 = qryCreate2();
	Query *qry2 = qryCreate2();
	FILE *matchFilePtr = fopen(matchFile, "w");
	FILE *qryFilePtr1 = fopen(qryFile1, "r");
	FILE *qryFilePtr2 = fopen(qryFile2, "r");
	HitList **result = (HitList **) malloc(NUM_RESULTS * sizeof(HitList *));


	/* Find matching reference positions for each query pair and their
	 * reverse complements. */
	fprintf(stderr, "   Searching for queries in references...\n");
	for (i = 0; i < numQrs; ++i)
	{
		qry1 = qryGetNext2(qryFilePtr1, qry1);
		qry2 = qryGetNext2(qryFilePtr2, qry2);
		refSearchQuery_paired(hashTable, qry1, qry2, qryLen, seedLen,
				maxNumHits, minFragSize, maxFragSize, result);

		printHeuristicMatches_paired(qry1, qry2, result, matchFilePtr);

		if ((i + 1) % 100000 == 0)
		{
			fprintf(stderr, "   Searched %d queries...\n", i + 1);
			break;
		}
	}

	/* Clean up resources */
	fclose(matchFilePtr);
	fclose(qryFilePtr1);
	fclose(qryFilePtr2);
	qryDelete(qry1);
	qryDelete(qry2);
	hitListDelete(result[0]);
	hitListDelete(result[1]);
	hitListDelete(result[2]);
	hitListDelete(result[3]);
	hitListDelete(result[4]);
	hitListDelete(result[5]);
	hitListDelete(result[6]);
	hitListDelete(result[7]);
	free(result);

	/* Free memory occupied by hashTable. */
	if (hashTable != NULL)
	{
		int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE, seedLen);
		for (i = 0; i < numRefTuples; ++i)
			refPosListDelete(hashTable[i]);
		free(hashTable);
	}
}


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param qryFile1 Query file name.
 * @param qryFile2 Paired query file name.
 * @param numQrs Number of queries in the query file.
 * @param qryLen Length of each query in the query file.
 * @param refFile Reference file name.
 * @param seedLen Seed length.
 * @param matchFile File in which filtered reference sequence fragment
 * positions are stored.
 * @param maxNumHits Max number of hits per query.
 * @param[out] numMatches Total number of hits.
 */
void searchQueries_paired2(char *qryFile1, char *qryFile2, int numQrs,
		int qryLen, char *refFile, int seedLen, char *matchFile,
		uint maxNumHits, int *numMatches)
{
	if (qryFile1 == NULL || qryFile2 == NULL || qryLen == 0
			|| refFile == NULL || seedLen <= 0 || matchFile == NULL
			|| numQrs <= 0)
		return;

	/* Preprocess reference. */
	fprintf(stderr, "Preprocessing reference sequences...\n");
//	lookupTableCreate(refFile, seedLen, maxNumHits);
	lookupTable2Create(refFile, seedLen, maxNumHits, TUPLE_IGNORE_THRES);

	/* Preprocess query. */
	fprintf(stderr, "Preprocessing query sequences...");
	int i;
	Query *qry1 = qryCreate2();
	Query *qry2 = qryCreate2();
	FILE *matchFilePtr = fopen(matchFile, "w");
	FILE *qryFilePtr1 = fopen(qryFile1, "r");
	FILE *qryFilePtr2 = fopen(qryFile2, "r");
	fprintf(stderr, "done.\n");


	/* Find matching reference positions for each query pair and their
	 * reverse complements. */
	HitList *bestMatches = NULL;
	int id = 0;
	time_t startTime, endTime;
	double diffTime = 0.0;
	*numMatches = 0;
	fprintf(stderr, "Mapping queries to references...\n");
	for (i = 0; i < numQrs; ++i)
	{
		time(&startTime);
		qry1 = qryGetNext2(qryFilePtr1, qry1);
		qry2 = qryGetNext2(qryFilePtr2, qry2);
		++id;

		/* Find matching reference positions for query 1. */
//		bestMatches = lookupTableMapQry(qry1, qryLen, 0);
		bestMatches = lookupTable2MapQry(qry1, qryLen, 0);
		if (bestMatches != NULL)
		{
			printHeuristicMatches3(id, qry1, 0, 1, bestMatches, matchFilePtr,
					seedLen, numMatches);
			hitListDelete(bestMatches);
		}

		/* Find matching reference position for reverse complement of query 1. */
//		bestMatches = lookupTableMapQry(qry1, qryLen, 1);
		bestMatches = lookupTable2MapQry(qry1, qryLen, 1);
		if (bestMatches != NULL)
		{
			printHeuristicMatches3(id, qry1, 1, 1, bestMatches, matchFilePtr,
					seedLen, numMatches);
			hitListDelete(bestMatches);
		}

		/* Find matching reference positions for query 2. */
//		bestMatches = lookupTableMapQry(qry2, qryLen, 0);
		bestMatches = lookupTable2MapQry(qry2, qryLen, 0);
		if (bestMatches != NULL)
		{
			printHeuristicMatches3(id, qry2, 0, 2, bestMatches, matchFilePtr,
					seedLen, numMatches);
			hitListDelete(bestMatches);
		}

		/* Find matching reference position for reverse complement of query 2. */
//		bestMatches = lookupTableMapQry(qry2, qryLen, 1);
		bestMatches = lookupTable2MapQry(qry2, qryLen, 1);
		if (bestMatches != NULL)
		{
			printHeuristicMatches3(id, qry2, 1, 2, bestMatches, matchFilePtr,
					seedLen, numMatches);
			hitListDelete(bestMatches);
		}

		time(&endTime);
		diffTime += difftime(endTime, startTime);

		if ((i + 1) % 100000 == 0)
		{
			fprintf(stderr, "   Searched %d queries...\n", i + 1);
			fprintf(stderr, "      Time = %.2lf secs\n", diffTime);
			diffTime = 0.0;
//			break;
		}
	}

	/* Clean up resources */
	fclose(matchFilePtr);
	fclose(qryFilePtr1);
	fclose(qryFilePtr2);
	qryDelete(qry1);
	qryDelete(qry2);

	/* Free memory occupied by lookup table. */
	fprintf(stderr, "Deleting lookup table...");
//	lookupTableDelete();
	lookupTable2Delete();
	fprintf(stderr, "done.\n");
}


/**
 * Prints matches obtained from heuristics to the given file.
 *
 * @param qry Query for which the matches will be printed.
 * @param isReverseComplement A value of '1' indicates that the matches
 * are for reverse complement; a value of '0' indicates otherwise.
 * @param matches Best-matching reference coordinates obtained from heuristic
 * method.
 * @param filePtr Pointer to the file where output would be written. The file
 * must be already open for writing.
 * @param seedLen Seed length.
 */
static void printHeuristicMatches(Query *qry, int isReverseComplement,
		HitList *matches, FILE *filePtr, int seedLen, int *numMatches)
{
	while (matches != NULL)
	{
		++(*numMatches);
		fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\n", qry->name,
				isReverseComplement, matches->index, matches->shift,
				matches->offset, qry->nameOffset, qry->seqOffset);
		matches = matches->next;
	}
}


/**
 * Prints matches obtained from heuristics to the given file.
 *
 * @param 	qry			Query for which the matches will be printed.
 * @param 	isRevComp	A value of '1' indicates that the matches
 * are for reverse complement; a value of '0' indicates otherwise.
 * @param	refIdx		Reference indexes.
 * @param	shift		Shifts.
 * @param	refPos		Reference positions.
 * @param	numMatches	Number of matches in @a refIdx.
 * @param 	filePtr		Pointer to the file where output would be written. The file
 * must be already open for writing.
 */
static void printHeuristicMatches4(Query *qry, int isRevComp, char *refIdx,
		int *shift, int *refPos, int numMatches, FILE *filePtr)
{
	static int i;
	for (i = 0; i < numMatches; ++i)
	{
		fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\n", qry->name,
				isRevComp, refIdx[i], shift[i], refPos[i], qry->nameOffset,
				qry->seqOffset);
	}
}


/**
 * Prints matches obtained from heuristics to the given file.
 *
 * @param 		qry			Query objects.
 * @param 		refIdx		Reference indices.
 * @param 		refPos		Reference positions.
 * @param 		numResults	Number of elements in @a qry.
 * @param[out]	numMatches	Number of valid matches.
 * @param 		filePtr		Pointer to the output file.
 */
static void printHeuristicMatches5(Query **qry, char *refIdx, int *refPos,
		int numResults, int *numMatches, FILE *filePtr)
{
	static int i, j;
	j = 0;
	for (i = 0; i < numResults; ++i)
	{
		if (refIdx[i] != CHAR_MAX)
		{
			fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\n", qry[j]->name,
					FORWARD_STRAND, refIdx[i], 0, refPos[i], qry[j]->nameOffset,
					qry[j]->seqOffset);
			++(*numMatches);
		}
		++i;

		if (refIdx[i] != CHAR_MAX)
		{
			fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\n", qry[j]->name,
					REV_STRAND, refIdx[i], 0, refPos[i], qry[j]->nameOffset,
					qry[j]->seqOffset);
			++(*numMatches);
		}
		++j;
	}
}


/**
 * Prints matches obtained from heuristics to the given file.
 *
 * @param 		qry			Query objects.
 * @param 		refIdx		Reference indices.
 * @param 		refPos		Reference positions.
 * @param 		numQrs		Total number of forward and reverse query sequences.
 * @param		maxHits		Max allowed hits per query sequence.
 * @param[out]	numMatches	Number of valid matches.
 * @param 		filePtr		Pointer to the output file.
 */
static void printHeuristicMatches5_2(Query **qry, char *refIdx, int *refPos,
		int numQrs, int maxHits, int *numMatches, FILE *filePtr)
{
	static int i, j, k, l;
	j = 0;
	for (i = 0; i < numQrs; ++i)
	{
		for (k = 0; k < maxHits; ++k)
		{
			l = (i * maxHits) + k;
			if (refIdx[l] != -1)
			{
				fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\n", qry[j]->name,
						FORWARD_STRAND, refIdx[l], 0, refPos[l],
						qry[j]->nameOffset, qry[j]->seqOffset);
				++(*numMatches);
			}
		}
		++i;

		for (k = 0; k < maxHits; ++k)
		{
			l = (i * maxHits) + k;
			if (refIdx[l] != -1)
			{
				fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\n", qry[j]->name,
						REV_STRAND, refIdx[l], 0, refPos[l], qry[j]->nameOffset,
						qry[j]->seqOffset);
				++(*numMatches);
			}
		}
		++j;
	}
}


/**
 * Prints matches obtained from heuristics to the given file.
 *
 * @param	qry			Query objects.
 * @param	strand		Indicates whether the hits belong to a forward
 * or reverse strand.
 * @param 	refIdx		Reference indices.
 * @param 	refPos		Reference positions.
 * @param 	numHits		Number of hits for the give query.
 * @param 	filePtr		Pointer to the output file.
 */
static void printHeuristicMatches6(Query *qry, int strand, char *refIdx,
		int *refPos, int numHits, FILE *filePtr)
{
	static int i;
	static char *qrySeq;

	if (strand == 0)
		qrySeq = qry->seq;
	else
		qrySeq = qry->revCompSeq;
	for (i = 0; i < numHits; ++i)
	{
		fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\n", qrySeq,
				strand, refIdx[i], 0, refPos[i], qry->nameOffset,
				qry->seqOffset);
	}
}


/**
 * Prints matches obtained from heuristics to the given file.
 *
 * @param	qry			Query objects.
 * @param	strand		Indicates whether the hits belong to a forward
 * or reverse strand.
 * @param 	refIdx		Reference indices.
 * @param 	refPos		Reference positions.
 * @param 	numHits		Number of hits for the given query.
 * @param[out]	numMatches	Keeps track of the total number of hits.
 * @param 	filePtr		Pointer to the output file.
 */
static void printHeuristicMatches10(Query *qry, int strand, char *refIdx,
		int *refPos, int *refSize, int numHits, int *numMatches, FILE *filePtr)
{
	static int i, pos, clustEnd;
	for (i = 0; i < numHits; ++i)
	{
		pos = refPos[i];
		clustEnd = min(refSize[i], pos + CLUST_SPAN);
		while (pos < clustEnd)
		{
			fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\n", qry->name,
					strand, refIdx[i], 0, pos, qry->nameOffset, qry->seqOffset);
			++(*numMatches);
			pos += 100;
		}
	}
}


/**
 * Prints matches obtained from heuristics to the given file.
 *
 * @param	qry			Query objects.
 * @param	strand		Indicates whether the hits belong to a forward
 * or reverse strand.
 * @param 	refIdx		Reference indices.
 * @param 	refPos		Reference positions.
 * @param 	numHits		Number of hits for the given query.
 * @param[out]	numMatches	Keeps track of the total number of hits.
 * @param 	filePtr		Pointer to the output file.
 */
static void printHeuristicMatches11(Query *qry, int strand, char *refIdx,
		int *refPos, int numHits, int *numMatches, FILE *filePtr)
{
	static int i;
	for (i = 0; i < numHits; ++i)
	{
		fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\n", qry->name,
				strand, refIdx[i], 0, refPos[i], qry->nameOffset, qry->seqOffset);
		++(*numMatches);
	}
}


/**
 * Prints heuristic matches obtained from paired-end queries.
 *
 * @param qry1 First @a Query object.
 * @param qry2 Second @a Query object.
 * @param matches Array of @a HitList containing 8 elements. Element 1 contains
 * a list of hits for query 1; Element 2 contains the corresponding paired list
 * of hits for query 2. Element 3 contains a list of hits for reverse complement
 * of query 1; Element 4 contains the corresponding paired list of hits for
 * query 2. Element 5 contains a list of hits for query 1; Element 6 contains
 * the corresponding paired list of hits for the reverse complement of query 2.
 * Element 6 contains a list of hits for reverse complement of query 1;
 * Element 7 contains the corresponding paired list of hits for query 2.
 * @param filePtr Output file pointer.
 */
static void printHeuristicMatches_paired(const Query *qry1, const Query *qry2,
		HitList **matches, FILE *filePtr)
{
	if (matches == NULL)
		return;

	HitList *qry1Hit, *qry2Hit;
	int isQry1RevComp, isQry2RevComp;

	/* Query 1 and query 2. */
	if (matches[0] != NULL && matches[1] != NULL)
	{
		qry1Hit = matches[0];
		qry2Hit = matches[1];
		isQry1RevComp = 0;
		isQry2RevComp = 0;
		while (qry1Hit != NULL && qry2Hit != NULL)
		{
			fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\n"
				"%s\t%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\n", qry1->name,
					isQry1RevComp, qry1Hit->index, qry1Hit->shift,
					qry1Hit->offset, qry1->nameOffset, qry1->seqOffset,
					qry2->name, isQry1RevComp, qry2Hit->index,
					qry2Hit->shift, qry2Hit->offset, qry2->nameOffset,
					qry2->seqOffset);
			qry1Hit = qry1Hit->next;
			qry2Hit = qry2Hit->next;
		}
	}

	/* Query 1 reverse complement and query 2. */
	if (matches[2] != NULL && matches[3] != NULL)
	{
		qry1Hit = matches[2];
		qry2Hit = matches[3];
		isQry1RevComp = 1;
		isQry2RevComp = 0;
		while (qry1Hit != NULL && qry2Hit != NULL)
		{
			fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\n"
				"%s\t%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\n", qry1->name,
					isQry1RevComp, qry1Hit->index, qry1Hit->shift,
					qry1Hit->offset,qry1->nameOffset, qry1->seqOffset,
					qry2->name, isQry2RevComp, qry2Hit->index, qry2Hit->shift,
					qry2Hit->offset, qry2->nameOffset, qry2->seqOffset);
			qry1Hit = qry1Hit->next;
			qry2Hit = qry2Hit->next;
		}
	}

	/* Query 1 and reverse complement of query 2. */
	if (matches[4] != NULL && matches[5] != NULL)
	{
		qry1Hit = matches[4];
		qry2Hit = matches[5];
		isQry1RevComp = 0;
		isQry2RevComp = 1;
		while (qry1Hit != NULL && qry2Hit != NULL)
		{
			fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\n"
				"%s\t%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\n", qry1->name,
					isQry1RevComp, qry1Hit->index, qry1Hit->shift,
					qry1Hit->offset, qry1->nameOffset, qry1->seqOffset,
					qry2->name, isQry2RevComp, qry2Hit->index, qry2Hit->shift,
					qry2Hit->offset, qry2->nameOffset, qry2->seqOffset);
			qry1Hit = qry1Hit->next;
			qry2Hit = qry2Hit->next;
		}
	}

	/* Query 1 reverse complement and query 2 reverse complement. */
	if (matches[6] != NULL && matches[7] != NULL)
	{
		qry1Hit = matches[6];
		qry2Hit = matches[7];
		isQry1RevComp = 1;
		isQry2RevComp = 1;
		while (qry1Hit != NULL && qry2Hit != NULL)
		{
			fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\n"
				"%s\t%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\n", qry1->name,
					isQry1RevComp, qry1Hit->index, qry1Hit->shift,
					qry1Hit->offset, qry1->nameOffset, qry1->seqOffset,
					qry2->name, isQry2RevComp, qry2Hit->index, qry2Hit->shift,
					qry2Hit->offset, qry2->nameOffset, qry2->seqOffset);
			qry1Hit = qry1Hit->next;
			qry2Hit = qry2Hit->next;
		}
	}
}


/**
 * This is a wrapper function that wraps @a getNumBlocks function. It has
 * been created so that @a getNumBlocks can be unit-tested.
 *
 * @param[out]	blocksX			Holds the number of blocks in the X dimension.
 * @param[out]	blocksY			Holds the number of blocks in the Y dimension.
 * @param		availBlocksX	Number of available blocks in the X dimension.
 * @param		availBlocksY	Number if available blocks in the Y dimension.
 * @param		numQrsPerIter	Number of queries per iteration.
 */
void getNumBlocks_wrap(int *blocksX, int *blocksY, long availBlocksX,
		long availBlocksY, long numQrsPerIter)
{
	getNumBlocks(blocksX, blocksY, availBlocksX, availBlocksY, numQrsPerIter);
}


/**
 * Calculates the appropriate number of blocks in the X and Y dimension.
 *
 * @param[out]	blocksX			Holds the number of blocks in the X dimension.
 * @param[out]	blocksY			Holds the number of blocks in the Y dimension.
 * @param		availBlocksX	Number of available blocks in the X dimension.
 * @param		availBlocksY	Number if available blocks in the Y dimension.
 * @param		numQrsPerIter	Number of queries per iteration.
 */
static void getNumBlocks(int *blocksX, int *blocksY, long availBlocksX,
		long availBlocksY, long numQrsPerIter)
{
	long x = availBlocksX;
	long y = availBlocksY;
	long product = x * y;

	if (product <= numQrsPerIter)
	{
		*blocksX = (int) x;
		*blocksY = (int) y;
		return;
	}

	while (product >= numQrsPerIter)
	{
		--x;
		--y;
		product = x * y;
	}
	*blocksX = (int) (x + 1);
	*blocksY = (int) (y + 1);
}
