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

/**
 * @file
 *
 * @section DESCRIPTION
 *
 * Implements GPU-based filtering using seed-and-extend.
 */

#include "search2.h"
#include "common.h"
#include "lookupTable13.h"

#define	FORWARD_STRAND	0
#define REV_STRAND		1


/**
 * Finds reference sequence fragments that contain query tuples and saves
 * in the given file the coordinates of those reference sequence fragments
 * that have a higher density of query tuples.
 *
 * @param qryFile			Query file name.
 * @param refFile			Reference file name.
 * @param matchFile			File in which filtered reference sequence fragment
 * positions are stored.
 * @param seedLen			Seed length.
 * @param maxNumHits		Max number of hits per query.
 * @param numMatches		Number of matches obtained.
 * @param tupleIgnoreThres	Tuple ignore threshold.
 */
void search2Queries(char *qryFile, char *refFile, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches, int tupleIgnoreThres)
{
#define MAX_CHAR_VAL	127
#define SHR_MEM_MARGIN	256

	/* Query the GPU for configuration. */
	cudaSetDevice(SELECTED_DEVICE_INDEX);
	PRINT_CUDA_ERROR()
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, SELECTED_DEVICE_INDEX);
	PRINT_CUDA_ERROR()

	/* Create lookup table. */
	fprintf(stderr, "Preprocessing reference sequences...\n");
	int numTotalTuples;
	lookupTable13Create(refFile, seedLen, tupleIgnoreThres, &numTotalTuples);
	maxNumHits = 3;
//	lookupTable13Create2("/home/pgupta/data/gpusw/keys1000.txt",
//			"/home/pgupta/data/gpusw/vals1000.txt",
//			"/home/pgupta/data/gpusw/refPosHashMap1000.txt", maxNumHits,
//			seedLen, tupleIgnoreThres, &numTotalTuples);
	fprintf(stderr, "done.\n");

	/* Copy lookup table to GPU. */
	fprintf(stderr, "Copying lookup table to GPU...\n");
	int *keys, *values, numKeys, numValues, *numTuplesPerHash, *refPosHashMap;
	int maxRefPosComposite;
	lookupTable13CpyHashTableToGPU(&keys, &numKeys, &values, &numValues,
			&numTuplesPerHash, &refPosHashMap, &maxRefPosComposite);
	fprintf(stderr, "done.\n");

	/* Find out how many queries can fit on the GPU. */
	int numQrsPerIter = getNumQrsPerIter(seedLen, numTotalTuples, maxNumHits,
			maxRefPosComposite);
//	int numQrsPerIter = 1;
	fprintf(stderr, "Number of queries per iteration = %d\n", numQrsPerIter);

	/* Allocate memory on the CPU. */
	int i;
	int halfNumQrsPerIter = numQrsPerIter / 2;
	if (numQrsPerIter == 1)
		halfNumQrsPerIter = 1;
	Query **qry = (Query **) malloc(halfNumQrsPerIter * sizeof(Query *));
	for (i = 0; i < halfNumQrsPerIter; ++i)
		qry[i] = qryCreate2();
	int memPerQryInput1 = (MAX_QRY_SEQ_LENGTH + 1) * sizeof(char);
	int memPerQryInput2 = sizeof(uchar);
	char *qrs = (char *) malloc(numQrsPerIter * memPerQryInput1);
	uchar *qryLen = (uchar *) malloc(numQrsPerIter * sizeof(uchar));
	int *qryDistance = (int *) malloc(numQrsPerIter * sizeof(int));
	char *refIdx = (char *) malloc(maxNumHits * numQrsPerIter * sizeof(char));
	int *refPos = (int *) malloc(maxNumHits * numQrsPerIter * sizeof(int));

	/* Copy memory from host to constant memory on GPU. */
	lookupTable13CpyConstMemToGPU();

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
	int *qryDistance_d;
	cudaMalloc((void **) &qryDistance_d, numQrsPerIter * sizeof(int));
	PRINT_CUDA_ERROR()

	/* Find the ideal array size for the arrays in shared memory. */
	int availShrMem = deviceProp.sharedMemPerBlock - SHR_MEM_MARGIN;
	int currentShrMem = maxNumHits * sizeof(int);
	int idealArrSize = 0;
	int intSize = sizeof(int), charSize = sizeof(char);
	int sizePerElement = (5 * intSize) + charSize;
	do
	{
		currentShrMem += sizePerElement;
		++idealArrSize;
	} while (currentShrMem <= availShrMem);
	idealArrSize = idealArrSize - (idealArrSize % WARP_SIZE);
	currentShrMem = (idealArrSize * sizePerElement) + (maxNumHits * intSize);
	currentShrMem += (intSize - (currentShrMem % intSize));
	fprintf(stderr, "Array size on shared memory = %d\n", idealArrSize);
	fprintf(stderr, "Dynamically allocated shared memory = %d\n", currentShrMem);

	/* Filter queries on the GPU. */
	FILE *qryFilePtr = fopen(qryFile, "r");
	FILE *outputFilePtr = fopen(matchFile, "w");
	int isDone = FALSE, j = 0, totalQrs = 0, randNum;
	int numResultsPerIter, memResultsPerIter_char, memResultsPerIter_int;
	time_t startTime, endTime;
	double diffTime = 0.0;
	dim3 dimGrid, dimBlock;
	dimGrid.x = numQrsPerIter;
	dimGrid.y = 1;
	dimBlock.x = 64;
	dimBlock.y = 1;
	dimBlock.z = 1;
	*numMatches = 0;
	int k;
	while (isDone == FALSE)
	{
		time(&startTime);

		/* Fetch query sequences from query file. */
		fprintf(stderr, "   Fetching query sequences...");
		j = 0;
		for (i = 0; i < halfNumQrsPerIter; ++i)
		{
			if (qryGetNext2(qryFilePtr, qry[i]) == NULL)
			{
				isDone = TRUE;
				break;
			}
//			for (k = 0; k < 34; ++k)
//				qryGetNext2(qryFilePtr, qry[i]);
			strcpy(qrs + (j * MAX_QRY_SEQ_LENGTH), qry[i]->seq);
			qryLen[j] = qry[i]->length;
			++j;
			strcpy(qrs + (j * MAX_QRY_SEQ_LENGTH), qry[i]->revCompSeq);
			qryLen[j] = qry[i]->length;
			++j;
		}
		fprintf(stderr, "done.\n");
		numQrsPerIter = j;
		numResultsPerIter = maxNumHits * numQrsPerIter;
		memResultsPerIter_char = numResultsPerIter * charSize;
		memResultsPerIter_int = numResultsPerIter * intSize;

		/* Copy queries and their length to the GPU. */
		fprintf(stderr, "   Copying queries to GPU...");
		cudaMemcpy(qrs_d, qrs, numQrsPerIter * memPerQryInput1,
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qryLen_d, qryLen, numQrsPerIter * memPerQryInput2,
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		qryDistance[0] = 0;
		for (i = 1; i < numQrsPerIter; ++i)
			qryDistance[i] = qryDistance[i - 1] + qryLen[i - 1];
		cudaMemcpy(qryDistance_d, qryDistance, numQrsPerIter * intSize,
				cudaMemcpyHostToDevice);
		fprintf(stderr, "done.\n");

		/* Initialize output array on the GPU. */
		fprintf(stderr, "   Initializing output array on GPU...");
		cudaMemset(refIdx_d, -1, memResultsPerIter_char);
		PRINT_CUDA_ERROR()
		cudaMemset(refPos_d, -1, memResultsPerIter_int);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Call kernel function. */
		fprintf(stderr, "   Calling kernel function...");
		randNum = rand();
		lookupTable13MapQry2_gpu<<<dimGrid, dimBlock, currentShrMem>>>(keys,
				values, numTuplesPerHash, refPosHashMap, qrs_d, qryLen_d,
				qryDistance_d, MAX_QRY_SEQ_LENGTH, refIdx_d, refPos_d,
				maxNumHits, seedLen, randNum, idealArrSize);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Copy results from GPU to CPU. */
		fprintf(stderr, "   Copying results to CPU...");
		cudaMemcpy(refIdx, refIdx_d, memResultsPerIter_char,
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos, refPos_d, memResultsPerIter_int,
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		fprintf(stderr, "done.\n");

		/* Print results to file. */
		printHeuristicMatches(qry, refIdx, refPos, numQrsPerIter, maxNumHits,
				numMatches, outputFilePtr);

		/* Show progress. */
		time(&endTime);
		diffTime = difftime(endTime, startTime);
		totalQrs += numQrsPerIter;
		fprintf(stderr, "   Searched %d queries (total: %d)...%.2lf secs\n\n",
				numQrsPerIter, totalQrs, diffTime);
//		break;
	}

	/* Free memory. */
	fclose(qryFilePtr);
	fclose(outputFilePtr);
	free(qrs);
	free(qryLen);
	fprintf(stderr, "here0\n");
//	for (i = 0; i < halfNumQrsPerIter; ++i)
//		qryDelete(qry[i]);
	fprintf(stderr, "here1\n");
	free(qry);
	free(refIdx);
	free(refPos);
	cudaFree(qrs_d);
	cudaFree(qryLen_d);
	cudaFree(refIdx_d);
	cudaFree(refPos_d);
	cudaFree(keys);
	cudaFree(values);
	cudaFree(numTuplesPerHash);
	cudaFree(refPosHashMap);
}


/**
 * Returns the number of queries that can fit on the GPU within the available
 * amount of memory.
 *
 * @param	seedLen				Seed length.
 * @param	numTotalTuples		Total number of reference tuples.
 * @param	maxNumHits			Max number of allowed hits per query.
 * @param	maxRefPosComposite	Max composite reference position.
 * @return	The number of queries that can fit on the GPU.
 */
static int getNumQrsPerIter(int seedLen, int numTotalTuples, int maxNumHits,
		int maxRefPosComposite)
{
#define MIN_QRS_FOR_GPU 1000

	/* Query the GPU for configuration. */
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, SELECTED_DEVICE_INDEX);
	PRINT_CUDA_ERROR()

	/* Find how many queries can fit on the GPU. */
	long globalMem = deviceProp.totalGlobalMem;
	long memHashTableKeys = ((long) pow((float) DNA_ALPHABET_SIZE, seedLen))
			* sizeof(int);
	long memHashTableVals = numTotalTuples * sizeof(int);
	long memHashTableNumTuples = memHashTableKeys;
	long memRefPosHashMap = maxRefPosComposite * sizeof(int);
	long memHashTable = memHashTableKeys + memHashTableVals
			+ memHashTableNumTuples + memRefPosHashMap;
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

	return numQrsPerIter;
}


/**
 * Prints matches obtained from heuristics to the given file.
 *
 * @param 		qry				Query objects.
 * @param 		refIdx			Reference indices.
 * @param 		refPos			Reference positions.
 * @param 		numQrsPerIter	Number of elements in @a qry.
 * @param		maxNumHits		Max number of hits per query.
 * @param[out]	numMatches		Number of valid matches.
 * @param 		filePtr			Pointer to the output file.
 */
static void printHeuristicMatches(Query **qry, char *refIdx, int *refPos,
		int numQrsPerIter, int maxNumHits, int *numMatches, FILE *filePtr)
{
	static int i, j, k, l;
	for (i = 0; i < numQrsPerIter; ++i)
	{
		l = i / 2;
		if (i % 2 == 0)
		{
			for (j = 0; j < maxNumHits; ++j)
			{
				k = (i * maxNumHits) + j;
				if (refIdx[k] != -1)
				{
					fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\n",
							qry[l]->name, FORWARD_STRAND, refIdx[k], 0,
							refPos[k], qry[l]->nameOffset, qry[l]->seqOffset);
					++(*numMatches);
				}
			}
		}
		else
		{
			for (j = 0; j < maxNumHits; ++j)
			{
				k = (i * maxNumHits) + j;
				if (refIdx[k] != -1)
				{
					fprintf(filePtr, "%s\t%d\t%d\t%d\t%d\t%lu\t%lu\n",
							qry[l]->name, REV_STRAND, refIdx[k], 0, refPos[k],
							qry[l]->nameOffset, qry[l]->seqOffset);
					++(*numMatches);
				}
			}
		}
	}
}
