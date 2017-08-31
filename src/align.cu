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

#include <stdio.h>
#include "align.h"
#include "common.h"
#include "smithWaterman.h"
#include "smithWaterman2.h"
#include "memory.h"
#include "output.h"
#include "input.h"
#include "query.h"
#include "refMap.h"
#include <time.h>
#include <pthread.h>

#define NUM_SEQS			50
#define	_NUM_THREADS		8


/**
 * Aligns queries to references.
 *
 * @param queryFileName Path to the query file
 * @param querySeqLength Length of each query sequence
 * @param refFileName Path to the reference file
 * @param matchFile Temporary file in which file offsets of each query and
 * its best-matching reference sequences have been stored.
 * @param outFormat The output format in which the alignments would be printed.
 * @param match Match score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param outputFile Output file where alignments will be stored.
 * @param numMatches Total number of hits obtained in Phase 1.
 */
void alignQueries(char *queryFileName, int querySeqLength, char *refFileName,
		char *matchFile, int outFormat, float match, float mismatch,
		float gapOpenPenalty, float gapExtPenalty, const char *outputFile,
		int numMatches)
{
	fprintf(stderr, "Number of matches = %d\n", numMatches);
	int refSeqLength = 3 * querySeqLength;


	/*
	 * Get GPU specifications (max number of threads per block,
	 * max number of blocks per grid dimension, etc.) using CUDA API
	 */
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	PRINT_CUDA_ERROR()
	cudaDeviceProp deviceProp[deviceCount];
	PRINT_CUDA_ERROR()

	int device;
	for (device = 0; device < deviceCount; ++device)
		cudaGetDeviceProperties(&deviceProp[device], device);
	device = -1;
	cudaGetDevice(&device);
	if (device == -1)
	{
		device = SELECTED_DEVICE_INDEX;
		cudaSetDevice(device);
		PRINT_CUDA_ERROR()
	}
	PRINT_CUDA_ERROR()

	int maxRegsPerBlock = deviceProp[device].regsPerBlock;
	int maxThreadsPerBlock = (int) floor((float) maxRegsPerBlock / NUM_REGS_USED);
	if (maxThreadsPerBlock > deviceProp[device].maxThreadsPerBlock)
		maxThreadsPerBlock = deviceProp[device].maxThreadsPerBlock;


	int memPerRefQueryPair = getMemPerRefQryPair(refSeqLength, querySeqLength);
	fprintf(stderr, "Calculate max number of queries per batch...");
	int maxQueriesPerIteration = getNumQueriesPerIteration(deviceProp,
			memPerRefQueryPair, numMatches);
	fprintf(stderr, "(%d)...done.\n", maxQueriesPerIteration);
	int numQryIterations = (int) ceil((float) numMatches
			/ maxQueriesPerIteration);


	/* Allocate memory on host */
	fprintf(stderr, "Allocating memory on host...");
	int combinedLength = refSeqLength + querySeqLength + 2;
	float *score;
	char *refSeq;
	char *querySeq;
	int *alignStart;
	int *alignEnd;
	int *alignLength;
	int *refDistance;
	char *refConsensus;
	char *queryConsensus;
	char *queryName;
	char *refName;
	allocateHostMem(&score, &refDistance, &alignStart, &alignEnd, &alignLength,
			&queryConsensus, &refConsensus, &queryName, &querySeq, &refName,
			&refSeq, maxQueriesPerIteration, MAX_QRY_NAME_LENGTH,
			MAX_REF_NAME_LENGTH, querySeqLength, refSeqLength, combinedLength);
	fprintf(stderr, "done.\n");

	/* Allocate memory on device */
	fprintf(stderr, "Allocating memory on device...");
	char *refSeq_d;
	char *querySeq_d;
	float *score_d;
	float *hMatrix_d;
	int *cellBacktracker_d;
	int *alignStart_d;
	int *alignEnd_d;
	int *alignLength_d;
	char *refConsensus_d;
	char *queryConsensus_d;
	allocateDeviceMem(&querySeq_d, &refSeq_d, &score_d, &hMatrix_d,
			&cellBacktracker_d, &alignStart_d, &alignEnd_d, &alignLength_d,
			&queryConsensus_d, &refConsensus_d, maxQueriesPerIteration,
			querySeqLength, refSeqLength, combinedLength);
	fprintf(stderr, "done.\n");

	uint refNameOffset;
	uint refSeqOffset;
	uint refSeqFragOffset;
	uint queryNameOffset;
	uint querySeqOffset;
	uint isReverseComplement;
	FILE *tmpFilePtr = fopen(matchFile, "r");
	FILE *outputFilePtr = fopen(outputFile, "w");
	int i, j, refIdx;
	int refsPerIteration = 0, queriesPerIteration = 0;
	char tmpQrySeq[querySeqLength + 1];
	int numShiftBases = refSeqLength / 3;
	time_t startTime, endTime;
	double diffTime = 0.0;
	FILE *queryFilePtr = fopen(queryFileName, "r");
	FILE *refFilePtr = fopen(refFileName, "r");

//	fprintf(stderr, "Creating reference map...");
//	time(&startTime);
//	refMapCreate(refFileName);
//	time(&endTime);
//	diffTime = difftime(endTime, startTime);
//	fprintf(stderr, "done (time = %.2lf secs).\n", diffTime);
	fprintf(stderr, "Computing alignments:\n");
	for (i = 0; i < numQryIterations; ++i)
	{
		time(&startTime);

		memset(refSeq, 0, maxQueriesPerIteration * (refSeqLength + 1)
				* sizeof(char));

		/* Get number of queries for this iteration */
		queriesPerIteration = getNumSeqsInCurrentIteration(i, numMatches,
				numQryIterations, maxQueriesPerIteration);


		/* Fetch reference and query sequences */
		for (j = 0; j < queriesPerIteration; ++j)
		{
			getOffsets(tmpFilePtr, &queryNameOffset, &querySeqOffset,
					&refNameOffset, &refSeqOffset, &refSeqFragOffset,
					refDistance + j, &isReverseComplement, &refIdx);

			refDistance[j] = max(0, refDistance[j] - numShiftBases);
			refMapGetOffsets(refIdx, refDistance[j], &refNameOffset,
					&refSeqFragOffset);

//			getSeq(refFilePtr, refSeqFragOffset,
//					refSeq + (j * (refSeqLength + 1)), refSeqLength);
			getRefSeq2(refFilePtr, refSeqFragOffset,
					refSeq + (j * (refSeqLength + 1)), refSeqLength);
			getSeqName(refFilePtr, refNameOffset, refName + (j
					* MAX_SEQ_NAME_LENGTH));

			if (isReverseComplement == 1)
			{
				getSeq(queryFilePtr, querySeqOffset, tmpQrySeq, querySeqLength);
				qryGetReverseComplement(tmpQrySeq, querySeqLength,
						querySeq + (j * (querySeqLength + 1)));
				getSeqName(queryFilePtr, queryNameOffset, queryName + (j
									* MAX_SEQ_NAME_LENGTH));
				strcat(queryName + (j * MAX_SEQ_NAME_LENGTH), "*");
			}
			else
			{
				getSeq(queryFilePtr, querySeqOffset, querySeq
						+ (j * (querySeqLength + 1)), querySeqLength);
				getSeqName(queryFilePtr, queryNameOffset, queryName + (j
						* MAX_SEQ_NAME_LENGTH));
			}
		}

		/* Copy queries from host to device */
		copyMemHostToDevice(querySeq_d, querySeq, queriesPerIteration
				* (querySeqLength + 1) * sizeof(char));

		/* Copy references from host to device */
		copyMemHostToDevice(refSeq_d, refSeq, queriesPerIteration
				* (refSeqLength + 1) * sizeof(char));

		cudaMemset(score_d, 0, queriesPerIteration * sizeof(float));
		PRINT_CUDA_ERROR()

		/* Call device function */
		alignSequences(querySeq_d, refSeq_d,
				querySeqLength, refSeqLength, match, mismatch, gapOpenPenalty,
				gapExtPenalty, score_d, alignStart_d, alignEnd_d,
				alignLength_d, queryConsensus_d, refConsensus_d,
				queriesPerIteration, refsPerIteration, maxThreadsPerBlock);

		/* Copy device memory back to host memory */
		copyMemDeviceToHost(score, score_d, alignStart, alignStart_d, alignEnd,
				alignEnd_d, alignLength, alignLength_d, queryConsensus,
				queryConsensus_d, refConsensus, refConsensus_d,
				queriesPerIteration, refsPerIteration, combinedLength);

		/* Print results */
		if (outFormat == SAM_OUTPUT_FORMAT)
		{
			printOutputSam(queriesPerIteration, queryName, refName, score,
					alignStart, alignEnd, alignLength, combinedLength,
					queryConsensus, refConsensus, outputFilePtr);
		}
		else
		{
			chooseBestAlignment(queryName, MAX_SEQ_NAME_LENGTH, score,
					queriesPerIteration);
			printOutput(queriesPerIteration, queryName, refName, score,
					refDistance, alignStart, alignEnd, alignLength,
					combinedLength, queryConsensus, refConsensus,
					outputFilePtr);
		}

		time(&endTime);
		diffTime = difftime(endTime, startTime);
		fprintf(stderr, "   Completed %d / %d query iterations...%.2lf secs\n",
				(i + 1), numQryIterations, diffTime);
	}
	fclose(queryFilePtr);
	fclose(refFilePtr);
	fclose(tmpFilePtr);
	fclose(outputFilePtr);
//	refMapFree();

	/* Free resources on device */
	fprintf(stderr, "Freeing device resources...\n");
	cudaFree(refSeq_d);
	cudaFree(querySeq_d);
	cudaFree(score_d);
	cudaFree(alignStart_d);
	cudaFree(alignEnd_d);
	cudaFree(alignLength_d);
	cudaFree(refConsensus_d);
	cudaFree(queryConsensus_d);
	PRINT_CUDA_ERROR()

	/* Free resources on host */
	fprintf(stderr, "Freeing host resources...\n");
	free(queryName);
	free(querySeq);
	free(refName);
	free(refSeq);
	free(score);
	free(alignStart);
	free(alignEnd);
	free(alignLength);
	free(refConsensus);
	free(queryConsensus);
	free(refDistance);
}


/**
 * Aligns queries to references.
 *
 * @param queryFileName Path to the query file
 * @param querySeqLength Length of each query sequence
 * @param refFileName Path to the reference file
 * @param matchFile Temporary file in which file offsets of each query and
 * its best-matching reference sequences have been stored.
 * @param outFormat The output format in which the alignments would be printed.
 * @param match Match score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param outputFile Output file where alignments will be stored.
 * @param numMatches Total number of hits obtained in Phase 1.
 */
void alignQueries2(char *queryFileName, int querySeqLength, char *refFileName,
		char *matchFile, int outFormat, float match, float mismatch,
		float gapOpenPenalty, float gapExtPenalty, const char *outputFile,
		int numMatches)
{
	fprintf(stderr, "Number of matches = %d\n", numMatches);
	int refSeqLength = 3 * querySeqLength;


	/*
	 * Get GPU specifications (max number of threads per block,
	 * max number of blocks per grid dimension, etc.) using CUDA API
	 */
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	PRINT_CUDA_ERROR()
	cudaDeviceProp deviceProp[deviceCount];
	PRINT_CUDA_ERROR()

	int device;
	for (device = 0; device < deviceCount; ++device)
		cudaGetDeviceProperties(&deviceProp[device], device);
	device = -1;
	cudaGetDevice(&device);
	if (device == -1)
	{
		device = SELECTED_DEVICE_INDEX;
		cudaSetDevice(device);
		PRINT_CUDA_ERROR()
	}
	PRINT_CUDA_ERROR()

	int maxRegsPerBlock = deviceProp[device].regsPerBlock;
	int maxThreadsPerBlock = (int) floor((float) maxRegsPerBlock / NUM_REGS_USED);
	if (maxThreadsPerBlock > deviceProp[device].maxThreadsPerBlock)
		maxThreadsPerBlock = deviceProp[device].maxThreadsPerBlock;


	int memPerRefQueryPair = getMemPerRefQryPair(refSeqLength, querySeqLength);
	fprintf(stderr, "Calculate max number of queries per batch...");
//	int maxQueriesPerIteration = getNumQueriesPerIteration(deviceProp,
//			memPerRefQueryPair, numMatches);
	int maxQueriesPerIteration = numMatches;
	fprintf(stderr, "(%d)...done.\n", maxQueriesPerIteration);
	int numQryIterations = (int) ceil((float) numMatches
			/ maxQueriesPerIteration);


	/* Allocate memory on host */
	fprintf(stderr, "Allocating memory on host...");
	int combinedLength = refSeqLength + querySeqLength + 2;
	float *score;
	char *refSeq;
	char *querySeq;
	int *alignStart;
	int *alignEnd;
	int *alignLength;
	int *refDistance;
	char *refConsensus;
	char *queryConsensus;
	char *queryName;
	char *refName;
	allocateHostMem(&score, &refDistance, &alignStart, &alignEnd, &alignLength,
			&queryConsensus, &refConsensus, &queryName, &querySeq, &refName,
			&refSeq, maxQueriesPerIteration, MAX_QRY_NAME_LENGTH,
			MAX_REF_NAME_LENGTH, querySeqLength, refSeqLength, combinedLength);
	fprintf(stderr, "done.\n");

	/* Allocate memory on device */
	fprintf(stderr, "Allocating memory on device...");
	char *refSeq_d;
	char *querySeq_d;
	float *score_d;
	float *hMatrix_d;
	int *cellBacktracker_d;
	int *alignStart_d;
	int *alignEnd_d;
	int *alignLength_d;
	char *refConsensus_d;
	char *queryConsensus_d;
	allocateDeviceMem(&querySeq_d, &refSeq_d, &score_d, &hMatrix_d,
			&cellBacktracker_d, &alignStart_d, &alignEnd_d, &alignLength_d,
			&queryConsensus_d, &refConsensus_d, maxQueriesPerIteration,
			querySeqLength, refSeqLength, combinedLength);
	fprintf(stderr, "done.\n");

	uint refNameOffset;
	uint refSeqOffset;
	uint refSeqFragOffset;
	uint queryNameOffset;
	uint querySeqOffset;
	uint isReverseComplement;
	FILE *tmpFilePtr = fopen(matchFile, "r");
	FILE *outputFilePtr = fopen(outputFile, "w");
	int i, j, refIdx;
	int refsPerIteration = 0, queriesPerIteration = 0;
	char tmpQrySeq[querySeqLength + 1];
//	int numShiftBases = refSeqLength / 3;
	int numShiftBases = 0;
	time_t startTime, endTime;
	double diffTime = 0.0;
	FILE *queryFilePtr = fopen(queryFileName, "r");
	FILE *refFilePtr = fopen(refFileName, "r");

//	fprintf(stderr, "Creating reference map...");
//	time(&startTime);
//	refMapCreate(refFileName);
//	time(&endTime);
//	diffTime = difftime(endTime, startTime);
//	fprintf(stderr, "done (time = %.2lf secs).\n", diffTime);
	fprintf(stderr, "Computing alignments:\n");
	for (i = 0; i < numQryIterations; ++i)
	{
		time(&startTime);

		memset(refSeq, 0, maxQueriesPerIteration * (refSeqLength + 1)
				* sizeof(char));

		/* Get number of queries for this iteration */
		queriesPerIteration = getNumSeqsInCurrentIteration(i, numMatches,
				numQryIterations, maxQueriesPerIteration);


		/* Fetch reference and query sequences */
		fprintf(stderr, "   Fetching sequences...");
		for (j = 0; j < queriesPerIteration; ++j)
		{
			getOffsets(tmpFilePtr, &queryNameOffset, &querySeqOffset,
					&refNameOffset, &refSeqOffset, &refSeqFragOffset,
					refDistance + j, &isReverseComplement, &refIdx);

			refDistance[j] = max(0, refDistance[j] - numShiftBases);
			refMapGetOffsets(refIdx, refDistance[j], &refNameOffset,
					&refSeqFragOffset);

//			getSeq(refFilePtr, refSeqFragOffset,
//					refSeq + (j * (refSeqLength + 1)), refSeqLength);
//			getRefSeq(refFilePtr, refSeqFragOffset,
//					refSeq + (j * (refSeqLength + 1)), refSeqLength);
			getRefSeq2(refFilePtr, refSeqFragOffset,
					refSeq + (j * (refSeqLength + 1)), refSeqLength);
			getSeqName(refFilePtr, refNameOffset, refName + (j
					* MAX_SEQ_NAME_LENGTH));


			if (isReverseComplement == 1)
			{
				getSeq(queryFilePtr, querySeqOffset, tmpQrySeq, querySeqLength);
				qryGetReverseComplement(tmpQrySeq, querySeqLength,
						querySeq + (j * (querySeqLength + 1)));
				getSeqName(queryFilePtr, queryNameOffset, queryName + (j
									* MAX_SEQ_NAME_LENGTH));
				strcat(queryName + (j * MAX_SEQ_NAME_LENGTH), "*");
			}
			else
			{
				getSeq(queryFilePtr, querySeqOffset, querySeq
						+ (j * (querySeqLength + 1)), querySeqLength);
				getSeqName(queryFilePtr, queryNameOffset, queryName + (j
						* MAX_SEQ_NAME_LENGTH));
			}
		}
		fprintf(stderr, "done.\n");

		/* Copy queries from host to device */
		fprintf(stderr, "   Copying query sequences from host to device...");
		copyMemHostToDevice(querySeq_d, querySeq, queriesPerIteration
				* (querySeqLength + 1) * sizeof(char));
		fprintf(stderr, "done.\n");

		/* Copy references from host to device */
		fprintf(stderr, "   Copying ref sequences from host to device...");
		copyMemHostToDevice(refSeq_d, refSeq, queriesPerIteration
				* (refSeqLength + 1) * sizeof(char));
		fprintf(stderr, "done.\n");

		cudaMemset(score_d, 0, queriesPerIteration * sizeof(float));
		PRINT_CUDA_ERROR()

		/* Call device function */
		fprintf(stderr, "   Executing kernel function...");
		alignSequences(querySeq_d, refSeq_d,
				querySeqLength, refSeqLength, match, mismatch, gapOpenPenalty,
				gapExtPenalty, score_d, alignStart_d, alignEnd_d,
				alignLength_d, queryConsensus_d, refConsensus_d,
				queriesPerIteration, refsPerIteration, maxThreadsPerBlock);
		fprintf(stderr, "done.\n");

		/* Copy device memory back to host memory */
		fprintf(stderr, "   Copying results from device to host...");
		copyMemDeviceToHost(score, score_d, alignStart, alignStart_d, alignEnd,
				alignEnd_d, alignLength, alignLength_d, queryConsensus,
				queryConsensus_d, refConsensus, refConsensus_d,
				queriesPerIteration, refsPerIteration, combinedLength);
		fprintf(stderr, "done.\n");

		/* Print results */
		fprintf(stderr, "   Printing results to output file...");
		if (outFormat == SAM_OUTPUT_FORMAT)
		{
			printOutputSam(queriesPerIteration, queryName, refName, score,
					alignStart, alignEnd, alignLength, combinedLength,
					queryConsensus, refConsensus, outputFilePtr);
		}
		else
		{
			chooseBestAlignment(queryName, MAX_SEQ_NAME_LENGTH, score,
					queriesPerIteration);
			printOutput(queriesPerIteration, queryName, refName, score,
					refDistance, alignStart, alignEnd, alignLength,
					combinedLength, queryConsensus, refConsensus,
					outputFilePtr);
		}
		fprintf(stderr, "done.\n");

		time(&endTime);
		diffTime = difftime(endTime, startTime);
		fprintf(stderr, "   Completed %d / %d query iterations...%.2lf secs\n\n",
				(i + 1), numQryIterations, diffTime);
	}
	fclose(queryFilePtr);
	fclose(refFilePtr);
	fclose(tmpFilePtr);
	fclose(outputFilePtr);
//	refMapFree();

	/* Free resources on device */
	fprintf(stderr, "Freeing device resources...\n");
	cudaFree(refSeq_d);
	cudaFree(querySeq_d);
	cudaFree(score_d);
	cudaFree(alignStart_d);
	cudaFree(alignEnd_d);
	cudaFree(alignLength_d);
	cudaFree(refConsensus_d);
	cudaFree(queryConsensus_d);
	PRINT_CUDA_ERROR()

	/* Free resources on host */
	fprintf(stderr, "Freeing host resources...\n");
	free(queryName);
	free(querySeq);
	free(refName);
	free(refSeq);
	free(score);
	free(alignStart);
	free(alignEnd);
	free(alignLength);
	free(refConsensus);
	free(queryConsensus);
	free(refDistance);
}


/**
 * Aligns queries to references.
 *
 * @param queryFileName Path to the query file
 * @param maxQryLen Length of each query sequence
 * @param refFileName Path to the reference file
 * @param matchFile Temporary file in which file offsets of each query and
 * its best-matching reference sequences have been stored.
 * @param outFormat The output format in which the alignments would be printed.
 * @param match Match score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param outputFile Output file where alignments will be stored.
 * @param numMatches Total number of hits obtained in Phase 1.
 */
void alignQueries3(char *queryFileName, int maxQryLen, char *refFileName,
		char *matchFile, int outFormat, float match, float mismatch,
		float gapOpenPenalty, float gapExtPenalty, const char *outputFile,
		int numMatches)
{
	fprintf(stderr, "Number of matches = %d\n", numMatches);
	int refSeqLength = REF_LENGTH;

	/* Get GPU specifications (max number of threads per block,
	 * max number of blocks per grid dimension, etc.) using CUDA API. */
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	PRINT_CUDA_ERROR()
	cudaDeviceProp deviceProp[deviceCount];
	PRINT_CUDA_ERROR()

	int device;
	for (device = 0; device < deviceCount; ++device)
		cudaGetDeviceProperties(&deviceProp[device], device);
	device = -1;
	cudaGetDevice(&device);
	if (device == -1)
	{
		device = SELECTED_DEVICE_INDEX;
		cudaSetDevice(device);
		PRINT_CUDA_ERROR()
	}
	PRINT_CUDA_ERROR()

	int maxRegsPerBlock = deviceProp[device].regsPerBlock;
	int maxThreadsPerBlock = (int) floor((float) maxRegsPerBlock / NUM_REGS_USED);
	if (maxThreadsPerBlock > deviceProp[device].maxThreadsPerBlock)
		maxThreadsPerBlock = deviceProp[device].maxThreadsPerBlock;


	int memPerRefQueryPair = getMemPerRefQryPair(refSeqLength, maxQryLen);
	fprintf(stderr, "Calculate max number of queries per batch...");
	int maxQueriesPerIteration = getNumQueriesPerIteration(deviceProp,
			memPerRefQueryPair, numMatches);
	fprintf(stderr, "(%d)...done.\n", maxQueriesPerIteration);
	int numQryIterations = (int) ceil((float) numMatches
			/ maxQueriesPerIteration);


	/* Allocate memory on host */
	fprintf(stderr, "Allocating memory on host...");
	int combinedLength = refSeqLength + maxQryLen + 2;
	float *score;
	char *refSeq;
	char *querySeq;
	int *alignStart;
	int *alignEnd;
	int *alignLength;
	int *refDistance;
	char *refConsensus;
	char *queryConsensus;
	char *queryName;
	char *refName;
	allocateHostMem(&score, &refDistance, &alignStart, &alignEnd, &alignLength,
			&queryConsensus, &refConsensus, &queryName, &querySeq, &refName,
			&refSeq, maxQueriesPerIteration, MAX_QRY_NAME_LENGTH,
			MAX_REF_NAME_LENGTH, maxQryLen, refSeqLength, combinedLength);
	fprintf(stderr, "done.\n");

	/* Allocate memory on device */
	fprintf(stderr, "Allocating memory on device...");
	char *refSeq_d;
	char *querySeq_d;
	float *score_d;
	float *hMatrix_d;
	int *cellBacktracker_d;
	int *alignStart_d;
	int *alignEnd_d;
	int *alignLength_d;
	char *refConsensus_d;
	char *queryConsensus_d;
	allocateDeviceMem(&querySeq_d, &refSeq_d, &score_d, &hMatrix_d,
			&cellBacktracker_d, &alignStart_d, &alignEnd_d, &alignLength_d,
			&queryConsensus_d, &refConsensus_d, maxQueriesPerIteration,
			maxQryLen, refSeqLength, combinedLength);
	fprintf(stderr, "done.\n");

	uint refNameOffset;
	uint refSeqOffset;
	uint refSeqFragOffset;
	uint queryNameOffset;
	uint querySeqOffset;
	uint isReverseComplement;
	FILE *tmpFilePtr = fopen(matchFile, "r");
	FILE *outputFilePtr = fopen(outputFile, "w");
	int i, j, refIdx;
	int refsPerIteration = 0, queriesPerIteration = 0;
	char tmpQrySeq[maxQryLen + 1];
	int numShiftBases = refSeqLength / 3;
	time_t startTime, endTime;
	double diffTime = 0.0;
	FILE *queryFilePtr = fopen(queryFileName, "r");
	FILE *refFilePtr = fopen(refFileName, "r");

//	fprintf(stderr, "Creating reference map...");
//	time(&startTime);
//	refMapCreate(refFileName);
//	time(&endTime);
//	diffTime = difftime(endTime, startTime);
//	fprintf(stderr, "done (time = %.2lf secs).\n", diffTime);
	fprintf(stderr, "Computing alignments:\n");
	for (i = 0; i < numQryIterations; ++i)
	{
		time(&startTime);

		memset(refSeq, 0, maxQueriesPerIteration * (refSeqLength + 1)
				* sizeof(char));

		/* Get number of queries for this iteration */
		queriesPerIteration = getNumSeqsInCurrentIteration(i, numMatches,
				numQryIterations, maxQueriesPerIteration);


		/* Fetch reference and query sequences */
		for (j = 0; j < queriesPerIteration; ++j)
		{
			getOffsets(tmpFilePtr, &queryNameOffset, &querySeqOffset,
					&refNameOffset, &refSeqOffset, &refSeqFragOffset,
					refDistance + j, &isReverseComplement, &refIdx);

			refDistance[j] = max(0, refDistance[j] - numShiftBases);
			refMapGetOffsets(refIdx, refDistance[j], &refNameOffset,
					&refSeqFragOffset);

//			getSeq(refFilePtr, refSeqFragOffset,
//					refSeq + (j * (refSeqLength + 1)), refSeqLength);
			getRefSeq2(refFilePtr, refSeqFragOffset,
					refSeq + (j * (refSeqLength + 1)), refSeqLength);
			getSeqName(refFilePtr, refNameOffset, refName + (j
					* MAX_SEQ_NAME_LENGTH));

			if (isReverseComplement == 1)
			{
//				getSeq(queryFilePtr, querySeqOffset, tmpQrySeq, querySeqLength);
				getSeq2(queryFilePtr, querySeqOffset, tmpQrySeq, maxQryLen);
				qryGetReverseComplement(tmpQrySeq, maxQryLen,
						querySeq + (j * (maxQryLen + 1)));
				getSeqName(queryFilePtr, queryNameOffset, queryName + (j
									* MAX_SEQ_NAME_LENGTH));
				strcat(queryName + (j * MAX_SEQ_NAME_LENGTH), "*");
			}
			else
			{
//				getSeq(queryFilePtr, querySeqOffset, querySeq
//						+ (j * (querySeqLength + 1)), querySeqLength);
				getSeq2(queryFilePtr, querySeqOffset, querySeq
						+ (j * (maxQryLen + 1)), maxQryLen);
				getSeqName(queryFilePtr, queryNameOffset, queryName + (j
						* MAX_SEQ_NAME_LENGTH));
			}
		}

		/* Copy queries from host to device */
		copyMemHostToDevice(querySeq_d, querySeq, queriesPerIteration
				* (maxQryLen + 1) * sizeof(char));

		/* Copy references from host to device */
		copyMemHostToDevice(refSeq_d, refSeq, queriesPerIteration
				* (refSeqLength + 1) * sizeof(char));

		cudaMemset(score_d, 0, queriesPerIteration * sizeof(float));
		PRINT_CUDA_ERROR()

		/* Call device function */
		alignSequences(querySeq_d, refSeq_d,
				maxQryLen, refSeqLength, match, mismatch, gapOpenPenalty,
				gapExtPenalty, score_d, alignStart_d, alignEnd_d,
				alignLength_d, queryConsensus_d, refConsensus_d,
				queriesPerIteration, refsPerIteration, maxThreadsPerBlock);

		/* Copy device memory back to host memory */
		copyMemDeviceToHost(score, score_d, alignStart, alignStart_d, alignEnd,
				alignEnd_d, alignLength, alignLength_d, queryConsensus,
				queryConsensus_d, refConsensus, refConsensus_d,
				queriesPerIteration, refsPerIteration, combinedLength);

		/* Print results */
		if (outFormat == SAM_OUTPUT_FORMAT)
		{
			printOutputSam(queriesPerIteration, queryName, refName, score,
					alignStart, alignEnd, alignLength, combinedLength,
					queryConsensus, refConsensus, outputFilePtr);
		}
		else
		{
			chooseBestAlignment(queryName, MAX_SEQ_NAME_LENGTH, score,
					queriesPerIteration);
			printOutput(queriesPerIteration, queryName, refName, score,
					refDistance, alignStart, alignEnd, alignLength,
					combinedLength, queryConsensus, refConsensus,
					outputFilePtr);
		}

		time(&endTime);
		diffTime = difftime(endTime, startTime);
		fprintf(stderr, "   Completed %d / %d query iterations...%.2lf secs\n",
				(i + 1), numQryIterations, diffTime);
	}
	fclose(queryFilePtr);
	fclose(refFilePtr);
	fclose(tmpFilePtr);
	fclose(outputFilePtr);
//	refMapFree();

	/* Free resources on device */
	fprintf(stderr, "Freeing device resources...\n");
	cudaFree(refSeq_d);
	cudaFree(querySeq_d);
	cudaFree(score_d);
	cudaFree(alignStart_d);
	cudaFree(alignEnd_d);
	cudaFree(alignLength_d);
	cudaFree(refConsensus_d);
	cudaFree(queryConsensus_d);
	PRINT_CUDA_ERROR()

	/* Free resources on host */
	fprintf(stderr, "Freeing host resources...\n");
	free(queryName);
	free(querySeq);
	free(refName);
	free(refSeq);
	free(score);
	free(alignStart);
	free(alignEnd);
	free(alignLength);
	free(refConsensus);
	free(queryConsensus);
	free(refDistance);
}


/**
 * This is a wrapper function that wraps @a chooseBestAlignment. This function
 * has been added so that @a chooseBestAlignment can be unit-tested.
 *
 * Finds the best alignment score for a given query and sets the score
 * of the queries with lower scores to -1.
 *
 * @param			qryName		Query names.
 * @param			qryNameLen	Max length of each query name.
 * @param[in,out]	score		Alignment scores.
 * @param			size		The number of elements in each of @a qryName
 * and @a score.
 */
void chooseBestAlignment_wrap(char *qryName, int qryNameLen, float *score,
		int size)
{
	chooseBestAlignment(qryName, qryNameLen, score, size);
}


/**
 * Finds the best alignment score for a given query and sets the score
 * of the queries with lower scores to -1.
 *
 * It is assumed that the query names and scores of the same query are
 * consecutive.
 *
 * @param			qryName		Query names.
 * @param			qryNameLen	Max length of each query name.
 * @param[in,out]	score		Alignment scores.
 * @param			size		The number of elements in each of @a qryName
 * and @a score.
 */
static void chooseBestAlignment(char *qryName, int qryNameLen, float *score,
		int size)
{
#define	MAX_QRY_HITS	20	/* Maximum number of hits a query will have. */

	static float maxScore;
	static int i, j, randNum, tmp, start, end;
	static int idxArr[MAX_QRY_HITS];

	start = 0;
	end = size - 1;
	while (start < size)
	{
		/* Find the starting and ending index for a given query. */
		for (i = start + 1; i < size; ++i)
		{
			if (qryNameCmp(qryName + (i * qryNameLen),
					qryName + ((i - 1) * qryNameLen), qryNameLen) != 0)
			{
				end = i - 1;
				break;
			}
		}
		if (end < start)
			end = size - 1;


		/* Find the max score. */
		maxScore = 0.0;
		for (i = start; i <= end; ++i)
		{
			if (maxScore < score[i])
				maxScore = score[i];
		}


		/* Keep track of the indexes of those scores that are equal to max
		 * score. */
		j = 0;
		for (i = start; i <= end; ++i)
		{
			if (score[i] == maxScore)
			{
				idxArr[j] = i;
				++j;
			}
		}
		randNum = rand() % j;


		/* Set the score of queries with indexes not equal to 'randNum' to
		 * -1.0. */
		tmp = idxArr[randNum];
		for (i = start; i <= end; ++i)
		{
			if (i != tmp)
				score[i] = -1.0;
		}

		start = end + 1;
	}
}


/**
 * Finds the best alignment score for a given query and sets the score
 * of the queries with lower scores to -1.
 *
 * It is assumed that the query names and scores of the same query are
 * consecutive.
 *
 * @param			qryName		Query names.
 * @param			qryNameLen	Max length of each query name.
 * @param[in,out]	score		Alignment scores.
 * @param			size		The number of elements in each of @a qryName
 * and @a score.
 */
static void chooseBestAlignment2(char *qryName, int qryNameLen, float *score,
		int size)
{
#define	MAX_QRY_HITS	20	/* Maximum number of hits a query will have. */

	float maxScore;
	int i, j, randNum, tmp, start, end;
	int idxArr[MAX_QRY_HITS];

	start = 0;
	end = size - 1;
	while (start < size)
	{
		/* Find the starting and ending index for a given query. */
		for (i = start + 1; i < size; ++i)
		{
			if (qryNameCmp(qryName + (i * qryNameLen),
					qryName + ((i - 1) * qryNameLen), qryNameLen) != 0)
			{
				end = i - 1;
				break;
			}
		}
		if (end < start)
			end = size - 1;


		/* Find the max score. */
		maxScore = 0.0;
		for (i = start; i <= end; ++i)
		{
			if (maxScore < score[i])
				maxScore = score[i];
		}


		/* Keep track of the indexes of those scores that are equal to max
		 * score. */
		j = 0;
		for (i = start; i <= end; ++i)
		{
			if (score[i] == maxScore)
			{
				idxArr[j] = i;
				++j;
			}
		}
		randNum = rand() % j;


		/* Set the score of queries with indexes not equal to 'randNum' to
		 * -1.0. */
		tmp = idxArr[randNum];
		for (i = start; i <= end; ++i)
		{
			if (i != tmp)
				score[i] = -1.0;
		}

		start = end + 1;
	}
}


/**
 * This is a wrapper function that wraps @a chooseBestAlignment_paired.
 * This function has been added so that @a chooseBestAlignment_paired can be
 * unit-tested.
 *
 * Finds the best alignment score for a given query-pair and sets the score
 * of the pairs with lower scores to -1.
 *
 * @param			qryName		Query names.
 * @param			qryNameLen	Max length of each query name.
 * @param[in,out]	score		Alignment scores.
 * @param			size		The number of elements in each of @a qryName
 * and @a score.
 */
void chooseBestAlignment_paired_wrap(char *qryName, int qryNameLen,
		float *score, int size)
{
	chooseBestAlignment_paired(qryName, qryNameLen, score, size);
}


/**
 * Finds the best-scoring query-pair and sets the score of the pairs with
 * lower scores to -1.
 *
 * It is assumed that the same query pairs are consecutive.
 *
 * @param			qryName		Query names.
 * @param			qryNameLen	Max length of each query name.
 * @param[in,out]	score		Alignment scores.
 * @param			size		The number of elements in each of @a qryName
 * and @a score.
 */
static void chooseBestAlignment_paired(char *qryName, int qryNameLen,
		float *score, int size)
{
#define	MAX_QRY_HITS	20	/* Maximum number of hits a query-pair will have. */
	static int start, end, i, j, randNum, tmp1, tmp2;
	static float maxScore;
	static int idxArr[MAX_QRY_HITS];
	start = 0;
	end = size - 1;
	while (start < size)
	{
		/* Find the starting and ending index for a given query. */
		for (i = start + 2; i < size; i = i + 2)
		{
			if (qryNameCmp(qryName + (i * qryNameLen),
					qryName + ((i - 2) * qryNameLen), qryNameLen) != 0)
			{
				end = i - 1;
				break;
			}
		}
		if (end < start)
			end = size - 1;

		/* Find the max score. */
		maxScore = 0.0;
		for (i = start; i <= end; i = i + 2)
		{
			if (maxScore < (score[i] + score[i + 1]))
				maxScore = (score[i] + score[i + 1]);
		}

		/* Keep track of the indexes of those scores that are equal to max
		 * score. */
		j = 0;
		for (i = start; i <= end; i = i + 2)
		{
			if ((score[i] + score[i + 1]) == maxScore)
			{
				idxArr[j] = i;
				++j;
			}
		}
		randNum = rand() % j;

		/* Set the score of queries with indexes not equal to 'randNum' to
		 * -1.0. */
		tmp1 = idxArr[randNum];
		tmp2 = tmp1 + 1;
		for (i = start; i <= end; ++i)
		{
			if (i != tmp1 && i != tmp2)
				score[i] = -1.0;
		}

		start = end + 1;
	}
}


/**
 * This is a wrapper function that wraps @a chooseBestAlignment_paired2.
 * This function has been added so that @a chooseBestAlignment_paired2 can be
 * unit-tested.
 *
 * Finds the best alignment score for a given query-pair and sets the score
 * of the pairs with lower scores to -1.
 *
 * @param			qryName		Query names.
 * @param			qryNameLen	Max length of each query name.
 * @param[in,out]	score		Alignment scores.
 * @param			alignStart	Alignment start positions on the reference.
 * @param			alignEnd	Alignment end positions on the reference.
 * @param			size		The number of elements in each of @a qryName
 * and @a score.
 * @param			minFragSize	Minimum fragment size.
 * @param			maxFragSize	Maximum fragment size.
 */
void chooseBestAlignment_paired2_wrap(char *qryName, int qryNameLen,
		float *score, int *alignStart, int *alignEnd, int size,
		uint minFragSize, uint maxFragSize)
{
	chooseBestAlignment_paired2(qryName, qryNameLen, score, alignStart,
			alignEnd, size, minFragSize, maxFragSize);
}


/**
 * Finds the best-scoring query-pair and sets the score of the pairs with
 * lower scores to -1.
 *
 * It is assumed that the same query pairs are consecutive.
 *
 * @param			qryName		Query names.
 * @param			qryNameLen	Max length of each query name.
 * @param[in,out]	score		Alignment scores.
 * @param			alignStart	Alignment start positions on the reference.
 * @param			alignEnd	Alignment end positions on the reference.
 * @param			size		The number of elements in each of @a qryName
 * @param			minFragSize	Minimum fragment size.
 * @param			maxFragSize	Maximum fragment size.
 * and @a score.
 */
static void chooseBestAlignment_paired2(char *qryName, int qryNameLen,
		float *score, int *alignStart, int *alignEnd, int size,
		uint minFragSize, uint maxFragSize)
{
#define	MAX_QRY_PAIR_HITS	(2 * MAX_NUM_HITS_PAIRED_END)	/* Maximum number of
hits a query-pair will have. */
	static int start, end, i, j, randNum, tmp1, tmp2, distance;
	static float maxScore;
	static int idxArr[MAX_QRY_PAIR_HITS];
	start = 0;
	end = size - 1;
	while (start < size)
	{
		/* Find the starting and ending index for a given query. */
		for (i = start + 2; i < size; i = i + 2)
		{
			if (qryNameCmp(qryName + (i * qryNameLen),
					qryName + ((i - 2) * qryNameLen), qryNameLen) != 0)
			{
				end = i - 1;
				break;
			}
		}
		if (end < start)
			end = size - 1;

		/* Find the max score. */
		maxScore = 0.0;
		for (i = start; i <= end; i = i + 2)
		{
			if (maxScore < (score[i] + score[i + 1]))
			{
				distance = alignEnd[i + 1] - alignStart[i];
				if (distance >= minFragSize && distance <= maxFragSize)
					maxScore = (score[i] + score[i + 1]);
			}
		}

		/* Keep track of the indexes of those scores that are equal to max
		 * score. */
		j = 0;
		for (i = start; i <= end; i = i + 2)
		{
			if ((score[i] + score[i + 1]) == maxScore)
			{
				distance = alignEnd[i + 1] - alignStart[i];
				if (distance >= minFragSize && distance <= maxFragSize)
				{
					idxArr[j] = i;
					++j;
				}
			}
		}

		/* Set the score of queries with indexes not equal to 'randNum' to
		 * -1.0. */
		if (j > 0)
		{
			randNum = rand() % j;
			tmp1 = idxArr[randNum];
			tmp2 = tmp1 + 1;
			for (i = start; i <= end; ++i)
			{
				if (i != tmp1 && i != tmp2)
					score[i] = -1.0;
			}
		}
		else
		{
			for (i = start; i <= end; ++i)
				score[i] = -1.0;
		}

		start = end + 1;
	}
}


/**
 * This is a wrapper function that wraps @a chooseBestAlignment_paired3.
 * It has been added so that @a chooseBestAlignment_paired3 can be unit-tested.
 *
 * @param			qryName		Query names.
 * @param			qryNameLen	Max length of each query name.
 * @param			refName		Reference name.
 * @param			revComp		Indicates whether the query is a reverse
 * complement.
 * @param			refNameLen	Reference name length.
 * @param[in,out]	score		Alignment scores.
 * @param			refDistance	Reference distances.
 * @param			alignStart	Alignment start positions on the reference.
 * @param			alignEnd	Alignment end positions on the reference.
 * @param			size		The number of elements in each of @a qryName
 * @param			minFragSize	Minimum fragment size.
 * @param			maxFragSize	Maximum fragment size.
 */
void chooseBestAlignment_paired3_wrap(char *qryName, int qryNameLen,
		uint *revComp, char *refName, int refNameLen, float *score,
		int *refDistance, int *alignStart, int *alignEnd, int size,
		uint minFragSize, uint maxFragSize)
{
	chooseBestAlignment_paired3(qryName, qryNameLen, revComp, refName,
			refNameLen, score, refDistance, alignStart, alignEnd, size,
			minFragSize, maxFragSize);
}


/**
 * Finds the best-scoring query-pair and sets the score of the pairs with
 * lower scores to -1.
 *
 * It is assumed that the same query pairs are consecutive.
 *
 * @param			qryName		Query names.
 * @param			qryNameLen	Max length of each query name.
 * @param			refName		Reference name.
 * @param			revComp		Indicates whether the query is a reverse
 * complement.
 * @param			refNameLen	Reference name length.
 * @param[in,out]	score		Alignment scores.
 * @param			refDistance	Reference distances.
 * @param			alignStart	Alignment start positions on the reference.
 * @param			alignEnd	Alignment end positions on the reference.
 * @param			size		The number of elements in each of @a qryName
 * @param			minFragSize	Minimum fragment size.
 * @param			maxFragSize	Maximum fragment size.
 */
static void chooseBestAlignment_paired3(char *qryName, int qryNameLen,
		uint *revComp, char *refName, int refNameLen, float *score,
		int *refDistance, int *alignStart, int *alignEnd, int size,
		uint minFragSize, uint maxFragSize)
{
#define	MAX_QRY_PAIR_HITS	(2 * MAX_NUM_HITS_PAIRED_END)	/* Maximum number of
hits a query-pair will have. */
	static int start, end, start2, end2, i, j, k, randNum, tmp1, tmp2, distance;
	static float maxScore;
	static int idxArr1[MAX_QRY_PAIR_HITS], idxArr2[MAX_QRY_PAIR_HITS];
	static int hasMatePair;
	static char mateName[MAX_QRY_NAME_LENGTH];
	start = 0;
	end = size - 1;
	while (start < size)
	{
		/* Find the starting and ending index for a given query. */
		for (i = start + 1; i < size; ++i)
		{
			if (qryNameCmp(qryName + (i * qryNameLen),
					qryName + ((i - 1) * qryNameLen), qryNameLen) != 0)
			{
				end = i - 1;
				break;
			}
		}
		if (end < start)
			end = size - 1;


		/* Find the starting and ending index of the mate pair. */
		getMatePairName(qryName + (start * qryNameLen), mateName);
		start2 = end + 1;
		end2 = size - 1;
		hasMatePair = 0;
		if (start2 < size)
		{
			for (i = start2; i < size; ++i)
			{
				if (qryNameCmp(mateName, qryName + (i * qryNameLen),
						qryNameLen) != 0)
				{
					end2 = i - 1;
					break;
				}
				else
					hasMatePair = 1;
			}

			if (end2 < start2)
			{
				/* Query has a mate. */
				if (hasMatePair == 1)
					end2 = size - 1;
				/* Query doesn't have a mate. */
				else
				{
					for (i = start; i <= end; ++i)
						score[i] = -1.0;
					start = start2;
					continue;
				}
			}
		}
		else
		{
			for (i = start; i <= end; ++i)
				score[i] = -1.0;
			start = start2;
			continue;
		}


		/* For each mate-pair, find the max score, subject to contraints. */
		maxScore = 0.0;
		for (i = start; i <= end; ++i)
		{
			for (j = start2; j <= end2; ++j)
			{
				if (((revComp[i] == 0 && revComp[j] == 1)
						|| (revComp[j] == 0 && revComp[i] == 1))
						&& maxScore < (score[i] + score[j]))
				{
					if (strcmp(refName + (i * refNameLen),
							refName + (j * refNameLen)) == 0)
					{
						if (revComp[i] == 0 && revComp[j] == 1)
						{
							distance = (refDistance[j] + alignEnd[j])
									- (refDistance[i] + alignStart[i]);
						}
						else
						{
							distance = (refDistance[i] + alignEnd[i])
									- (refDistance[j] + alignStart[j]);
						}
						if (distance >= minFragSize && distance <= maxFragSize)
							maxScore = score[i] + score[j];
					}
					else
						maxScore = score[i] + score[j];
				}
			}
		}

		/* Keep track of the indexes of mate-pairs that have the max scores,
		 * subject to contraints. */
		k = 0;
		for (i = start; i <= end; ++i)
		{
			for (j = start2; j <= end2; ++j)
			{
				if (((revComp[i] == 0 && revComp[j] == 1)
						|| (revComp[j] == 0 && revComp[i] == 1))
						&& (score[i] + score[j]) == maxScore)
				{
					if (strcmp(refName + (i * refNameLen),
							refName + (j * refNameLen)) == 0)
					{
						if (revComp[i] == 0 && revComp[j] == 1)
						{
							distance = (refDistance[j] + alignEnd[j])
									- (refDistance[i] + alignStart[i]);
						}
						else
						{
							distance = (refDistance[i] + alignEnd[i])
									- (refDistance[j] + alignStart[j]);
						}
						if (distance >= minFragSize && distance <= maxFragSize)
						{
							idxArr1[k] = i;
							idxArr2[k] = j;
							++k;
						}
					}
					else
					{
						idxArr1[k] = i;
						idxArr2[k] = j;
						++k;
					}
				}
			}
		}

		/* Set the score to -1.0 of those pairs that do not have max-scores. */
		if (k > 0)
		{
			randNum = rand() % k;
			tmp1 = idxArr1[randNum];
			tmp2 = idxArr2[randNum];
			for (i = start; i <= end; ++i)
			{
				if (i != tmp1)
					score[i] = -1.0;
			}

			for (i = start2; i <= end2; ++i)
			{
				if (i != tmp2)
					score[i] = -1.0;
			}
		}
		else
		{
			for (i = start; i <= end; ++i)
				score[i] = -1.0;

			for (i = start2; i <= end2; ++i)
				score[i] = -1.0;
		}

		start = end2 + 1;
	}
}


/**
 * Returns the mate pair name of the given query name. The mate pair names are
 * similar, except that the first query of the pair is suffixed with "/1"
 * and the second query of the pair is suffixed with "/2".
 *
 * @param		qryName		Name of a query in the pair.
 * @param[out]	mateName	Name of the mate pair.
 */
void getMatePairName(const char *qryName, char *mateName)
{
	static int length;
	length = strlen(qryName);
	--length;
	if (qryName[length] == '*')
	{
		if (qryName[length - 1] == '1')
		{
			strncpy(mateName, qryName, length - 1);
			mateName[length - 1] = '2';
			mateName[length] = '\0';
		}
		else if (qryName[length - 1] == '2')
		{
			strncpy(mateName, qryName, length - 1);
			mateName[length - 1] = '1';
			mateName[length] = '\0';
		}
	}
	else
	{
		if (qryName[length] == '1')
		{
			strncpy(mateName, qryName, length);
			mateName[length] = '2';
			mateName[length + 1] = '\0';
		}
		else if (qryName[length] == '2')
		{
			strncpy(mateName, qryName, length);
			mateName[length] = '1';
			mateName[length + 1] = '\0';
		}
	}
}


/**
 * This is a wrapper function that wraps @a qryNameCmp function. This function
 * has been added so that @a qryNameCmp can be unit-tested.
 *
 * Compares the two query names and returns 0 if they are same (asterisk at
 * the end of query name is ignored) and -1 if they are not same.
 *
 * @param	s		First query name.
 * @param	t		Second query name.
 * @param	size	Number of characters in each name.
 * @return	0 if both query names are same (asterisk '*' at the end of the
 * query name is ignored), -1 if they are not same.
 */
int qryNameCmp_wrap(char *s, char *t, int size)
{
	return qryNameCmp(s, t, size);
}


/**
 * Compares the two query names and returns 0 if they are same (asterisk at
 * the end of query name is ignored) and -1 if they are not same..
 *
 * @param	s		First query name.
 * @param	t		Second query name.
 * @param	size	Number of characters in each name.
 * @return	0 if both query names are same (asterisk '*' at the end of the
 * query name is ignored), -1 if they are not same.
 */
static int qryNameCmp(char *s, char *t, int size)
{
	static int i;

	for (i = 0; i < size; ++i)
	{
		if (s[i] == '\0' && t[i] == '\0')
			return 0;
		else if (s[i] == '\0' || t[i] == '\0')
			break;
		else if (s[i] != t[i])
			return -1;
	}

	if ((s[i] == '*' && s[i + 1] == '\0') || (t[i] == '*' && t[i + 1] == '\0'))
		return 0;
	else
		return -1;
}


/**
 * Aligns queries to references.
 *
 * @param qryFile1 Path to the query file1.
 * @param qryFile2 Path to the query file2.
 * @param numQrs Number of queries in the query file.
 * @param qryLen Length of each query sequence.
 * @param refFile Path to the reference file.
 * @param numRefs Number of reference sequences in reference file.
 * @param refLen Length of each reference sequence.
 * @param matchFile Temporary file in which file offsets of each query and
 * its best-matching reference sequences have been stored.
 * @param outFormat The output format in which the alignments would be printed.
 * @param match Match score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param outputFile Output file where alignments will be stored.
 * @param minFragSize Minimum fragment size.
 * @param maxFragSize Maximum fragment size.
 */
void alignQueries_paired(char *qry1File, char *qry2File, int numQrs,
		int qryLen, char *refFile, int numRefs, int refLen, char *matchFile,
		int outFormat, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, const char *outputFile, uint minFragSize,
		uint maxFragSize)
{
	FILE *qry1FilePtr = fopen(qry1File, "r");
	FILE *qry2FilePtr = fopen(qry2File, "r");
	FILE *refFilePtr = fopen(refFile, "r");

	int numMatches = getNumMatches(matchFile);
	numRefs = numMatches;
	refLen = 3 * qryLen;


	/*
	 * Get GPU specifications (max number of threads per block,
	 * max number of blocks per grid dimension, etc.) using CUDA API
	 */
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	int device;
	cudaDeviceProp deviceProp[deviceCount];
	for (device = 0; device < deviceCount; ++device)
		cudaGetDeviceProperties(&deviceProp[device], device);
	cudaSetDevice(SELECTED_DEVICE_INDEX);
	PRINT_CUDA_ERROR()

	int maxRegsPerBlock = deviceProp[SELECTED_DEVICE_INDEX].regsPerBlock;
	int maxThreadsPerBlock = (int) floor((float) maxRegsPerBlock / NUM_REGS_USED);
	if (maxThreadsPerBlock > deviceProp[SELECTED_DEVICE_INDEX].maxThreadsPerBlock)
		maxThreadsPerBlock = deviceProp[SELECTED_DEVICE_INDEX].maxThreadsPerBlock;

	int memPerRefQueryPair = getMemPerRefQryPair(refLen, qryLen);
	int maxQrsPerIteration = getNumQueriesPerIteration_paired(deviceProp,
			memPerRefQueryPair, numMatches);
	int numQryIterations = (int) ceil((float) numMatches / maxQrsPerIteration);


	/* Allocate memory on host */
	fprintf(stderr, "Allocating memory on host...\n");
	int combinedLength = refLen + qryLen + 2;
	float *score;
	char *refSeq;
	char *qrySeq;
	uint *revComp;
	int *alignStart;
	int *alignEnd;
	int *alignLen;
	int *refDistance;
	char *refConsensus;
	char *qryConsensus;
	char *qryName;
	char *refName;
	allocateHostMem(&score, &refDistance, &alignStart, &alignEnd, &alignLen,
			&qryConsensus, &refConsensus, &qryName, &qrySeq, &refName,
			&refSeq, maxQrsPerIteration, MAX_QRY_NAME_LENGTH,
			MAX_REF_NAME_LENGTH, qryLen, refLen, combinedLength);
	revComp = (uint *) malloc(maxQrsPerIteration * sizeof(uint));

	/* Allocate memory on device */
	fprintf(stderr, "Allocating memory on device...\n");
	char *refSeq_d;
	char *qrySeq_d;
	float *score_d;
	float *hMatrix_d;
	int *cellBacktracker_d;
	int *alignStart_d;
	int *alignEnd_d;
	int *alignLen_d;
	char *refConsensus_d;
	char *qryConsensus_d;
	allocateDeviceMem(&qrySeq_d, &refSeq_d, &score_d, &hMatrix_d,
			&cellBacktracker_d, &alignStart_d, &alignEnd_d, &alignLen_d,
			&qryConsensus_d, &refConsensus_d, maxQrsPerIteration,
			qryLen, refLen, combinedLength);

	fprintf(stderr, "Computing alignments...\n");
	ulong refNameOffset;
	ulong refSeqOffset;
	ulong refSeqFragOffset;
	ulong qryNameOffset;
	ulong qrySeqOffset;
	FILE *tmpFilePtr = fopen(matchFile, "r");
	FILE *outputFilePtr = fopen(outputFile, "w");
	int i, j, refIdx;
	int refsPerIteration = 0, qrsPerIteration = 0;
	char tmpQrySeq[qryLen];
	int numShiftBases = refLen / 3;
	int mateIdx;

	refMapCreate(refFile);

	for (i = 0; i < numQryIterations; ++i)
	{
		/* Get number of queries for this iteration */
		qrsPerIteration = getNumSeqsInCurrentIteration(i, numMatches,
				numQryIterations, maxQrsPerIteration);

		/* Fetch reference and query sequences */
		for (j = 0; j < qrsPerIteration; ++j)
		{
			getOffsets_paired(tmpFilePtr, &qryNameOffset, &qrySeqOffset,
					&mateIdx, &refNameOffset, &refSeqOffset, &refSeqFragOffset,
					refDistance + j, revComp + j, &refIdx);

			refDistance[j] = max(0, refDistance[j] - numShiftBases);
			refSeqFragOffset = refMapGetFileOffset(refNameOffset,
					refDistance[j]);

			getSeq(refFilePtr, refSeqFragOffset, refSeq + (j * refLen),
					refLen);
			getSeqName(refFilePtr, refNameOffset, refName + (j
					* MAX_SEQ_NAME_LENGTH));


			if (mateIdx == 1)
			{
				if (revComp[j] == 1)
				{
					getSeq(qry1FilePtr, qrySeqOffset, tmpQrySeq, qryLen);
					qryGetReverseComplement(tmpQrySeq, qryLen,
							qrySeq + (j * qryLen));
					getSeqName(qry1FilePtr, qryNameOffset, qryName + (j
										* MAX_SEQ_NAME_LENGTH));
					strcat(qryName + (j * MAX_SEQ_NAME_LENGTH), "*");
				}
				else
				{
					getSeq(qry1FilePtr, qrySeqOffset, qrySeq
							+ (j * qryLen), qryLen);
					getSeqName(qry1FilePtr, qryNameOffset, qryName + (j
							* MAX_SEQ_NAME_LENGTH));
				}
			}
			else
			{
				if (revComp[j] == 1)
				{
					getSeq(qry2FilePtr, qrySeqOffset, tmpQrySeq, qryLen);
					qryGetReverseComplement(tmpQrySeq, qryLen,
							qrySeq + (j * qryLen));
					getSeqName(qry2FilePtr, qryNameOffset, qryName + (j
										* MAX_SEQ_NAME_LENGTH));
					strcat(qryName + (j * MAX_SEQ_NAME_LENGTH), "*");
				}
				else
				{
					getSeq(qry2FilePtr, qrySeqOffset, qrySeq
							+ (j * qryLen), qryLen);
					getSeqName(qry2FilePtr, qryNameOffset, qryName + (j
							* MAX_SEQ_NAME_LENGTH));
				}
			}
		}

		/* Copy queries from host to device */
		fprintf(stderr, "   Copying query sequences from host to device...\n");
		copyMemHostToDevice(qrySeq_d, qrySeq, qrsPerIteration
				* qryLen * sizeof(char));

		/* Copy references from host to device */
		fprintf(stderr, "   Copying ref sequences from host to device...\n");
		copyMemHostToDevice(refSeq_d, refSeq, qrsPerIteration
				* refLen * sizeof(char));

		cudaMemset(score_d, 0, qrsPerIteration * sizeof(float));
		PRINT_CUDA_ERROR()

		/* Call device function */
		alignSequences(qrySeq_d, refSeq_d, qryLen, refLen, match, mismatch,
				gapOpenPenalty, gapExtPenalty, score_d, alignStart_d,
				alignEnd_d, alignLen_d, qryConsensus_d, refConsensus_d,
				qrsPerIteration, refsPerIteration, maxThreadsPerBlock);

		/* Copy device memory back to host memory */
		fprintf(stderr, "   Copying results from device to host...\n");
		copyMemDeviceToHost(score, score_d, alignStart, alignStart_d, alignEnd,
				alignEnd_d, alignLen, alignLen_d, qryConsensus,
				qryConsensus_d, refConsensus, refConsensus_d,
				qrsPerIteration, refsPerIteration, combinedLength);

		/* Print results */
		fprintf(stderr, "   Printing results to output file...\n");
		if (outFormat == SAM_OUTPUT_FORMAT)
		{
			printOutputSam(qrsPerIteration, qryName, refName, score,
					alignStart, alignEnd, alignLen, combinedLength,
					qryConsensus, refConsensus, outputFilePtr);
		}
		else
		{
			chooseBestAlignment_paired3(qryName, MAX_QRY_NAME_LENGTH, revComp,
					refName, MAX_REF_NAME_LENGTH, score, refDistance,
					alignStart, alignEnd, qrsPerIteration, minFragSize,
					maxFragSize);

			printOutput_paired(qrsPerIteration, qryName, refName, score,
					refDistance, alignStart, alignEnd, alignLen,
					combinedLength, qryConsensus, refConsensus,
					outputFilePtr);

//			printOutput(qrsPerIteration, qryName, refName, score,
//					refDistance, alignStart, alignEnd, alignLen,
//					combinedLength, qryConsensus, refConsensus,
//					outputFilePtr);
		}

		fprintf(stderr, "   Completed %d / %d query iterations...\n\n", (i + 1),
				numQryIterations);
	}
	fclose(tmpFilePtr);
	fclose(outputFilePtr);
	refMapFree();

	/* Free resources on device */
	fprintf(stderr, "Freeing device resources...\n");
	cudaFree(refSeq_d);
	cudaFree(qrySeq_d);
	cudaFree(score_d);
	cudaFree(alignStart_d);
	cudaFree(alignEnd_d);
	cudaFree(alignLen_d);
	cudaFree(refConsensus_d);
	cudaFree(qryConsensus_d);
	PRINT_CUDA_ERROR()

	/* Free resources on host */
	fprintf(stderr, "Freeing host resources...\n");
	free(qryName);
	free(qrySeq);
	free(refName);
	free(refSeq);
	free(score);
	free(alignStart);
	free(alignEnd);
	free(alignLen);
	free(refConsensus);
	free(qryConsensus);
	free(refDistance);
	fprintf(stderr, "DONE!\n");
}


/**
 * This is a wrapper function that wraps @a getNumQrsPerIter_paired function.
 * It has been added so that @a getNumQrsPerIter_paired can be unit-tested.
 *
 * @param 		file				Map file name.
 * @param 		maxNumQrsPerIter	Maximum number of queries per iteration.
 * @param[out] 	numQrsPerIter		Array where number of queries per iteration
 * will be stored.
 * @param[out]	arrSize				Number of elements in @a numQrsPerIter.
 */
void getNumQrsPerIter_paired_wrap(const char *file, int maxNumQrsPerIter,
		int **numQrsPerIter, int *arrSize)
{
	getNumQrsPerIter_paired(file, maxNumQrsPerIter, numQrsPerIter, arrSize);
}


/**
 * Calculates the number of queries for each iteration. This function makes
 * sure that all hits belonging to a pair are in the same iteration.
 *
 * @param 		file				Map file name.
 * @param 		maxNumQrsPerIter	Maximum number of queries per iteration.
 * @param[out] 	numQrsPerIter		Array where number of queries per iteration
 * will be stored.
 * @param[out]	arrSize				Number of elements in @a numQrsPerIter.
 */
static void getNumQrsPerIter_paired(const char *file, int maxNumQrsPerIter,
		int **numQrsPerIter, int *arrSize)
{
	char line[MAX_LINE_LENGTH];
	int oldCount = 0, count = 0, id, oldId = 0;
	FILE *filePtr = fopen(file, "r");
	int isNewIter = FALSE;
	*arrSize = 0;
	while (fgets(line, MAX_LINE_LENGTH, filePtr) != NULL)
	{
		++count;
		sscanf(line, "%d\t%*s\t%*d\t%*d\t%*d\t%*d\t%*d\t%*u\t%*u\t%*u\t%*u", &id);

		if (isNewIter == TRUE)
		{
			isNewIter = FALSE;
			if (id == oldId)
				(*numQrsPerIter)[*arrSize - 1] = oldCount;
			else
				count = 1;
		}
		else if (id != oldId)
			oldCount = count - 1;

		if (count == maxNumQrsPerIter)
		{
			*numQrsPerIter = (int *) realloc((*numQrsPerIter),
					((*arrSize) + 1) * sizeof(int));
			(*numQrsPerIter)[*arrSize] = count;
			++(*arrSize);
			count = count - oldCount;
			isNewIter = TRUE;
		}

		oldId = id;
	}
	if (isNewIter == FALSE && count <= maxNumQrsPerIter)
	{
		*numQrsPerIter = (int *) realloc(*numQrsPerIter,
				((*arrSize) + 1) * sizeof(int));
		(*numQrsPerIter)[*arrSize] = count;
		++(*arrSize);
	}
}


/**
 * Aligns queries to references.
 *
 * @param qryFile1 Path to the query file1.
 * @param qryFile2 Path to the query file2.
 * @param qryLen Length of each query sequence.
 * @param refFile Path to the reference file.
 * @param matchFile Temporary file in which file offsets of each query and
 * its best-matching reference sequences have been stored.
 * @param outFormat The output format in which the alignments would be printed.
 * @param match Match score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param outputFile Output file where alignments will be stored.
 * @param minFragSize Minimum fragment size.
 * @param maxFragSize Maximum fragment size.
 * @param numMatches Total number of hits obtained in Phase 1.
 */
void alignQueries_paired2(char *qry1File, char *qry2File, int qryLen,
		char *refFile, char *matchFile, int outFormat, float match,
		float mismatch, float gapOpenPenalty, float gapExtPenalty,
		const char *outputFile, uint minFragSize, uint maxFragSize,
		int numMatches)
{
	FILE *qry1FilePtr = fopen(qry1File, "r");
	FILE *qry2FilePtr = fopen(qry2File, "r");
	FILE *refFilePtr = fopen(refFile, "r");

	printf("Number of matches = %d\n", numMatches);
	int refLen = 3 * qryLen;


	/*
	 * Get GPU specifications (max number of threads per block,
	 * max number of blocks per grid dimension, etc.) using CUDA API
	 */
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	int device;
	cudaDeviceProp deviceProp[deviceCount];
	for (device = 0; device < deviceCount; ++device)
		cudaGetDeviceProperties(&deviceProp[device], device);
	cudaSetDevice(SELECTED_DEVICE_INDEX);
	PRINT_CUDA_ERROR()

	int maxRegsPerBlock = deviceProp[SELECTED_DEVICE_INDEX].regsPerBlock;
	int maxThreadsPerBlock = (int) floor((float) maxRegsPerBlock / NUM_REGS_USED);
	if (maxThreadsPerBlock > deviceProp[SELECTED_DEVICE_INDEX].maxThreadsPerBlock)
		maxThreadsPerBlock = deviceProp[SELECTED_DEVICE_INDEX].maxThreadsPerBlock;

	int memPerRefQueryPair = getMemPerRefQryPair(refLen, qryLen);
	int maxQrsPerIteration = getNumQueriesPerIteration_paired(deviceProp,
			memPerRefQueryPair, numMatches);
	int numQryIterations;

	fprintf(stderr, "Calculating number of queries per batch...");
	int *numQrsPerIter = NULL;
	getNumQrsPerIter_paired(matchFile, maxQrsPerIteration, &numQrsPerIter,
			&numQryIterations);
	fprintf(stderr, "done.\n");


	/* Allocate memory on host */
	fprintf(stderr, "Allocating memory on host...");
	int combinedLength = refLen + qryLen + 2;
	float *score;
	char *refSeq;
	char *qrySeq;
	uint *revComp;
	int *alignStart;
	int *alignEnd;
	int *alignLen;
	int *refDistance;
	char *refConsensus;
	char *qryConsensus;
	char *qryName;
	char *refName;
	allocateHostMem(&score, &refDistance, &alignStart, &alignEnd, &alignLen,
			&qryConsensus, &refConsensus, &qryName, &qrySeq, &refName,
			&refSeq, maxQrsPerIteration, MAX_QRY_NAME_LENGTH,
			MAX_REF_NAME_LENGTH, qryLen, refLen, combinedLength);
	int *qryIdx = (int *) malloc(maxQrsPerIteration * sizeof(int));
	revComp = (uint *) malloc(maxQrsPerIteration * sizeof(uint));
	fprintf(stderr, "done.\n");

	/* Allocate memory on device */
	fprintf(stderr, "Allocating memory on device...");
	char *refSeq_d;
	char *qrySeq_d;
	float *score_d;
	float *hMatrix_d;
	int *cellBacktracker_d;
	int *alignStart_d;
	int *alignEnd_d;
	int *alignLen_d;
	char *refConsensus_d;
	char *qryConsensus_d;
	allocateDeviceMem(&qrySeq_d, &refSeq_d, &score_d, &hMatrix_d,
			&cellBacktracker_d, &alignStart_d, &alignEnd_d, &alignLen_d,
			&qryConsensus_d, &refConsensus_d, maxQrsPerIteration,
			qryLen, refLen, combinedLength);
	fprintf(stderr, "done.\n");

	uint refNameOffset;
	uint refSeqOffset;
	uint refSeqFragOffset;
	uint qryNameOffset;
	uint qrySeqOffset;
	FILE *tmpFilePtr = fopen(matchFile, "r");
	FILE *outputFilePtr = fopen(outputFile, "w");
	int i, j, refIdx;
	int refsPerIteration = 0, qrsPerIteration = 0;
	char tmpQrySeq[qryLen + 1];
	int numShiftBases = refLen / 3;
	int mateIdx;

//	fprintf(stderr, "Creating a map of the reference...");
//	refMapCreate(refFile);
//	fprintf(stderr, "done.\n");

	time_t startTime, endTime;
	double diffTime = 0.0;
	fprintf(stderr, "Computing alignments:\n");
	for (i = 0; i < numQryIterations; ++i)
	{
		time(&startTime);

		/* Get number of queries for this iteration */
		qrsPerIteration = numQrsPerIter[i];

		/* Fetch reference and query sequences */
		fprintf(stderr, "   Gathering reference and query sequences...");
		for (j = 0; j < qrsPerIteration; ++j)
		{
			getOffsets_paired2(tmpFilePtr, &qryNameOffset, &qrySeqOffset,
					&mateIdx, &refNameOffset, &refSeqOffset, &refSeqFragOffset,
					refDistance + j, revComp + j, &refIdx, qryIdx + j);

			refDistance[j] = max(0, refDistance[j] - numShiftBases);
			refMapGetOffsets(refIdx, refDistance[j], &refNameOffset,
					&refSeqFragOffset);

			getSeq(refFilePtr, refSeqFragOffset, refSeq + (j * (refLen + 1)),
					refLen);
			getSeqName(refFilePtr, refNameOffset, refName + (j
					* MAX_SEQ_NAME_LENGTH));

			if (mateIdx == 1)
			{
				if (revComp[j] == 1)
				{
					getSeq(qry1FilePtr, qrySeqOffset, tmpQrySeq, qryLen);
					qryGetReverseComplement(tmpQrySeq, qryLen,
							qrySeq + (j * (qryLen + 1)));
					getSeqName(qry1FilePtr, qryNameOffset, qryName + (j
										* MAX_SEQ_NAME_LENGTH));
					strcat(qryName + (j * MAX_SEQ_NAME_LENGTH), "*");
				}
				else
				{
					getSeq(qry1FilePtr, qrySeqOffset, qrySeq
							+ (j * (qryLen + 1)), qryLen);
					getSeqName(qry1FilePtr, qryNameOffset, qryName + (j
							* MAX_SEQ_NAME_LENGTH));
				}
			}
			else
			{
				if (revComp[j] == 1)
				{
					getSeq(qry2FilePtr, qrySeqOffset, tmpQrySeq, qryLen);
					qryGetReverseComplement(tmpQrySeq, qryLen,
							qrySeq + (j * (qryLen + 1)));
					getSeqName(qry2FilePtr, qryNameOffset, qryName + (j
										* MAX_SEQ_NAME_LENGTH));
					strcat(qryName + (j * MAX_SEQ_NAME_LENGTH), "*");
				}
				else
				{
					getSeq(qry2FilePtr, qrySeqOffset, qrySeq
							+ (j * (qryLen + 1)), qryLen);
					getSeqName(qry2FilePtr, qryNameOffset, qryName + (j
							* MAX_SEQ_NAME_LENGTH));
				}
			}
		}
		fprintf(stderr, "done.\n");

		/* Copy queries from host to device */
		fprintf(stderr, "   Copying query sequences from host to device...");
		copyMemHostToDevice(qrySeq_d, qrySeq, qrsPerIteration
				* (qryLen + 1) * sizeof(char));
		fprintf(stderr, "done.\n");

		/* Copy references from host to device */
		fprintf(stderr, "   Copying ref sequences from host to device...");
		copyMemHostToDevice(refSeq_d, refSeq, qrsPerIteration
				* (refLen + 1) * sizeof(char));
		fprintf(stderr, "done.\n");

		cudaMemset(score_d, 0, qrsPerIteration * sizeof(float));
		PRINT_CUDA_ERROR()

		/* Call device function */
		fprintf(stderr, "   Executing kernel function...");
		alignSequences(qrySeq_d, refSeq_d, qryLen, refLen, match, mismatch,
				gapOpenPenalty, gapExtPenalty, score_d, alignStart_d,
				alignEnd_d, alignLen_d, qryConsensus_d, refConsensus_d,
				qrsPerIteration, refsPerIteration, maxThreadsPerBlock);
		fprintf(stderr, "done.\n");

		/* Copy device memory back to host memory */
		fprintf(stderr, "   Copying results from device to host...");
		copyMemDeviceToHost(score, score_d, alignStart, alignStart_d, alignEnd,
				alignEnd_d, alignLen, alignLen_d, qryConsensus,
				qryConsensus_d, refConsensus, refConsensus_d,
				qrsPerIteration, refsPerIteration, combinedLength);
		fprintf(stderr, "done.\n");

		/* Print results */
		fprintf(stderr, "   Printing results to output file...");
		if (outFormat == SAM_OUTPUT_FORMAT)
		{
			printOutputSam(qrsPerIteration, qryName, refName, score,
					alignStart, alignEnd, alignLen, combinedLength,
					qryConsensus, refConsensus, outputFilePtr);
		}
		else
		{
			chooseBestAlignment_paired3(qryName, MAX_QRY_NAME_LENGTH, revComp,
					refName, MAX_REF_NAME_LENGTH, score, refDistance,
					alignStart, alignEnd, qrsPerIteration, minFragSize,
					maxFragSize);

			printOutput_paired2(qrsPerIteration, qryIdx, qryName, refName, score,
					refDistance, alignStart, alignEnd, alignLen,
					combinedLength, qryConsensus, refConsensus,
					outputFilePtr);

//			printOutput(qrsPerIteration, qryName, refName, score,
//					refDistance, alignStart, alignEnd, alignLen,
//					combinedLength, qryConsensus, refConsensus,
//					outputFilePtr);
		}
		fprintf(stderr, "done.\n");

		time(&endTime);
		diffTime = difftime(endTime, startTime);
		fprintf(stderr, "   Time = %.2lf secs\n", diffTime);
		fprintf(stderr, "   Completed %d / %d query iterations\n\n", (i + 1),
				numQryIterations);
	}
	fclose(tmpFilePtr);
	fclose(outputFilePtr);
//	refMapFree();

	/* Free resources on device */
	fprintf(stderr, "Freeing device resources...\n");
	cudaFree(refSeq_d);
	cudaFree(qrySeq_d);
	cudaFree(score_d);
	cudaFree(alignStart_d);
	cudaFree(alignEnd_d);
	cudaFree(alignLen_d);
	cudaFree(refConsensus_d);
	cudaFree(qryConsensus_d);
	PRINT_CUDA_ERROR()

	/* Free resources on host */
	fprintf(stderr, "Freeing host resources...\n");
	free(qryName);
	free(qrySeq);
	free(refName);
	free(refSeq);
	free(score);
	free(alignStart);
	free(alignEnd);
	free(alignLen);
	free(refConsensus);
	free(qryConsensus);
	free(refDistance);
	free(numQrsPerIter);
	free(qryIdx);
	free(revComp);
}


/**
 * Aligns queries to references using a CPU-based Smith-Waterman implementation.
 *
 * @param qryFile Query sequence file.
 * @param numQrs Total number of queries in @a qryFile.
 * @param qryLen Query sequence length.
 * @param refFile Reference sequence file.
 * @param matchFile File from where data is to be input.
 * @param outFormat Output format.
 * @param match Match score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param outputFile Output file.
 */
void alignQueries_cpu(char *qryFile, int numQrs, int qryLen, char *refFile,
		char *matchFile, int outFormat, float match, float mismatch,
		float gapOpenPenalty, float gapExtPenalty, const char *outputFile)
{
	FILE *qryFilePtr = fopen(qryFile, "r");
	FILE *refFilePtr = fopen(refFile, "r");

	int numMatches = getNumMatches(matchFile);
	fprintf(stderr, "**num matches = %d\n", numMatches);
	int refLen = 3 * qryLen;


	fprintf(stderr, "**************Phase II*******************\n");

	/* Allocate memory on host */
	fprintf(stderr, "Allocating memory on host...\n");
	int combinedLen = refLen + qryLen + 2;
	float *alignScore = (float *) malloc(sizeof(float));
	char *refSeq = (char *) malloc(refLen * sizeof(char));
	char *qrySeq = (char *) malloc(qryLen * sizeof(char));
	int *alignStart = (int *) malloc(sizeof(int));
	int *alignEnd = (int *) malloc(sizeof(int));
	int *alignLen = (int *) malloc(sizeof(int));
	int *refDistance = (int *) malloc(sizeof(int));
	char *refConsensus = (char *) malloc(combinedLen * sizeof(char));
	char *qryConsensus = (char *) malloc(combinedLen * sizeof(char));
	char *qryName = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
	char *refName = (char *) malloc(MAX_REF_NAME_LENGTH * sizeof(char));
	float *HMatrix = (float *) malloc((refLen + 1) * (qryLen + 1)
			* sizeof(float));
	int *rowBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
			* sizeof(int));
	int *colBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
			* sizeof(int));
	int row, col, index;
	for (row = 0; row <= qryLen; ++row)
	{
		for (col = 0; col <= refLen; ++col)
		{
			index = (row * refLen) + col;
			HMatrix[index] = 0.0f;
			rowBacktracker[index] = 0;
			colBacktracker[index] = 0;
		}
	}

	fprintf(stderr, "Computing alignments...\n");
	uint refNameOffset;
	uint refSeqOffset;
	uint refSeqFragOffset;
	uint qryNameOffset;
	uint qrySeqOffset;
	uint isReverseComplement;
	FILE *tmpFilePtr = fopen(matchFile, "r");
	FILE *outputFilePtr = fopen(outputFile, "w");
	int qrsPerIteration = 1;
	char tmpQrySeq[qryLen];
	int i, j, refIdx;
	int numShiftBases = refLen / 3;
	time_t startTime, endTime;
	double diffTime = 0.0, totalTime = 0.0;
	for (i = 0; i < numMatches; ++i)
	{
		time(&startTime);
		getOffsets(tmpFilePtr, &qryNameOffset, &qrySeqOffset, &refNameOffset,
				&refSeqOffset, &refSeqFragOffset, refDistance,
				&isReverseComplement, &refIdx);

		*refDistance = max(0, *refDistance - numShiftBases);
		refMapGetOffsets(refIdx, *refDistance, &refNameOffset, &refSeqFragOffset);

//		getSeq(refFilePtr, refSeqFragOffset, refSeq, refLen);
		getRefSeq2(refFilePtr, refSeqFragOffset, refSeq, refLen);
		getSeqName(refFilePtr, refNameOffset, refName);

		if (isReverseComplement == 1)
		{
			getSeq(qryFilePtr, qrySeqOffset, tmpQrySeq, qryLen);
			qryGetReverseComplement(tmpQrySeq, qryLen, qrySeq);
			getSeqName(qryFilePtr, qryNameOffset, qryName);
			strcat(qryName, "*");
		}
		else
		{
			getSeq(qryFilePtr, qrySeqOffset, qrySeq, qryLen);
			getSeqName(qryFilePtr, qryNameOffset, qryName);
		}

		*alignScore = 0.0f;
		*alignStart = 0;
		*alignEnd = 0;
		*alignLen = 0;
		for (row = 0; row <= qryLen; ++row)
		{
			for (col = 0; col <= refLen; ++col)
			{
				index = (row * refLen) + col;
				HMatrix[index] = 0.0f;
				rowBacktracker[index] = 0;
				colBacktracker[index] = 0;
			}
		}
		for (j = 0; j < combinedLen; ++j)
		{
			refConsensus[j] = '\0';
			qryConsensus[j] = '\0';
		}
		smithWaterman_cpu(HMatrix, rowBacktracker, colBacktracker, refSeq,
				refLen, qrySeq, qryLen, match, mismatch, gapOpenPenalty,
				gapExtPenalty, alignScore, alignStart, alignEnd, alignLen,
				refConsensus, qryConsensus);

		/* Print results */
		if (outFormat == SAM_OUTPUT_FORMAT)
		{
			printOutputSam(qrsPerIteration, qryName, refName, alignScore,
					alignStart, alignEnd, alignLen, combinedLen, qryConsensus,
					refConsensus, outputFilePtr);
		}
		else
		{
			chooseBestAlignment(qryName, MAX_SEQ_NAME_LENGTH, alignScore,
					qrsPerIteration);
			printOutput(qrsPerIteration, qryName, refName, alignScore,
					refDistance, alignStart, alignEnd, alignLen, combinedLen,
					qryConsensus, refConsensus, outputFilePtr);
		}

		time(&endTime);
		diffTime = difftime(endTime, startTime);
		totalTime += diffTime;
		if (i % 50000 == 0)
		{
			fprintf(stderr, "   Completed %d / %d query matches...%.2lf secs\n",
					(i + 1), numMatches, totalTime);
			totalTime = 0.0;
//			break;
		}
	}
	fclose(tmpFilePtr);
	fclose(outputFilePtr);

	/* Free resources */
	fprintf(stderr, "Freeing host resources...\n");
	free(qryName);
	free(qrySeq);
	free(refName);
	free(refSeq);
	free(alignScore);
	free(alignStart);
	free(alignEnd);
	free(alignLen);
	free(refConsensus);
	free(qryConsensus);
	free(refDistance);
}


typedef struct data
{
	float *HMatrix;
	int *rowBacktracker;
	int *colBacktracker;
	char *refName;
	char *refSeq;
	int refLen;
	char *qryName;
	char *qrySeq;
	int qryLen;
	float match;
	float mismatch;
	float gapOpenPenalty;
	float gapExtPenalty;
	float alignScore;
	int alignStart;
	int alignEnd;
	int alignLen;
	char *refConsensus;
	char *qryConsensus;
	int combinedLen;
	int refDistance;
} Data;


/**
 * Aligns queries to references using a multi-threaded CPU-based Smith-Waterman
 * implementation.
 *
 * @param qryFile Query sequence file.
 * @param numQrs Total number of queries in @a qryFile.
 * @param qryLen Query sequence length.
 * @param refFile Reference sequence file.
 * @param matchFile File from where data is to be input.
 * @param outFormat Output format.
 * @param match Match score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param outputFile Output file.
 */
void alignQueries_cpu_mt(char *qryFile, int numQrs, int qryLen, char *refFile,
		char *matchFile, int outFormat, float match, float mismatch,
		float gapOpenPenalty, float gapExtPenalty, const char *outputFile)
{
	FILE *qryFilePtr = fopen(qryFile, "r");
	FILE *refFilePtr = fopen(refFile, "r");

	int numMatches = getNumMatches(matchFile);
	fprintf(stderr, "**num matches = %d\n", numMatches);
	int refLen = REF_LENGTH;

	fprintf(stderr, "**************Phase II*******************\n");

	fprintf(stderr, "Computing alignments...\n");
	uint refNameOffset;
	uint refSeqOffset;
	uint refSeqFragOffset;
	uint qryNameOffset;
	uint qrySeqOffset;
	uint isReverseComplement;
	FILE *tmpFilePtr = fopen(matchFile, "r");
	FILE *outputFilePtr = fopen(outputFile, "w");
	int qrsPerIteration = 1;
	char tmpQrySeq[qryLen];
	int i, j, refIdx;
	int numShiftBases = refLen / 3;
	time_t startTime, endTime;
	double diffTime = 0.0, totalTime = 0.0;
	int numThreads = _NUM_THREADS;

	pthread_t thread[numThreads];
	Data **data = (Data **) malloc(numThreads * sizeof(Data *));
	for (i = 0; i < numThreads; ++i)
	{
		data[i] = (Data *) malloc(sizeof(Data));
		data[i]->HMatrix = (float *) malloc((refLen + 1) * (qryLen + 1)
				* sizeof(float));
		data[i]->rowBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
				* sizeof(int));
		data[i]->colBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
				* sizeof(int));
		data[i]->refName = (char *) malloc(MAX_REF_NAME_LENGTH * sizeof(char));
		data[i]->refSeq = (char *) malloc(refLen * sizeof(char));
		data[i]->refLen = refLen;
		data[i]->qryName = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
		data[i]->qrySeq = (char *) malloc(qryLen * sizeof(char));
		data[i]->qryLen = qryLen;
		data[i]->match = match;
		data[i]->mismatch = mismatch;
		data[i]->gapOpenPenalty = gapOpenPenalty;
		data[i]->gapExtPenalty = gapExtPenalty;
		data[i]->alignScore = 0.0f;
		data[i]->alignStart = 0;
		data[i]->alignEnd = 0;
		data[i]->alignLen = 0;
		data[i]->combinedLen = data[i]->refLen + data[i]->qryLen + 2;
		data[i]->refDistance = 0;
		data[i]->refConsensus = (char *) malloc(data[i]->combinedLen * sizeof(char));
		data[i]->qryConsensus = (char *) malloc(data[i]->combinedLen * sizeof(char));
	}

	for (i = 0; i < numMatches; i += numThreads)
	{
		time(&startTime);

		/* Create threads. */
		for (j = 0; j < numThreads; ++j)
		{
			getOffsets(tmpFilePtr, &qryNameOffset, &qrySeqOffset, &refNameOffset,
					&refSeqOffset, &refSeqFragOffset, &data[j]->refDistance,
					&isReverseComplement, &refIdx);

			data[j]->refDistance = max(0, data[j]->refDistance - numShiftBases);
			refMapGetOffsets(refIdx, data[j]->refDistance, &refNameOffset,
					&refSeqFragOffset);

			getRefSeq2(refFilePtr, refSeqFragOffset, data[j]->refSeq, refLen);
			getSeqName(refFilePtr, refNameOffset, data[j]->refName);

			if (isReverseComplement == 1)
			{
				getSeq(qryFilePtr, qrySeqOffset, tmpQrySeq, qryLen);
				qryGetReverseComplement(tmpQrySeq, qryLen, data[j]->qrySeq);
				getSeqName(qryFilePtr, qryNameOffset, data[j]->qryName);
				strcat(data[j]->qryName, "*");
			}
			else
			{
				getSeq(qryFilePtr, qrySeqOffset, data[j]->qrySeq, qryLen);
				getSeqName(qryFilePtr, qryNameOffset, data[j]->qryName);
			}

			pthread_create(&(thread[j]), NULL, getAlignments, data[j]);
		}

		/* Wait for threads to complete, then print the results. */
		for (j = 0; j < numThreads; ++j)
		{
			pthread_join(thread[j], NULL);
			if (outFormat == SAM_OUTPUT_FORMAT)
			{
				printOutputSam(qrsPerIteration, data[j]->qryName,
						data[j]->refName, &data[j]->alignScore,
						&data[j]->alignStart, &data[j]->alignEnd,
						&data[j]->alignLen, data[j]->combinedLen,
						data[j]->qryConsensus, data[j]->refConsensus,
						outputFilePtr);
			}
			else
			{
				chooseBestAlignment2(data[j]->qryName, MAX_SEQ_NAME_LENGTH,
						&data[j]->alignScore, qrsPerIteration);
				printOutput(qrsPerIteration, data[j]->qryName, data[j]->refName,
						&data[j]->alignScore, &data[j]->refDistance,
						&data[j]->alignStart, &data[j]->alignEnd,
						&data[j]->alignLen, data[j]->combinedLen,
						data[j]->qryConsensus, data[j]->refConsensus,
						outputFilePtr);
			}
		}

		time(&endTime);
		diffTime = difftime(endTime, startTime);
		totalTime += diffTime;
//		if ((i + numThreads) % 50000 == 0)
		if (i % 50000 == 0)
		{
			fprintf(stderr, "   Completed %d / %d query matches...%.2lf secs\n",
					(i + 1), numMatches, totalTime);
			totalTime = 0.0;
//			break;
		}
	}
	fclose(tmpFilePtr);
	fclose(outputFilePtr);

	/* Free resources */
	for (i = 0; i < numThreads; ++i)
	{
		free(data[i]->HMatrix);
		free(data[i]->rowBacktracker);
		free(data[i]->colBacktracker);
		free(data[i]->refName);
		free(data[i]->refSeq);
		free(data[i]->qryName);
		free(data[i]->qrySeq);
	}
	free(data);
}


/**
 * Get alignments using Smith-Waterman algorithm.
 *
 * @param arg	An argument.
 */
void *getAlignments(void *arg)
{
	Data *data = (Data *) arg;
	data->alignScore = 0.0f;
	data->alignStart = 0;
	data->alignEnd = 0;
	data->alignLen = 0;
	int i, j, index;
	for (i = 0; i <= data->qryLen; ++i)
	{
		for (j = 0; j <= data->refLen; ++j)
		{
			index = (i * data->refLen) + j;
			data->HMatrix[index] = 0.0f;
			data->rowBacktracker[index] = 0;
			data->colBacktracker[index] = 0;
		}
	}
	for (j = 0; j < data->combinedLen; ++j)
	{
		data->refConsensus[j] = '\0';
		data->qryConsensus[j] = '\0';
	}
	smithWaterman_cpu(data->HMatrix, data->rowBacktracker, data->colBacktracker,
			data->refSeq, data->refLen, data->qrySeq, data->qryLen, data->match,
			data->mismatch, data->gapOpenPenalty, data->gapExtPenalty,
			&data->alignScore, &data->alignStart, &data->alignEnd,
			&data->alignLen, data->refConsensus, data->qryConsensus);
	return NULL;
}


/**
 * Aligns queries to references using a CPU-based Smith-Waterman implementation.
 *
 * @param qry1File Query sequence file1.
 * @param qry2File Query sequence file2.
 * @param numQrs Total number of queries in @a qryFile.
 * @param qryLen Query sequence length.
 * @param refFile Reference sequence file.
 * @param numRefs Total number of references in @a refFile.
 * @param refLen Reference sequence length.
 * @param matchFile File from where data is to be input.
 * @param outFormat Output format.
 * @param match Match score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param outputFile Output file.
 */
void alignQueries_paired_cpu(char *qry1File, char *qry2File, int numQrs,
		int qryLen, char *refFile, int numRefs, int refLen, char *matchFile,
		int outFormat, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, const char *outputFile)
{
	FILE *qry1FilePtr = fopen(qry1File, "r");
	FILE *qry2FilePtr = fopen(qry2File, "r");
	FILE *refFilePtr = fopen(refFile, "r");

	int numMatches = getNumMatches(matchFile);
	fprintf(stderr, "**num matches = %d\n", numMatches);
	refLen = 3 * qryLen;


	fprintf(stderr, "**************Phase II*******************\n");

	/* Allocate memory on host */
	fprintf(stderr, "Allocating memory on host...\n");
	int combinedLen = refLen + qryLen + 2;
	float *alignScore = (float *) malloc(sizeof(float));
	char *refSeq = (char *) malloc(refLen * sizeof(char));
	char *qrySeq = (char *) malloc(qryLen * sizeof(char));
	int *alignStart = (int *) malloc(sizeof(int));
	int *alignEnd = (int *) malloc(sizeof(int));
	int *alignLen = (int *) malloc(sizeof(int));
	int *refDistance = (int *) malloc(sizeof(int));
	char *refConsensus = (char *) malloc(combinedLen * sizeof(char));
	char *qryConsensus = (char *) malloc(combinedLen * sizeof(char));
	char *qryName = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
	char *refName = (char *) malloc(MAX_REF_NAME_LENGTH * sizeof(char));
	float *HMatrix = (float *) malloc((refLen + 1) * (qryLen + 1)
			* sizeof(float));
	int *rowBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
			* sizeof(int));
	int *colBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
			* sizeof(int));
	int row, col, index;
	for (row = 0; row <= qryLen; ++row)
	{
		for (col = 0; col <= refLen; ++col)
		{
			index = (row * refLen) + col;
			HMatrix[index] = 0.0f;
			rowBacktracker[index] = 0;
			colBacktracker[index] = 0;
		}
	}

	fprintf(stderr, "Computing alignments...\n");
	uint refNameOffset;
	uint refSeqOffset;
	uint refSeqFragOffset;
	uint qryNameOffset;
	uint qrySeqOffset;
	uint isReverseComplement;
	FILE *tmpFilePtr = fopen(matchFile, "r");
	FILE *outputFilePtr = fopen(outputFile, "w");
	int qrsPerIteration = 1;
	char tmpQrySeq[qryLen];
	int i, j, isFirstQry = 0, refIdx;
	for (i = 0; i < numMatches; ++i)
	{
		getOffsets(tmpFilePtr, &qryNameOffset, &qrySeqOffset, &refNameOffset,
				&refSeqOffset, &refSeqFragOffset, refDistance,
				&isReverseComplement, &refIdx);

		getSeq(refFilePtr, refSeqFragOffset, refSeq, refLen);
		getSeqName(refFilePtr, refNameOffset, refName);

		if ((i % 2) == 0)
		{
			if (isReverseComplement == 1)
			{
				getSeq(qry1FilePtr, qrySeqOffset, tmpQrySeq, qryLen);
				qryGetReverseComplement(tmpQrySeq, qryLen, qrySeq);
				getSeqName(qry1FilePtr, qryNameOffset, qryName);
				strcat(qryName, "*");
			}
			else
			{
				getSeq(qry1FilePtr, qrySeqOffset, qrySeq, qryLen);
				getSeqName(qry1FilePtr, qryNameOffset, qryName);
			}
		}
		else
		{
			if (isReverseComplement == 1)
			{
				getSeq(qry2FilePtr, qrySeqOffset, tmpQrySeq, qryLen);
				qryGetReverseComplement(tmpQrySeq, qryLen, qrySeq);
				getSeqName(qry2FilePtr, qryNameOffset, qryName);
				strcat(qryName, "*");
			}
			else
			{
				getSeq(qry2FilePtr, qrySeqOffset, qrySeq, qryLen);
				getSeqName(qry2FilePtr, qryNameOffset, qryName);
			}
		}

		*alignScore = 0.0f;
		*alignStart = 0;
		*alignEnd = 0;
		*alignLen = 0;
		for (row = 0; row <= qryLen; ++row)
		{
			for (col = 0; col <= refLen; ++col)
			{
				index = (row * refLen) + col;
				HMatrix[index] = 0.0f;
				rowBacktracker[index] = 0;
				colBacktracker[index] = 0;
			}
		}
		for (j = 0; j < combinedLen; ++j)
		{
			refConsensus[j] = '\0';
			qryConsensus[j] = '\0';
		}
		smithWaterman_cpu(HMatrix, rowBacktracker, colBacktracker, refSeq,
				refLen, qrySeq, qryLen, match, mismatch, gapOpenPenalty,
				gapExtPenalty, alignScore, alignStart, alignEnd, alignLen,
				refConsensus, qryConsensus);

		/* Print results */
		if (outFormat == SAM_OUTPUT_FORMAT)
		{
			printOutputSam(qrsPerIteration, qryName, refName, alignScore,
					alignStart, alignEnd, alignLen, combinedLen, qryConsensus,
					refConsensus, outputFilePtr);
		}
		else
		{
			if ((i % 2) == 0)
				isFirstQry = 1;
			else
				isFirstQry = 0;
			printOutput_paired_cpu(qrsPerIteration, qryName, refName,
					alignScore, refDistance, alignStart, alignEnd, alignLen,
					combinedLen, qryConsensus, refConsensus, outputFilePtr,
					isFirstQry);
		}
	}
	fclose(tmpFilePtr);
	fclose(outputFilePtr);

	/* Free resources */
	fprintf(stderr, "Freeing host resources...\n");
	free(qryName);
	free(qrySeq);
	free(refName);
	free(refSeq);
	free(alignScore);
	free(alignStart);
	free(alignEnd);
	free(alignLen);
	free(refConsensus);
	free(qryConsensus);
	free(refDistance);

}


/**
 * Returns memory required on the GPU per reference-query pair.
 *
 * @param refLength Length of each reference sequence.
 * @param queryLength Length of each query sequence.
 * @return Memory required on the GPU per reference-sequence pair.
 */
int getMemPerRefQryPair(int refLength, int queryLength)
{
	int memHMatrix = 0;
	int memBacktrackMatrix = 0;
	int memRefSeq = (refLength + 1) * sizeof(char);
	int memQuerySeq = (queryLength + 1) * sizeof(char);
	int memScore = sizeof(float);
	int memAlignStart = sizeof(int);
	int memAlignEnd = sizeof(int);
	int memAlignLength = sizeof(ushort);
	int memRefConsensus = (refLength + queryLength + 2) * sizeof(char);
	int memQueryConsensus = (refLength + queryLength + 2) * sizeof(char);
	int total = memHMatrix + memBacktrackMatrix + memRefSeq + memQuerySeq
			+ memScore + memAlignStart + memAlignEnd + memAlignLength
			+ memRefConsensus + memQueryConsensus;

	return total;
}


/**
 * Returns the number of queries suitable for each iteration
 *
 * @param deviceProp Pointer to CUDA device properties
 * @param memPerRefQueryPair Memory required per reference query pair
 * @param numTotalQueries Total number of queries
 * @return Number of queries suitable for each iteration
 */
int getNumQueriesPerIteration(cudaDeviceProp *deviceProp,
		int memPerRefQueryPair, int numTotalQueries)
{
	long globalMemory = deviceProp->totalGlobalMem - GLOBAL_MEM_MARGIN;

	int numQueriesPerIteration = (int) floor(
			((float) globalMemory / memPerRefQueryPair));

	/* Make sure the number of queries per iteration is not greater
	 * than the max number of blocks */
	int maxBlocksX = deviceProp->maxGridSize[0];
	if (numQueriesPerIteration > maxBlocksX)
		numQueriesPerIteration = maxBlocksX;

	/* Make sure the number of queries per iteration is not greater
	 * than the total number of queries */
	if (numQueriesPerIteration > numTotalQueries)
		numQueriesPerIteration = numTotalQueries;

	/* Make sure there are even number of queries per iteration. */
//	if (numQueriesPerIteration % 2 == 1)
//		--numQueriesPerIteration;
	if (numQueriesPerIteration >= 2 && numQueriesPerIteration % 2 == 1)
		--numQueriesPerIteration;

	return numQueriesPerIteration;
}


/**
 * Returns the number of queries suitable for each iteration
 *
 * @param deviceProp Pointer to CUDA device properties
 * @param memPerRefQueryPair Memory required per reference query pair
 * @param numTotalQueries Total number of queries
 * @return Number of queries suitable for each iteration
 */
int getNumQueriesPerIteration_paired(cudaDeviceProp *deviceProp,
		int memPerRefQueryPair, int numTotalQueries)
{
	long globalMemory = deviceProp->totalGlobalMem - GLOBAL_MEM_MARGIN;

	int numQueriesPerIteration = (int) floor(
			((float) globalMemory / memPerRefQueryPair));

	/* Make sure the number of queries per iteration is not greater
	 * than the max number of blocks */
	int maxBlocksX = deviceProp->maxGridSize[0];
	if (numQueriesPerIteration > maxBlocksX)
		numQueriesPerIteration = maxBlocksX;

	/* Make sure the number of queries per iteration is not greater
	 * than the total number of queries */
	if (numQueriesPerIteration > numTotalQueries)
		numQueriesPerIteration = numTotalQueries;

	return numQueriesPerIteration;
}


/**
 * Returns number of sequences for the given iteration
 *
 * @param iterationIndex Index of the current iteration
 * @param numTotalSeqs Total number of sequences
 * @param numIterations Number of iterations
 * @param maxSeqsPerIteration Maximum number of sequences per iteration
 * @return Number of sequences for current iteration
 */
int getNumSeqsInCurrentIteration(
		int iterationIndex,
		int numTotalSeqs,
		int numIterations,
		int maxSeqsPerIteration)
{
	/* If this is the last iteration, then adjust the number of
	 * sequences for this iteration */
	if (iterationIndex == numIterations - 1)
	{
		return getNumSeqsInLastIteration(
				numTotalSeqs,
				numIterations,
				maxSeqsPerIteration);
	}
	/* Otherwise, just return the maximum number of sequences
	 * per iteration */
	else
		return maxSeqsPerIteration;
}


/**
 * Returns number of references in last iteration
 *
 * @param numTotalSeqs Total number of sequences
 * @param numIterations Number of iterations
 * @param maxSeqsPerIteration Maximum number of sequences per iteration
 * @return Number of sequences in last iteration
 */
int getNumSeqsInLastIteration(
		int numTotalSeqs,
		int numIterations,
		int maxSeqsPerIteration)
{
	int maxSeqs = numIterations * maxSeqsPerIteration;
    if (numTotalSeqs < maxSeqs)
    {
    	int numSeqsUntilLastIteration =
    			((numIterations - 1) * maxSeqsPerIteration);
    	return (numTotalSeqs - numSeqsUntilLastIteration);
    }
    else
    	return maxSeqsPerIteration;
}


/**
 * Aligns reference and query sequences
 *
 * @param querySeq_d Array containing query sequences
 * @param refSeq_d Array containing the reference sequences
 * @param querySeqLength Max query sequence length
 * @param refSeqLength Max reference sequence length
 * @param match Match value
 * @param mismatch Mismatch value
 * @param score_d Array where alignment scores will be stored
 * @param alignStart_d Array where alignment start positions will be stored
 * @param alignEnd_d Array where alignment end positions will be stored
 * @param alignLength_d Array where alignment lengths will be stored
 * @param queryConsensus_d Array where the aligned query sequences will be
 * stored
 * @param refConsensus_d Array where the aligned reference sequences will be
 * stored
 * @param queriesPerIteration Number of queries per iteration or partition
 * @param refsPerIteration Number of references per iteration or partition
 * @param maxThreadsPerBlock Max threads per block.
 */
void alignSequences(
		char *querySeq_d,
		char *refSeq_d,
		int querySeqLength,
		int refSeqLength,
		float match,
		float mismatch,
		float gapOpenPenalty,
		float gapExtPenalty,
		float *score_d,
		int *alignStart_d,
		int *alignEnd_d,
		int *alignLength_d,
		char *queryConsensus_d,
		char *refConsensus_d,
		int queriesPerIteration,
		int refsPerIteration,
		int maxThreadsPerBlock)
{
	dim3 dimGrid, dimBlock;
	dimBlock.x = maxThreadsPerBlock;
	dimBlock.y = 1;
	dimBlock.z = 1;
	dimGrid.x = queriesPerIteration;
	dimGrid.y = 1;
	smithWaterman<<<dimGrid, dimBlock>>>(
				refSeq_d,
				refSeqLength,
				querySeq_d,
				querySeqLength,
				match,
				mismatch,
				gapOpenPenalty,
				gapExtPenalty,
				score_d,
				alignStart_d,
				alignEnd_d,
				alignLength_d,
				refConsensus_d,
				queryConsensus_d,
				queriesPerIteration,
				refsPerIteration);
	PRINT_CUDA_ERROR()
}


/**
 * Returns the number of match records in the given file.
 *
 * @param filteredResultsFile File containing filtered matches.
 * @return Number of match records in the given file.
 */
static int getNumMatches(const char *filteredResultsFile)
{
	FILE *fp = fopen(filteredResultsFile, "r");
	char line[MAX_LINE_LENGTH];
	int count = 0;
	while (fgets(line, MAX_LINE_LENGTH, fp) != NULL)
		++count;
	fclose(fp);
	return count;
}


/**
 * Returns 1 if the query name has a '/2' in the end, otherwise returns 0.
 *
 * @param	qryName	Query name.
 * @return	Returns 1 if the given query name has a '/2' in the end, otherwise
 * returns 0.
 */
//static int isMate(const char *qryName, int qryNameLen)
//{
//
//}
