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

#include <stdlib.h>
#include <stdio.h>
#include "memory.h"
#include "common.h"


/**
 * Copies memory from device to host
 *
 * @param score Array on host where alignment scores will be copied to
 * @param score_d Array on device where alignment scores will be copied from
 * @param alignStart Array on host where alignment start positions will be
 * copied to
 * @param alignStart_d Array on device where alignment start positions will
 * be copied from
 * @param alignEnd Array on host where alignment end positions will be copied to
 * @param alignEnd_d Array on device where alignment end positions will be
 * copied from
 * @param alignLength Array on host where alignment lengths will be copied to
 * @param alignLength_d Array on device where alignment lengths will be copied
 * from
 * @param queryConsensus Array on host where aligned query sequences will be
 * copied to
 * @param queryConsensus_d Array on device where aligned query sequences will
 * be copied from
 * @param refConsensus Array on host where aligned reference sequences will
 * be copied to
 * @param refConsensus_d Array on device where aligned reference sequences
 * will be copied from
 * @param numQueries Number of queries
 * @param numRefs Number of references
 * @param combinedLength Combined length of reference and query sequence
 */
void copyMemDeviceToHost(
		float *score,
		float *score_d,
		int *alignStart,
		int *alignStart_d,
		int *alignEnd,
		int *alignEnd_d,
		int *alignLength,
		int *alignLength_d,
		char *queryConsensus,
		char *queryConsensus_d,
		char *refConsensus,
		char *refConsensus_d,
		int numQueries,
		int numRefs,
		int combinedLength)
{
	cudaMemcpy(score, score_d, numQueries * sizeof(float),
			cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	cudaMemcpy(alignStart, alignStart_d, numQueries * sizeof(int),
			cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	cudaMemcpy(alignEnd, alignEnd_d, numQueries * sizeof(int),
			cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	cudaMemcpy(alignLength, alignLength_d, numQueries * sizeof(int),
			cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	cudaMemcpy(refConsensus, refConsensus_d,
			numQueries * combinedLength * sizeof(char), cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	cudaMemcpy(queryConsensus, queryConsensus_d,
			numQueries * combinedLength * sizeof(char), cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()
}


/**
 * Allocates memory on device.
 *
 * @param querySeq Array that will store the query sequences.
 * @param refSeq Array that will store the reference sequences.
 * @param score Array that will store the alignment scores.
 * @param hMatrix Matrix that will be used to compute alignment scores.
 * @param cellBacktracker Matrix that will be used for backtracking.
 * @param alignStart Array that will store alignment start positions on the
 * reference.
 * @param alignEnd Array that will store alignment end positions on the
 * reference.
 * @param alignLength Array that will store alignment lengths.
 * @param queryConsensus Array that will store aligned query sequences.
 * @param refConsensus Array that will store aligned reference sequences.
 * @param numQueries Number of query sequences that will be stored in querySeq.
 * This is also the number of reference sequences that will be stored in
 * refSeq.
 * @param querySeqLength Max query sequence length.
 * @param refSeqLength Max reference sequence length.
 * @param combinedLength Combined length of reference and query sequence.
 */
void allocateDeviceMem(char **querySeq, char **refSeq, float **score,
		float **hMatrix, int **cellBacktracker, int **alignStart,
		int **alignEnd, int **alignLength, char **queryConsensus,
		char **refConsensus, int numQueries, int querySeqLength,
		int refSeqLength, int combinedLength)
{
	cudaMalloc((void **) querySeq, numQueries * (querySeqLength + 1) * sizeof(char));
	PRINT_CUDA_ERROR()

	cudaMalloc((void **) refSeq, numQueries * (refSeqLength + 1) * sizeof(char));
	PRINT_CUDA_ERROR()

	cudaMalloc((void **) score, numQueries * sizeof(float));
	PRINT_CUDA_ERROR()

    cudaMalloc((void **) alignStart, numQueries * sizeof(int));
    PRINT_CUDA_ERROR()

	cudaMalloc((void **) alignEnd, numQueries * sizeof(int));
	PRINT_CUDA_ERROR()

	cudaMalloc((void **) alignLength, numQueries * sizeof(int));
	PRINT_CUDA_ERROR()

	cudaMalloc((void **) refConsensus,
			numQueries * combinedLength * sizeof(char));
	PRINT_CUDA_ERROR()

	cudaMalloc((void **) queryConsensus,
			numQueries * combinedLength * sizeof(char));
	PRINT_CUDA_ERROR()
}


/**
 * Allocates memory on host.
 *
 * @param score Array for alignment scores.
 * @param refDistance Array for reference distances.
 * @param alignStart Array for alignment start positions.
 * @param alignEnd Array for alignment end positions.
 * @param alignLength Array for alignment lengths.
 * @param queryConsensus Array for aligned query sequences.
 * @param refConsensus Array for alignmed reference sequences.
 * @param queryName Array for query names.
 * @param querySeq Array for query sequences.
 * @param refName Array for reference names.
 * @param refSeq Array for reference sequences.
 * @param numQueries Number of queries that will be stored in querySeq.
 * @param maxQryNameLength Maximum query name length.
 * @param maxRefNameLength Maximum reference name length.
 * @param querySeqLength Max length of query sequence.
 * @param refSeqLength Max length of reference sequence.
 * @param combinedLength Combined length of reference and query sequence.
 */
void allocateHostMem(float **score, int **refDistance, int **alignStart,
		int **alignEnd, int **alignLength, char **queryConsensus,
		char **refConsensus, char **queryName, char **querySeq, char **refName,
		char **refSeq, int numQueries, int maxQryNameLength,
		int maxRefNameLength, int querySeqLength, int refSeqLength,
		int combinedLength)
{
	*score = (float *) malloc(numQueries * sizeof(float));
	*refDistance = (int *) malloc(numQueries * sizeof(int));
	*alignStart = (int *) malloc(numQueries * sizeof(int));
	*alignEnd = (int *) malloc(numQueries * sizeof(int));
	*alignLength = (int *) malloc(numQueries * sizeof(int));
	*refConsensus = (char *) malloc(numQueries * combinedLength * sizeof(char));
	*queryConsensus = (char *) malloc(numQueries * combinedLength * sizeof(char));
	*refName = (char *)  malloc(numQueries * maxRefNameLength * sizeof(char));
	*refSeq = (char *)  malloc(numQueries * (refSeqLength + 1) * sizeof(char));
	*queryName = (char *) malloc(numQueries * maxQryNameLength * sizeof(char));
	*querySeq = (char *) malloc(numQueries * (querySeqLength + 1) * sizeof(char));
}


/**
 * Copies memory from host to device.
 *
 * @param seq_d Character array on the device where sequences will be copied to.
 * @param seq Character array on the host where sequences will be copied from.
 * @param numBytes Number of bytes to be copied from host to device.
 */
void copyMemHostToDevice(char *seq_d, char *seq, int numBytes)
{
	cudaMemcpy(seq_d, seq, numBytes, cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()
}

