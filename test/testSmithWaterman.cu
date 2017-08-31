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
#include <check.h>
#include "../src/smithWaterman.h"
#include "testSmithWaterman.h"
#include "../src/common.h"


/**
 * Tests smithWaterman function.
 */
START_TEST(smithWaterman)
{
	/* One reference and one query. */
	{
		char *qrySeq_d, *refSeq_d, *qryConsensus_d, *refConsensus_d;
		float *score_d;
		int *alignStart_d, *alignEnd_d, *alignLen_d;
		int qryLen = 10, refLen = 20;
		int combinedLen = qryLen + refLen + 2;

		cudaMalloc(&qrySeq_d, qryLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&refSeq_d, refLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&qryConsensus_d, combinedLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&refConsensus_d, combinedLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&score_d, sizeof(float));
		PRINT_CUDA_ERROR()
		cudaMalloc(&alignStart_d, sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMalloc(&alignEnd_d, sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMalloc(&alignLen_d, sizeof(int));
		PRINT_CUDA_ERROR()

		char *qrySeq_h = (char *) malloc(qryLen * sizeof(char));
		strcpy(qrySeq_h, "ACGTACGTAA");
		cudaMemcpy(qrySeq_d, qrySeq_h, qryLen * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		char *refSeq_h = (char *) malloc(refLen * sizeof(char));
		strcpy(refSeq_h, "CCCCAAAACGTCCGTTACGT");
		cudaMemcpy(refSeq_d, refSeq_h, refLen * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		float match = 2.0f, mismatch = -1.0f;
		float gapOpenPenalty = -10.0f, gapExtPenalty = -1.0f;
		int numQrs = 1, numRefs = 1;
		smithWaterman<<<1, 512>>>(refSeq_d, refLen, qrySeq_d, qryLen,
				match, mismatch, gapOpenPenalty, gapExtPenalty, score_d,
				alignStart_d, alignEnd_d, alignLen_d, refConsensus_d,
				qryConsensus_d, numQrs, numRefs);

		char *qryConsensus_h = (char *) malloc(combinedLen * sizeof(char));
		char *refConsensus_h = (char *) malloc(combinedLen * sizeof(char));
		float *score_h = (float *) malloc(sizeof(float));
		int *alignStart_h = (int *) malloc(sizeof(int));
		int *alignEnd_h = (int *) malloc(sizeof(int));
		int *alignLen_h = (int *) malloc(sizeof(int));
		cudaMemcpy(score_h, score_d, sizeof(float), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(alignStart_h, alignStart_d, sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(alignEnd_h, alignEnd_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(alignLen_h, alignLen_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(alignEnd_h, alignEnd_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qryConsensus_h, qryConsensus_d, combinedLen * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(refConsensus_h, refConsensus_d, combinedLen * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (*score_h != 14.0f || *alignStart_h != 8 || *alignEnd_h != 17
				|| *alignLen_h != 10
				|| strncmp(qryConsensus_h, "AATGCATGCA", 10) != 0
				|| strncmp(refConsensus_h, "ATTGCCTGCA", 10) != 0)
			fail("Incorrect behavior when there is 1 query sequence and "
					"1 reference sequence.\n");

		cudaFree(qrySeq_d);
		cudaFree(refSeq_d);
		cudaFree(qryConsensus_d);
		cudaFree(refConsensus_d);
		cudaFree(score_d);
		cudaFree(alignStart_d);
		cudaFree(alignEnd_d);
		cudaFree(alignLen_d);
		free(score_h);
		free(alignStart_h);
		free(alignEnd_h);
		free(alignLen_h);
		free(qrySeq_h);
		free(refSeq_h);
		free(qryConsensus_h);
		free(refConsensus_h);
	}
}
END_TEST


/**
 * Tests computeMaxScore function.
 */
START_TEST(computeMaxScore)
{
	/* One query and one reference (no gaps). */
	{
		int qryLen = 10, refLen = 20;
		int combinedLen = qryLen + refLen + 2;
		float match = 2.0f, mismatch = -1.0f;
		float gapOpenPenalty = -10.0f, gapExtPenalty = -1.0f;
		char *topLeft_d, *topLeft_h, *left_d, *left_h, *top_d, *top_h;
		char *qrySeq_h, *qrySeq_d, *refSeq_h, *refSeq_d;
		float *maxScore_h, *maxScore_d;
		int *cellMax_h, *cellMax_d;

		int numCells = (qryLen + 1) * (refLen + 1);
		int charSize = (int) (sizeof(char) * NUM_BITS_IN_A_BYTE);
		int numBytes = (int) ceil((float) numCells / charSize);
		cudaMalloc(&topLeft_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&left_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&top_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&qrySeq_d, combinedLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&refSeq_d, combinedLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&maxScore_d, sizeof(float));
		PRINT_CUDA_ERROR()
		cudaMalloc(&cellMax_d, sizeof(int));
		PRINT_CUDA_ERROR()

		qrySeq_h = (char *) malloc(qryLen * sizeof(char));
		strcpy(qrySeq_h, "ACGTACGTAA");
		refSeq_h = (char *) malloc(refLen * sizeof(char));
		strcpy(refSeq_h, "CCCCAAAACGTCCGTTACGT");
		maxScore_h = (float *) malloc(sizeof(float));
		cellMax_h = (int *) malloc(sizeof(int));
		topLeft_h = (char *) malloc(numBytes * sizeof(char));
		top_h = (char *) malloc(numBytes * sizeof(char));
		left_h = (char *) malloc(numBytes * sizeof(char));

		cudaMemcpy(qrySeq_d, qrySeq_h, qryLen * sizeof(char),
				cudaMemcpyHostToDevice);
		cudaMemcpy(refSeq_d, refSeq_h, refLen * sizeof(char),
				cudaMemcpyHostToDevice);
		cudaMemset(topLeft_d, 0, numBytes * sizeof(char));
		cudaMemset(top_d, 0, numBytes * sizeof(char));
		cudaMemset(left_d, 0, numBytes * sizeof(char));
		computeMaxScore_wrapFunc<<<1, 512>>>(qryLen, refLen, topLeft_d,
				left_d, top_d, qrySeq_d, refSeq_d, match, mismatch,
				gapOpenPenalty, gapExtPenalty, maxScore_d, cellMax_d);
		cudaMemcpy(maxScore_h, maxScore_d, sizeof(float),
				cudaMemcpyDeviceToHost);
		cudaMemcpy(cellMax_h, cellMax_d, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(topLeft_h, top_d, numBytes * sizeof(char),
				cudaMemcpyDeviceToHost);
		cudaMemcpy(top_h, top_d, numBytes * sizeof(char),
				cudaMemcpyDeviceToHost);
		cudaMemcpy(left_h, left_d, numBytes * sizeof(char),
				cudaMemcpyDeviceToHost);

		if (*maxScore_h != 14.0f || *cellMax_h != 227
				|| (top_h[3] & (1 << 2)) != 0 || (left_h[3] & (1 << 2)) != 0
				|| (top_h[6] & (1 << 4)) != 0 || (left_h[6] & (1 << 4)) != 0
				|| (top_h[9] & (1 << 6)) != 0 || (left_h[9] & (1 << 6)) != 0
				|| (top_h[11] & (1 << 0)) != 0 || (left_h[11] & (1 << 0)) != 0
				|| (top_h[14] & (1 << 2)) != 0 || (left_h[14] & (1 << 2)) != 0
				|| (top_h[17] & (1 << 4)) != 0 || (left_h[17] & (1 << 4)) != 0
				|| (top_h[20] & (1 << 6)) != 0 || (left_h[20] & (1 << 6)) != 0
				|| (top_h[22] & (1 << 0)) != 0 || (left_h[22] & (1 << 0)) != 0
				|| (top_h[25] & (1 << 2)) != 0 || (left_h[25] & (1 << 2)) != 0
				|| (top_h[28] & (1 << 4)) != 0 || (left_h[28] & (1 << 4)) != 0)
			fail("Incorrect behavior when one query is aligned with one "
					"reference.\n");
		cudaFree(left_d);
		cudaFree(top_d);
		cudaFree(qrySeq_d);
		cudaFree(refSeq_d);
		cudaFree(maxScore_d);
		cudaFree(cellMax_d);
		free(qrySeq_h);
		free(refSeq_h);
		free(maxScore_h);
		free(cellMax_h);
	}

	/* One query and one reference (with gaps). */
	{
		int qryLen = 14, refLen = 13;
		int combinedLen = qryLen + refLen + 2;
		float match = 2.0f, mismatch = -1.0f;
		float gapOpenPenalty = -10.0f, gapExtPenalty = -1.0f;
		char *topLeft_d, *topLeft_h, *left_d, *left_h, *top_d, *top_h;
		char *qrySeq_h, *qrySeq_d, *refSeq_h, *refSeq_d;
		float *maxScore_h, *maxScore_d;
		int *cellMax_h, *cellMax_d;

		int numCells = (qryLen + 1) * (refLen + 1);
		int charSize = (int) (sizeof(char) * NUM_BITS_IN_A_BYTE);
		int numBytes = (int) ceil((float) numCells / charSize);
		cudaMalloc(&topLeft_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&left_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&top_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&qrySeq_d, combinedLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&refSeq_d, combinedLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&maxScore_d, sizeof(float));
		PRINT_CUDA_ERROR()
		cudaMalloc(&cellMax_d, sizeof(int));
		PRINT_CUDA_ERROR()

		qrySeq_h = (char *) malloc(qryLen * sizeof(char));
		strcpy(qrySeq_h, "CTTGTCACCATCCA");
		refSeq_h = (char *) malloc(refLen * sizeof(char));
		strcpy(refSeq_h, "CTTGTCACCATCA");
		maxScore_h = (float *) malloc(sizeof(float));
		cellMax_h = (int *) malloc(sizeof(int));
		topLeft_h = (char *) malloc(numBytes * sizeof(char));
		top_h = (char *) malloc(numBytes * sizeof(char));
		left_h = (char *) malloc(numBytes * sizeof(char));

		cudaMemcpy(qrySeq_d, qrySeq_h, qryLen * sizeof(char),
				cudaMemcpyHostToDevice);
		cudaMemcpy(refSeq_d, refSeq_h, refLen * sizeof(char),
				cudaMemcpyHostToDevice);
		cudaMemset(topLeft_d, 0, numBytes * sizeof(char));
		cudaMemset(top_d, 0, numBytes * sizeof(char));
		cudaMemset(left_d, 0, numBytes * sizeof(char));
		computeMaxScore_wrapFunc<<<1, 512>>>(qryLen, refLen, topLeft_d,
				left_d, top_d, qrySeq_d, refSeq_d, match, mismatch,
				gapOpenPenalty, gapExtPenalty, maxScore_d, cellMax_d);
		cudaMemcpy(maxScore_h, maxScore_d, sizeof(float),
				cudaMemcpyDeviceToHost);
		cudaMemcpy(cellMax_h, cellMax_d, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(topLeft_h, top_d, numBytes * sizeof(char),
				cudaMemcpyDeviceToHost);
		cudaMemcpy(top_h, top_d, numBytes * sizeof(char),
				cudaMemcpyDeviceToHost);
		cudaMemcpy(left_h, left_d, numBytes * sizeof(char),
				cudaMemcpyDeviceToHost);

		if (*maxScore_h != 24.0f || *cellMax_h != 180
				|| (top_h[1] & (1 << 0)) != 0 || (left_h[1] & (1 << 0)) != 0
				|| (top_h[3] & (1 << 1)) != 0 || (left_h[3] & (1 << 1)) != 0
				|| (top_h[5] & (1 << 2)) != 0 || (left_h[5] & (1 << 2)) != 0
				|| (top_h[7] & (1 << 3)) != 0 || (left_h[7] & (1 << 3)) != 0
				|| (top_h[9] & (1 << 4)) != 0 || (left_h[9] & (1 << 4)) != 0
				|| (top_h[11] & (1 << 5)) != 0 || (left_h[11] & (1 << 5)) != 0
				|| (top_h[13] & (1 << 6)) != 0 || (left_h[13] & (1 << 6)) != 0
				|| (top_h[15] & (1 << 7)) != 0 || (left_h[15] & (1 << 7)) != 0
				|| (top_h[16] & (1 << 0)) != 0 || (left_h[16] & (1 << 0)) != 0
				|| (top_h[18] & (1 << 1)) != 0 || (left_h[18] & (1 << 1)) != 0
				|| (top_h[20] & (1 << 2)) != 0 || (left_h[20] & (1 << 2)) != 0
				|| (top_h[22] & (1 << 3)) != 0 || (left_h[22] & (1 << 3)) != 0
				)
			fail("Incorrect behavior when one query is aligned with one "
					"reference (with gaps).\n");
		cudaFree(left_d);
		cudaFree(top_d);
		cudaFree(qrySeq_d);
		cudaFree(refSeq_d);
		cudaFree(maxScore_d);
		cudaFree(cellMax_d);
		free(qrySeq_h);
		free(refSeq_h);
		free(maxScore_h);
		free(cellMax_h);
	}


	/* Reference length is more than number of bases in reference. */
	{
		int qryLen = 10, refLen = 25;
		int combinedLen = qryLen + refLen + 2;
		float match = 2.0f, mismatch = -1.0f;
		float gapOpenPenalty = -10.0f, gapExtPenalty = -1.0f;
		char *topLeft_d, *topLeft_h, *left_d, *left_h, *top_d, *top_h;
		char *qrySeq_h, *qrySeq_d, *refSeq_h, *refSeq_d;
		float *maxScore_h, *maxScore_d;
		int *cellMax_h, *cellMax_d;

		int numCells = (qryLen + 1) * (refLen + 1);
		int charSize = (int) (sizeof(char) * NUM_BITS_IN_A_BYTE);
		int numBytes = (int) ceil((float) numCells / charSize);
		cudaMalloc(&topLeft_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&left_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&top_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&qrySeq_d, combinedLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&refSeq_d, combinedLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&maxScore_d, sizeof(float));
		PRINT_CUDA_ERROR()
		cudaMalloc(&cellMax_d, sizeof(int));
		PRINT_CUDA_ERROR()

		qrySeq_h = (char *) malloc(qryLen * sizeof(char));
		strcpy(qrySeq_h, "ACGTACGTAA");
		refSeq_h = (char *) malloc(refLen * sizeof(char));
		strcpy(refSeq_h, "CCCCAAAACGTCCGTTACGT");
		maxScore_h = (float *) malloc(sizeof(float));
		cellMax_h = (int *) malloc(sizeof(int));
		topLeft_h = (char *) malloc(numBytes * sizeof(char));
		top_h = (char *) malloc(numBytes * sizeof(char));
		left_h = (char *) malloc(numBytes * sizeof(char));

		cudaMemcpy(qrySeq_d, qrySeq_h, qryLen * sizeof(char),
				cudaMemcpyHostToDevice);
		cudaMemcpy(refSeq_d, refSeq_h, refLen * sizeof(char),
				cudaMemcpyHostToDevice);
		cudaMemset(topLeft_d, 0, numBytes * sizeof(char));
		cudaMemset(top_d, 0, numBytes * sizeof(char));
		cudaMemset(left_d, 0, numBytes * sizeof(char));
		computeMaxScore_wrapFunc<<<1, 512>>>(qryLen, refLen, topLeft_d,
				left_d, top_d, qrySeq_d, refSeq_d, match, mismatch,
				gapOpenPenalty, gapExtPenalty, maxScore_d, cellMax_d);
		cudaMemcpy(maxScore_h, maxScore_d, sizeof(float),
				cudaMemcpyDeviceToHost);
		cudaMemcpy(cellMax_h, cellMax_d, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(topLeft_h, top_d, numBytes * sizeof(char),
				cudaMemcpyDeviceToHost);
		cudaMemcpy(top_h, top_d, numBytes * sizeof(char),
				cudaMemcpyDeviceToHost);
		cudaMemcpy(left_h, left_d, numBytes * sizeof(char),
				cudaMemcpyDeviceToHost);

		if (*maxScore_h != 14.0f || *cellMax_h != 277
				|| (top_h[3] & (1 << 2)) != 0 || (left_h[3] & (1 << 2)) != 0
				|| (top_h[6] & (1 << 4)) != 0 || (left_h[6] & (1 << 4)) != 0
				|| (top_h[9] & (1 << 6)) != 0 || (left_h[9] & (1 << 6)) != 0
				|| (top_h[11] & (1 << 0)) != 0 || (left_h[11] & (1 << 0)) != 0
				|| (top_h[14] & (1 << 2)) != 0 || (left_h[14] & (1 << 2)) != 0
				|| (top_h[17] & (1 << 4)) != 0 || (left_h[17] & (1 << 4)) != 0
				|| (top_h[20] & (1 << 6)) != 0 || (left_h[20] & (1 << 6)) != 0
				|| (top_h[22] & (1 << 0)) != 0 || (left_h[22] & (1 << 0)) != 0
				|| (top_h[25] & (1 << 2)) != 0 || (left_h[25] & (1 << 2)) != 0
				|| (top_h[28] & (1 << 4)) != 0 || (left_h[28] & (1 << 4)) != 0)
			fail("Incorrect behavior when reference length is more than "
					"number of bases in reference.\n");
		cudaFree(left_d);
		cudaFree(top_d);
		cudaFree(qrySeq_d);
		cudaFree(refSeq_d);
		cudaFree(maxScore_d);
		cudaFree(cellMax_d);
		free(qrySeq_h);
		free(refSeq_h);
		free(maxScore_h);
		free(cellMax_h);
	}
}
END_TEST


/**
 * Tests the backtrack function.
 */
START_TEST(backtrack)
{
	/* One query aligned with one reference. */
	{
		char *qrySeq_d, *refSeq_d, *qryConsensus_d, *refConsensus_d;
		char *topLeft_d, *top_d, *left_d;
		float *score_d;
		int *alignStart_d, *alignEnd_d, *alignLen_d, *cellMax_d;
		int qryLen = 10, refLen = 20;
		int combinedLen = qryLen + refLen + 2;

		int numCells = (qryLen + 1) * (refLen + 1);
		int charSize = (int) sizeof(char);
		int numBytes = (int) ceil((float) numCells / charSize);
		cudaMalloc(&topLeft_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&left_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&top_d, numBytes * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&qrySeq_d, qryLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&refSeq_d, refLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&qryConsensus_d, combinedLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&refConsensus_d, combinedLen * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&score_d, sizeof(float));
		PRINT_CUDA_ERROR()
		cudaMalloc(&alignStart_d, sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMalloc(&alignEnd_d, sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMalloc(&alignLen_d, sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMalloc(&cellMax_d, sizeof(int));
		PRINT_CUDA_ERROR()

		char *qrySeq_h = (char *) malloc(qryLen * sizeof(char));
		strcpy(qrySeq_h, "ACGTACGTAA");
		cudaMemcpy(qrySeq_d, qrySeq_h, qryLen * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		char *refSeq_h = (char *) malloc(refLen * sizeof(char));
		strcpy(refSeq_h, "CCCCAAAACGTCCGTTACGT");
		cudaMemcpy(refSeq_d, refSeq_h, refLen * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		float match = 2.0f, mismatch = -1.0f;
		float gapOpenPenalty = -10.0f, gapExtPenalty = -1.0f;
		char *qryConsensus_h = (char *) malloc(combinedLen * sizeof(char));
		char *refConsensus_h = (char *) malloc(combinedLen * sizeof(char));
		float *score_h = (float *) malloc(sizeof(float));
		int *alignStart_h = (int *) malloc(sizeof(int));
		int *alignEnd_h = (int *) malloc(sizeof(int));
		int *alignLen_h = (int *) malloc(sizeof(int));
		int *cellMax_h = (int *) malloc(sizeof(int));

		cudaMemcpy(qrySeq_d, qrySeq_h, qryLen * sizeof(char),
				cudaMemcpyHostToDevice);
		cudaMemcpy(refSeq_d, refSeq_h, refLen * sizeof(char),
				cudaMemcpyHostToDevice);
		cudaMemset(topLeft_d, 0, numBytes * sizeof(char));
		cudaMemset(top_d, 0, numBytes * sizeof(char));
		cudaMemset(left_d, 0, numBytes * sizeof(char));
		computeMaxScore_wrapFunc<<<1, 512>>>(qryLen, refLen, topLeft_d,
				left_d, top_d, qrySeq_d, refSeq_d, match, mismatch,
				gapOpenPenalty, gapExtPenalty, score_d, cellMax_d);
		cudaMemcpy(score_h, score_d, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(cellMax_h, cellMax_d, sizeof(int), cudaMemcpyDeviceToHost);
		backtrack_wrapFunc<<<1, 1>>>(qryConsensus_d, refConsensus_d,
				alignStart_d, alignEnd_d, alignLen_d, topLeft_d, left_d,
				top_d, qryLen, refLen, *cellMax_h, qrySeq_d, refSeq_d);
		cudaMemcpy(qryConsensus_h, qryConsensus_d, combinedLen * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(refConsensus_h, refConsensus_d, combinedLen * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(alignStart_h, alignStart_d, sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(alignEnd_h, alignEnd_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(alignLen_h, alignLen_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()
		cudaMemcpy(score_h, score_d, sizeof(float), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (*score_h != 14.0f
				|| *alignStart_h != 8
				|| *alignEnd_h != 17
				|| *alignLen_h != 10
				|| strncmp(qryConsensus_h, "AATGCATGCA", 10) != 0
				|| strncmp(refConsensus_h, "ATTGCCTGCA", 10) != 0
				)
			fail("Incorrect behavior when there is 1 query sequence and "
					"1 reference sequence.\n");

		cudaFree(qrySeq_d);
		cudaFree(refSeq_d);
		cudaFree(qryConsensus_d);
		cudaFree(refConsensus_d);
		cudaFree(score_d);
		cudaFree(alignStart_d);
		cudaFree(alignEnd_d);
		cudaFree(alignLen_d);
		PRINT_CUDA_ERROR()
		free(score_h);
		free(alignStart_h);
		free(alignEnd_h);
		free(alignLen_h);
		free(qrySeq_h);
		free(refSeq_h);
		free(qryConsensus_h);
		free(refConsensus_h);
	}
}
END_TEST


/**
 * Tests the CPU version of Smith-Waterman.
 */
START_TEST(smithWaterman_cpu)
{
	int qryLen = 10, refLen = 20;
	int combinedLen = qryLen + refLen + 2;
	float match = 2.0f, mismatch = -1.0f;
	float gapOpenPenalty = -10.0f, gapExtPenalty = -1.0f;
	float *alignScore = (float *) malloc(sizeof(float));
	int *alignStart = (int *) malloc(sizeof(int));
	int *alignEnd = (int *) malloc(sizeof(int));
	int *alignLen = (int *) malloc(sizeof(int));

	char *qrySeq = (char *) malloc(qryLen * sizeof(char));
	strcpy(qrySeq, "ACGTACGTAA");
	char *refSeq = (char *) malloc(refLen * sizeof(char));
	strcpy(refSeq, "CCCCAAAACGTCCGTTACGT");
	char *qryConsensus = (char *) malloc(combinedLen * sizeof(char));
	char *refConsensus = (char *) malloc(combinedLen * sizeof(char));
	float *HMatrix = (float *) malloc((refLen + 1) * (qryLen + 1)
			* sizeof(float));
	int *rowBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
			* sizeof(int));
	int *colBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
			* sizeof(int));

	smithWaterman_cpu(HMatrix, rowBacktracker, colBacktracker, refSeq,
				refLen, qrySeq, qryLen, match, mismatch, gapOpenPenalty,
				gapExtPenalty, alignScore, alignStart, alignEnd, alignLen,
				refConsensus, qryConsensus);

	if (*alignScore != 14.0f || *alignStart != 8 || *alignEnd != 17
			|| *alignLen != 10
			|| strncmp(qryConsensus, "AATGCATGCA", 10) != 0
			|| strncmp(refConsensus, "ATTGCCTGCA", 10) != 0)
		fail("Incorrect behavior when there is 1 query sequence and "
				"1 reference sequence.\n");

	free(qrySeq);
	free(refSeq);
	free(qryConsensus);
	free(refConsensus);
	free(HMatrix);
	free(rowBacktracker);
	free(colBacktracker);
	free(alignScore);
	free(alignStart);
	free(alignEnd);
	free(alignLen);
}
END_TEST


/**
 * Tests the CPU-based 'computeMaxScore' implementation.
 */
START_TEST(computeMaxScore_cpu)
{
	/* No gaps. */
	{
		int qryLen = 10, refLen = 20;
		float match = 2.0f, mismatch = -1.0f;
		float gapOpenPenalty = -10.0f, gapExtPenalty = -1.0f;
		char *qrySeq, *refSeq;
		float *alignScore = (float *) malloc(sizeof(float));
		int *maxRow = (int *) malloc(sizeof(int));
		int *maxCol = (int *) malloc(sizeof(int));
		float *H = (float *) malloc((refLen + 1) * (qryLen + 1) * sizeof(float));
		int *rowBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
				* sizeof(int));
		int *colBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
				* sizeof(int));

		qrySeq = (char *) malloc(qryLen * sizeof(char));
		strcpy(qrySeq, "ACGTACGTAA");
		refSeq = (char *) malloc(refLen * sizeof(char));
		strcpy(refSeq, "CCCCAAAACGTCCGTTACGT");

		computeMaxScore_cpu_wrapFunc(H, rowBacktracker, colBacktracker, refSeq,
					refLen, qrySeq, qryLen, match, mismatch, gapOpenPenalty,
					gapExtPenalty, alignScore, maxRow, maxCol);

		if (*alignScore != 14.0f && *maxRow != 10 && *maxCol != 17)
			fail("Incorrect behavior when one query is aligned with one "
					"reference.\n");
		free(qrySeq);
		free(refSeq);
		free(alignScore);
		free(maxRow);
		free(maxCol);
		free(H);
		free(rowBacktracker);
		free(colBacktracker);
	}

	/* With gaps. */
	{
		int qryLen = 14, refLen = 13;
		float match = 2.0f, mismatch = -1.0f;
		float gapOpenPenalty = -10.0f, gapExtPenalty = -1.0f;
		char *qrySeq, *refSeq;
		float *alignScore = (float *) malloc(sizeof(float));
		int *maxRow = (int *) malloc(sizeof(int));
		int *maxCol = (int *) malloc(sizeof(int));
		float *H = (float *) malloc((refLen + 1) * (qryLen + 1) * sizeof(float));
		int *rowBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
				* sizeof(int));
		int *colBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
				* sizeof(int));

		qrySeq = (char *) malloc(qryLen * sizeof(char));
		strcpy(qrySeq, "CTTGTCACCATCCA");
		refSeq = (char *) malloc(refLen * sizeof(char));
		strcpy(refSeq, "CTTGTCACCATCA");

		int i, size = (refLen + 1) * (qryLen + 1);
		for (i = 0; i < size; ++i)
			H[i] = 0.0f;

		computeMaxScore_cpu_wrapFunc(H, rowBacktracker, colBacktracker, refSeq,
					refLen, qrySeq, qryLen, match, mismatch, gapOpenPenalty,
					gapExtPenalty, alignScore, maxRow, maxCol);

		if (*alignScore != 24.0f && *maxRow != 12 && *maxCol != 12)
			fail("Incorrect behavior when alignment has gaps.\n");
		free(qrySeq);
		free(refSeq);
		free(alignScore);
		free(maxRow);
		free(maxCol);
		free(H);
		free(rowBacktracker);
		free(colBacktracker);
	}

	/* Reference length is more than number of bases in reference. */
	{
		int qryLen = 10, refLen = 25;
		float match = 2.0f, mismatch = -1.0f;
		float gapOpenPenalty = -10.0f, gapExtPenalty = -1.0f;
		char *qrySeq, *refSeq;
		float *alignScore = (float *) malloc(sizeof(float));
		int *maxRow = (int *) malloc(sizeof(int));
		int *maxCol = (int *) malloc(sizeof(int));
		float *H = (float *) malloc((refLen + 1) * (qryLen + 1) * sizeof(float));
		int *rowBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
				* sizeof(int));
		int *colBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
				* sizeof(int));

		qrySeq = (char *) malloc((qryLen + 1) * sizeof(char));
		strcpy(qrySeq, "ACGTACGTAA");
		refSeq = (char *) malloc(refLen * sizeof(char));
		strcpy(refSeq, "CCCCAAAACGTCCGTTACGT");

		computeMaxScore_cpu_wrapFunc(H, rowBacktracker, colBacktracker, refSeq,
					refLen, qrySeq, qryLen, match, mismatch, gapOpenPenalty,
					gapExtPenalty, alignScore, maxRow, maxCol);

		if (*alignScore != 14.0f && *maxRow != 10 && *maxCol != 17)
			fail("Incorrect behavior when reference length is more than "
					"number of bases in reference.\n");
		free(qrySeq);
		free(refSeq);
		free(alignScore);
		free(maxRow);
		free(maxCol);
		free(H);
		free(rowBacktracker);
		free(colBacktracker);
	}
}
END_TEST


/**
 * Tests the CPU implementation of @a backtrack function.
 */
START_TEST(backtrack_cpu)
{
	int qryLen = 10, refLen = 20;
	int combinedLen = qryLen + refLen + 2;

	char *qrySeq = (char *) malloc(qryLen * sizeof(char));
	strcpy(qrySeq, "ACGTACGTAA");
	char *refSeq = (char *) malloc(refLen * sizeof(char));
	strcpy(refSeq, "CCCCAAAACGTCCGTTACGT");

	float match = 2.0f, mismatch = -1.0f;
	float gapOpenPenalty = -10.0f, gapExtPenalty = -1.0f;
	char *qryConsensus = (char *) malloc(combinedLen * sizeof(char));
	char *refConsensus = (char *) malloc(combinedLen * sizeof(char));
	float *alignScore = (float *) malloc(sizeof(float));
	int *alignStart = (int *) malloc(sizeof(int));
	int *alignEnd = (int *) malloc(sizeof(int));
	int *alignLen = (int *) malloc(sizeof(int));
	int *maxRow = (int *) malloc(sizeof(int));
	int *maxCol = (int *) malloc(sizeof(int));
	float *H = (float *) malloc((refLen + 1) * (qryLen + 1) * sizeof(float));
	int *rowBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
			* sizeof(int));
	int *colBacktracker = (int *) malloc((refLen + 1) * (qryLen + 1)
			* sizeof(int));

	computeMaxScore_cpu_wrapFunc(H, rowBacktracker, colBacktracker, refSeq,
				refLen, qrySeq, qryLen, match, mismatch, gapOpenPenalty,
				gapExtPenalty, alignScore, maxRow, maxCol);
	backtrack_cpu_wrapFunc(rowBacktracker, colBacktracker, *maxRow, *maxCol,
				refSeq, refLen, qrySeq, qryLen, alignStart, alignEnd, alignLen,
				refConsensus, qryConsensus);

	if (*alignScore != 14.0f || *alignStart != 8 || *alignEnd != 17
			|| *alignLen != 10
			|| strncmp(qryConsensus, "AATGCATGCA", 10) != 0
			|| strncmp(refConsensus, "ATTGCCTGCA", 10) != 0)
		fail("Incorrect behavior when there is 1 query sequence and "
				"1 reference sequence.\n");

	free(alignScore);
	free(alignStart);
	free(alignEnd);
	free(alignLen);
	free(qrySeq);
	free(refSeq);
	free(qryConsensus);
	free(refConsensus);
	free(maxRow);
	free(maxCol);
}
END_TEST


/**
 * Tests @a isMatch function.
 */
START_TEST(isMatch)
{
	/* Both characters are same. */
	{
		char a = 'A';
		char b = 'A';
		int isMatch;
		char *a_d, *b_d;
		int *isMatch_d;

		cudaMalloc(&a_d, sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&b_d, sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&isMatch_d, sizeof(int));
		PRINT_CUDA_ERROR()

		cudaMemcpy(a_d, &a, sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		cudaMemcpy(b_d, &b, sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		isMatch_wrap<<<1, 1>>>(a_d, b_d, isMatch_d);

		cudaMemcpy(&isMatch, isMatch_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (isMatch != 0)
			fail("Incorrect behavior when both characters are same.\n");

		cudaFree(a_d);
		cudaFree(b_d);
		cudaFree(isMatch_d);
	}

	/* Both characters are different and none of them is 'N'. */
	{
		char a = 'A';
		char b = 'C';
		int isMatch;
		char *a_d, *b_d;
		int *isMatch_d;

		cudaMalloc(&a_d, sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&b_d, sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&isMatch_d, sizeof(int));
		PRINT_CUDA_ERROR()

		cudaMemcpy(a_d, &a, sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		cudaMemcpy(b_d, &b, sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		isMatch_wrap<<<1, 1>>>(a_d, b_d, isMatch_d);

		cudaMemcpy(&isMatch, isMatch_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (isMatch != -1)
			fail("Incorrect behavior when both characters are different "
					"and none of them is 'N'.\n");

		cudaFree(a_d);
		cudaFree(b_d);
		cudaFree(isMatch_d);
	}

	/* Both characters are 'N'. */
	{
		char a = 'N';
		char b = 'N';
		int isMatch;
		char *a_d, *b_d;
		int *isMatch_d;

		cudaMalloc(&a_d, sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&b_d, sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&isMatch_d, sizeof(int));
		PRINT_CUDA_ERROR()

		cudaMemcpy(a_d, &a, sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		cudaMemcpy(b_d, &b, sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		isMatch_wrap<<<1, 1>>>(a_d, b_d, isMatch_d);

		cudaMemcpy(&isMatch, isMatch_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (isMatch != 0)
			fail("Incorrect behavior when both characters are 'N'.\n");

		cudaFree(a_d);
		cudaFree(b_d);
		cudaFree(isMatch_d);
	}

	/* Both are different and first character is 'N'. */
	{
		char a = 'N';
		char b = 'A';
		int isMatch;
		char *a_d, *b_d;
		int *isMatch_d;

		cudaMalloc(&a_d, sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&b_d, sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&isMatch_d, sizeof(int));
		PRINT_CUDA_ERROR()

		cudaMemcpy(a_d, &a, sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		cudaMemcpy(b_d, &b, sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		isMatch_wrap<<<1, 1>>>(a_d, b_d, isMatch_d);

		cudaMemcpy(&isMatch, isMatch_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (isMatch != 0)
			fail("Incorrect behavior when both characters are different "
					"and first character is 'N'.\n");

		cudaFree(a_d);
		cudaFree(b_d);
		cudaFree(isMatch_d);
	}

	/* Both are different and second character is 'N'. */
	{
		char a = 'A';
		char b = 'N';
		int isMatch;
		char *a_d, *b_d;
		int *isMatch_d;

		cudaMalloc(&a_d, sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&b_d, sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMalloc(&isMatch_d, sizeof(int));
		PRINT_CUDA_ERROR()

		cudaMemcpy(a_d, &a, sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		cudaMemcpy(b_d, &b, sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		isMatch_wrap<<<1, 1>>>(a_d, b_d, isMatch_d);

		cudaMemcpy(&isMatch, isMatch_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (isMatch != 0)
			fail("Incorrect behavior when both characters are different "
					"and second character is 'N'.\n");

		cudaFree(a_d);
		cudaFree(b_d);
		cudaFree(isMatch_d);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *smithWatermanSuite(void)
{
	Suite *s = suite_create("smithWaterman");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, smithWaterman);
	tcase_add_test(testCaseCore, computeMaxScore);
	tcase_add_test(testCaseCore, backtrack);
//	tcase_add_test(testCaseCore, smithWaterman_cpu);
	tcase_add_test(testCaseCore, computeMaxScore_cpu);
	tcase_add_test(testCaseCore, backtrack_cpu);
	tcase_add_test(testCaseCore, isMatch);
	suite_add_tcase (s, testCaseCore);

	return s;
}
