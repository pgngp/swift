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
#include "../src/align.h"
#include "testAlign.h"
#include "../src/common.h"


/**
 * Tests alignSequences function.
 */
START_TEST(alignSequences)
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
		int numQrs = 1, numRefs = 1, maxThreadsPerBlock = 512;
		alignSequences(qrySeq_d, refSeq_d, qryLen, refLen, match, mismatch,
				gapOpenPenalty, gapExtPenalty, score_d, alignStart_d,
				alignEnd_d, alignLen_d, qryConsensus_d, refConsensus_d,
				numQrs, numRefs, maxThreadsPerBlock);

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
 * Tests @a chooseBestAlignment function.
 */
START_TEST(chooseBestAlignment)
{
	/* There is 1 query and size is 1. */
	{
		int size = 1;
		float score[] = {5.0};
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		sprintf(qryName, "qry0\0");

		chooseBestAlignment_wrap(qryName, qryNameLen, score, size);

		if (score[0] != 5.0)
			fail("Incorrect behavior when there is 1 query and size is 1.\n");

		free(qryName);
	}

	/* Size is 5, all belonging to the same query, and there is only 1 max
	 * score. */
	{
		int size = 5;
		float score[] = {5.0, 8.0, 1.0, 15.0, 5.0};
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		int i, tmp = 0;
		for (i = 0; i < size; ++i)
		{
			sprintf(qryName + tmp, "qry0\0");
			tmp += qryNameLen;
		}

		chooseBestAlignment_wrap(qryName, qryNameLen, score, size);

		if (score[0] != -1.0 || score[1] != -1.0 || score[2] != -1.0
				|| score[3] != 15.0 || score[4] != -1.0)
			fail("Incorrect behavior when size is 5, all belonging to the "
					"same query, and there is only 1 max score.\n");

		free(qryName);
	}

	/* Size is 5, all belonging to the same query, and there are 2 max
	 * scores. */
	{
		int size = 5;
		float score[] = {5.0, 15.0, 1.0, 15.0, 5.0};
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		int i, tmp = 0;
		for (i = 0; i < size; ++i)
		{
			sprintf(qryName + tmp, "qry0\0");
			tmp += qryNameLen;
		}

		chooseBestAlignment_wrap(qryName, qryNameLen, score, size);

		if (score[0] == -1.0 && score[2] == -1.0 && score[4] == -1.0
				&& ((score[1] == 15.0 && score[3] == -1.0)
						|| (score[1] == -1.0 && score[3] == 15.0)))
		{ }
		else
			fail("Incorrect behavior when size is 5, all belonging to the "
					"same query, and there are 2 max scores.\n");

		free(qryName);
	}

	/* There are 2 queries, size is 6. */
	{
		int size = 6;
		float score[] = {5.0, 15.0, 1.0, 15.0, 5.0, 3.0};
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		int i, tmp = 0;
		for (i = 0; i < 3; ++i)
		{
			sprintf(qryName + tmp, "qry0\0");
			tmp += qryNameLen;
		}
		for (i = 3; i < 6; ++i)
		{
			sprintf(qryName + tmp, "qry1\0");
			tmp += qryNameLen;
		}

		chooseBestAlignment_wrap(qryName, qryNameLen, score, size);

		if (score[0] == -1.0 && score[1] == 15.0 && score[2] == -1.0
				&& score[3] == 15.0 && score[4] == -1.0 && score[5] == -1.0)
		{ }
		else
			fail("Incorrect behavior when there are 2 queries and size is 6."
					"\n");

		free(qryName);
	}

	/* There are 2 queries, size is 6, and the second query has 2 max scores. */
	{
		int size = 6;
		float score[] = {5.0, 15.0, 1.0, 15.0, 15.0, 3.0};
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		int i, tmp = 0;
		for (i = 0; i < 3; ++i)
		{
			sprintf(qryName + tmp, "qry0\0");
			tmp += qryNameLen;
		}
		for (i = 3; i < 6; ++i)
		{
			sprintf(qryName + tmp, "qry1\0");
			tmp += qryNameLen;
		}

		chooseBestAlignment_wrap(qryName, qryNameLen, score, size);

		if (score[0] == -1.0 && score[1] == 15.0 && score[2] == -1.0
				&& score[5] == -1.0
				&& ((score[3] == 15.0 && score[4] == -1.0)
						|| (score[3] == -1.0 && score[4] == 15.0)))
		{ }
		else
			fail("Incorrect behavior when there are 2 queries and size is 6 "
					"and second query has 2 max scores.\n");

		free(qryName);
	}

	/* There are 3 queries, size is 3. */
	{
		int size = 3;
		float score[] = {5.0, 15.0, 1.0};
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		int tmp = 0;
		sprintf(qryName + tmp, "qry0\0");
		tmp += qryNameLen;
		sprintf(qryName + tmp, "qry1\0");
		tmp += qryNameLen;
		sprintf(qryName + tmp, "qry2\0");

		chooseBestAlignment_wrap(qryName, qryNameLen, score, size);

		if (score[0] == 5.0 && score[1] == 15.0 && score[2] == 1.0)
		{ }
		else
			fail("Incorrect behavior when there are 3 queries and size is 3. "
					"\n");

		free(qryName);
	}
}
END_TEST


/**
 * Tests @a chooseBestAlignment_paired function.
 */
START_TEST(chooseBestAlignment_paired)
{
	/* There is only 1 pair. */
	{
		int size = 2;
		float score[] = {5.0, 15.0};
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		sprintf(qryName, "qry0_1\0");
		sprintf(qryName + qryNameLen, "qry0_2\0");

		chooseBestAlignment_paired_wrap(qryName, qryNameLen, score, size);

		if (score[0] == 5.0 && score[1] == 15.0)
		{ }
		else
			fail("Incorrect behavior when there is only 1 pair.\n");

		free(qryName);
	}

	/* There is only 1 pair with 2 hits. */
	{
		int size = 4;
		float score[] = {5.0, 15.0, 3.0, 18.0};
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		sprintf(qryName, "qry0_1\0");
		sprintf(qryName + qryNameLen, "qry0_2\0");
		sprintf(qryName + (2 * qryNameLen), "qry0_1\0");
		sprintf(qryName + (3 * qryNameLen), "qry0_2\0");

		chooseBestAlignment_paired_wrap(qryName, qryNameLen, score, size);

		if (score[0] == -1.0 && score[1] == -1.0 && score[2] == 3.0
				&& score[3] == 18.0)
		{ }
		else
			fail("Incorrect behavior when there is only 1 pair with 2 hits.\n");

		free(qryName);
	}

	/* There are 2 pairs with 1 hit each. */
	{
		int size = 4;
		float score[] = {5.0, 15.0, 3.0, 18.0};
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		sprintf(qryName, "qry0_1\0");
		sprintf(qryName + qryNameLen, "qry0_2\0");
		sprintf(qryName + (2 * qryNameLen), "qry1_1\0");
		sprintf(qryName + (3 * qryNameLen), "qry1_2\0");

		chooseBestAlignment_paired_wrap(qryName, qryNameLen, score, size);

		if (score[0] == 5.0 && score[1] == 15.0 && score[2] == 3.0
				&& score[3] == 18.0)
		{ }
		else
			fail("Incorrect behavior when there are 2 pairs with 1 hit each.\n");

		free(qryName);
	}

	/* There are 2 pairs with 2 hits each. */
	{
		int size = 8;
		float score[] = {5.0, 15.0, 3.0, 18.0, 18.0, 3.0, 15.0, 5.0};
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		sprintf(qryName, "qry0_1\0");
		sprintf(qryName + qryNameLen, "qry0_2\0");
		sprintf(qryName + (2 * qryNameLen), "qry0_1\0");
		sprintf(qryName + (3 * qryNameLen), "qry0_2\0");
		sprintf(qryName + (4 * qryNameLen), "qry1_1\0");
		sprintf(qryName + (5 * qryNameLen), "qry1_2\0");
		sprintf(qryName + (6 * qryNameLen), "qry1_1\0");
		sprintf(qryName + (7 * qryNameLen), "qry1_2\0");

		chooseBestAlignment_paired_wrap(qryName, qryNameLen, score, size);

		if (score[0] == -1.0 && score[1] == -1.0 && score[2] == 3.0
				&& score[3] == 18.0 && score[4] == 18.0 && score[5] == 3.0
				&& score[6] == -1.0 && score[7] == -1.0)
		{ }
		else
			fail("Incorrect behavior when there are 2 pairs with 2 hits each."
					"\n");

		free(qryName);
	}

	/* There is 1 pair with 2 hits with the same combined score. */
	{
		int size = 4;
		float score[] = {5.0, 15.0, 18.0, 2.0};
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		sprintf(qryName, "qry0_1\0");
		sprintf(qryName + qryNameLen, "qry0_2\0");
		sprintf(qryName + (2 * qryNameLen), "qry0_1\0");
		sprintf(qryName + (3 * qryNameLen), "qry0_2\0");

		chooseBestAlignment_paired_wrap(qryName, qryNameLen, score, size);

		if ((score[0] == -1.0 && score[1] == -1.0 && score[2] == 18.0
				&& score[3] == 2.0) || (score[0] == 5.0 && score[1] == 15.0
						&& score[2] == -1.0 && score[3] == -1.0))
		{ }
		else
			fail("Incorrect behavior when there is 1 pair with 2 hits with "
					"the same combined score.\n");

		free(qryName);
	}
}
END_TEST


/**
 * Tests @a chooseBestAlignment_paired2 function.
 */
START_TEST(chooseBestAlignment_paired2)
{
	/* There is only 1 pair. */
	{
		int size = 2;
		float score[] = {5.0, 15.0};
		int alignStart[] = {1, 10};
		int alignEnd[] = {5, 15};
		uint minFragSize = 0;
		uint maxFragSize = 15;
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		sprintf(qryName, "qry0_1\0");
		sprintf(qryName + qryNameLen, "qry0_2\0");

		chooseBestAlignment_paired2_wrap(qryName, qryNameLen, score,
				alignStart, alignEnd, size, minFragSize, maxFragSize);

		if (score[0] == 5.0 && score[1] == 15.0)
		{ }
		else
			fail("Incorrect behavior when there is only 1 pair.\n");

		free(qryName);
	}

	/* There is only 1 pair with 2 hits. */
	{
		int size = 4;
		float score[] = {5.0, 15.0, 3.0, 18.0};
		int alignStart[] = {0, 12, 24, 28};
		int alignEnd[] = {11, 23, 35, 38};
		uint minFragSize = 0;
		uint maxFragSize = 20;
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		sprintf(qryName, "qry0_1\0");
		sprintf(qryName + qryNameLen, "qry0_2\0");
		sprintf(qryName + (2 * qryNameLen), "qry0_1\0");
		sprintf(qryName + (3 * qryNameLen), "qry0_2\0");

		chooseBestAlignment_paired2_wrap(qryName, qryNameLen, score, alignStart,
				alignEnd, size, minFragSize, maxFragSize);

		if (score[0] == -1.0 && score[1] == -1.0 && score[2] == 3.0
				&& score[3] == 18.0)
		{ }
		else
			fail("Incorrect behavior when there is only 1 pair with 2 hits.\n");

		free(qryName);
	}

	/* There are 2 pairs with 1 hit each. */
	{
		int size = 4;
		float score[] = {5.0, 15.0, 3.0, 18.0};
		int alignStart[] = {0, 5, 24, 28};
		int alignEnd[] = {11, 13, 35, 47};
		uint minFragSize = 0;
		uint maxFragSize = 15;
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		sprintf(qryName, "qry0_1\0");
		sprintf(qryName + qryNameLen, "qry0_2\0");
		sprintf(qryName + (2 * qryNameLen), "qry1_1\0");
		sprintf(qryName + (3 * qryNameLen), "qry1_2\0");

		chooseBestAlignment_paired2_wrap(qryName, qryNameLen, score,
				alignStart, alignEnd, size, minFragSize, maxFragSize);

		if (score[0] == 5.0 && score[1] == 15.0 && score[2] == -1.0
				&& score[3] == -1.0)
		{ }
		else
			fail("Incorrect behavior when there are 2 pairs with 1 hit each.\n");

		free(qryName);
	}

	/* There are 2 pairs with 2 hits each. */
	{
		int size = 8;
		float score[] = {5.0, 15.0, 3.0, 18.0, 18.0, 3.0, 15.0, 5.0};
		int alignStart[] = {0, 10, 20, 30, 40, 50, 60, 70};
		int alignEnd[] = {9, 19, 29, 42, 49, 59, 69, 79};
		uint minFragSize = 0;
		uint maxFragSize = 20;
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		sprintf(qryName, "qry0_1\0");
		sprintf(qryName + qryNameLen, "qry0_2\0");
		sprintf(qryName + (2 * qryNameLen), "qry0_1\0");
		sprintf(qryName + (3 * qryNameLen), "qry0_2\0");
		sprintf(qryName + (4 * qryNameLen), "qry1_1\0");
		sprintf(qryName + (5 * qryNameLen), "qry1_2\0");
		sprintf(qryName + (6 * qryNameLen), "qry1_1\0");
		sprintf(qryName + (7 * qryNameLen), "qry1_2\0");

		chooseBestAlignment_paired2_wrap(qryName, qryNameLen, score, alignStart,
				alignEnd, size, minFragSize, maxFragSize);

		if (score[0] == 5.0 && score[1] == 15.0 && score[2] == -1.0
				&& score[3] == -1.0 && score[4] == 18.0 && score[5] == 3.0
				&& score[6] == -1.0 && score[7] == -1.0)
		{ }
		else
			fail("Incorrect behavior when there are 2 pairs with 2 hits each."
					"\n");

		free(qryName);
	}

	/* There is 1 pair with 2 hits with the same combined score. */
	{
		int size = 4;
		float score[] = {5.0, 15.0, 18.0, 2.0};
		int alignStart[] = {0, 10, 20, 30};
		int alignEnd[] = {9, 19, 29, 35};
		uint minFragSize = 0;
		uint maxFragSize = 15;
		int qryNameLen = 10;
		char *qryName = (char *) malloc(size * qryNameLen * sizeof(char));
		sprintf(qryName, "qry0_1\0");
		sprintf(qryName + qryNameLen, "qry0_2\0");
		sprintf(qryName + (2 * qryNameLen), "qry0_1\0");
		sprintf(qryName + (3 * qryNameLen), "qry0_2\0");

		chooseBestAlignment_paired2_wrap(qryName, qryNameLen, score, alignStart,
				alignEnd, size, minFragSize, maxFragSize);

		if (score[0] == -1.0 && score[1] == -1.0 && score[2] == 18.0
				&& score[3] == 2.0)
		{ }
		else
			fail("Incorrect behavior when there is 1 pair with 2 hits with "
					"the same combined score.\n");

		free(qryName);
	}
}
END_TEST


/**
 * Tests @a qryNameCmp function.
 */
START_TEST(qryNameCmp)
{
	/* Both query names are same. */
	{
		char *s = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
		char *t = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));

		strcpy(s, "qry0");
		strcpy(t, "qry0");
		int result = qryNameCmp_wrap(s, t, MAX_QRY_NAME_LENGTH);
		if (result != 0)
			fail("Incorrect behavior when both query names are same.\n");

		free(s);
		free(t);
	}

	/* Both query names are same and the first query name has an asterisk
	 * at the end. */
	{
		char *s = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
		char *t = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));

		strcpy(s, "qry0*");
		strcpy(t, "qry0");
		int result = qryNameCmp_wrap(s, t, MAX_QRY_NAME_LENGTH);
		if (result != 0)
			fail("Incorrect behavior when both query names are same "
					"and the first query name has an asterisk at the end.\n");

		free(s);
		free(t);
	}

	/* Both query names are same and the second query name has an asterisk
	 * at the end. */
	{
		char *s = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
		char *t = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));

		strcpy(s, "qry0");
		strcpy(t, "qry0*");
		int result = qryNameCmp_wrap(s, t, MAX_QRY_NAME_LENGTH);
		if (result != 0)
			fail("Incorrect behavior when both query names are same "
					"and the second query name has an asterisk at the end.\n");

		free(s);
		free(t);
	}

	/* Query names are not same and their lengths are same. */
	{
		char *s = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
		char *t = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));

		strcpy(s, "qry0");
		strcpy(t, "qry1");
		int result = qryNameCmp_wrap(s, t, MAX_QRY_NAME_LENGTH);
		if (result != -1)
			fail("Incorrect behavior when query names are not same"
					" and their lengths are same.\n");

		free(s);
		free(t);
	}

	/* Query names are not same and their lengths are not same (case 1). */
	{
		char *s = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
		char *t = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));

		strcpy(s, "qry0");
		strcpy(t, "qry11");
		int result = qryNameCmp_wrap(s, t, MAX_QRY_NAME_LENGTH);
		if (result != -1)
			fail("Incorrect behavior when query names are not same"
					" and their lengths are not same (case 1).\n");

		free(s);
		free(t);
	}

	/* Query names are not same and their lengths are not same (case 2). */
	{
		char *s = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
		char *t = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));

		strcpy(s, "qry0");
		strcpy(t, "qry01");
		int result = qryNameCmp_wrap(s, t, MAX_QRY_NAME_LENGTH);
		if (result != -1)
			fail("Incorrect behavior when query names are not same"
					" and their lengths are not same (case 2).\n");

		free(s);
		free(t);
	}

	/* Query names are not same and their lengths are not same (case 3). */
	{
		char *s = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
		char *t = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));

		strcpy(s, "qry01");
		strcpy(t, "qry0");
		int result = qryNameCmp_wrap(s, t, MAX_QRY_NAME_LENGTH);
		if (result != -1)
			fail("Incorrect behavior when query names are not same"
					" and their lengths are not same (case 3).\n");

		free(s);
		free(t);
	}

	/* Query names are not same and there is an asterisk
	 * somewhere in the middle of the first query. */
	{
		char *s = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
		char *t = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));

		strcpy(s, "qry0");
		strcpy(t, "qry*0");
		int result = qryNameCmp_wrap(s, t, MAX_QRY_NAME_LENGTH);
		if (result != -1)
			fail("Incorrect behavior when query names are not same"
					" and there is an asterisk somewhere in the middle "
					"of the first query.\n");

		free(s);
		free(t);
	}
}
END_TEST


/**
 * Tests @a getMatePairName function.
 */
START_TEST(getMatePairName)
{
	/* When the mate pair is the second query and the first query is not a
	 * reverse complement. */
	{
		char qryName[MAX_QRY_NAME_LENGTH] = "chr1/1";
		char mateName[MAX_QRY_NAME_LENGTH];

		getMatePairName(qryName, mateName);

		if (strcmp(mateName, "chr1/2") != 0)
			fail("Incorrect behavior when the mate pair is the second query "
					"and the first query is not a reverse complement.\n");
	}

	/* When the mate pair is the first query and the second query is not a
	 * reverse complement. */
	{
		char qryName[MAX_QRY_NAME_LENGTH] = "chr1/2";
		char mateName[MAX_QRY_NAME_LENGTH];

		getMatePairName(qryName, mateName);

		if (strcmp(mateName, "chr1/1") != 0)
			fail("Incorrect behavior when the mate pair is the first query "
					"and the second query is not a reverse complement.\n");
	}

	/* When the mate pair is the second query and the first query is a
	 * reverse complement. */
	{
		char qryName[MAX_QRY_NAME_LENGTH] = "chr1/1*";
		char mateName[MAX_QRY_NAME_LENGTH];

		getMatePairName(qryName, mateName);

		if (strcmp(mateName, "chr1/2") != 0)
			fail("Incorrect behavior when the mate pair is the second query "
					"and the first query is a reverse complement.\n");
	}

	/* When the mate pair is the first query and the second query is a
	 * reverse complement. */
	{
		char qryName[MAX_QRY_NAME_LENGTH] = "chr1/2*";
		char mateName[MAX_QRY_NAME_LENGTH];

		getMatePairName(qryName, mateName);

		if (strcmp(mateName, "chr1/1") != 0)
			fail("Incorrect behavior when the mate pair is the first query "
					"and the second query is a reverse complement.\n");
	}
}
END_TEST


/**
 * Tests @a chooseBestAlignment_paired3 function.
 */
START_TEST(chooseBestAlignment_paired3)
{
	/* Two queries, no mate. */
	{
		int numQrs = 2;
		int qryNameLen = MAX_QRY_NAME_LENGTH;
		char *qryName = (char *) malloc(numQrs * qryNameLen * sizeof(char));
		strcpy(qryName, "qry1/1");
		strcpy(qryName + qryNameLen, "qry1/1*");
		int refNameLen = MAX_REF_NAME_LENGTH;
		char *refName = (char *) malloc(numQrs * refNameLen * sizeof(char));
		strcpy(refName, "ref1");
		strcpy(refName + MAX_REF_NAME_LENGTH, "ref1");
		float score[] = {5.0, 12.0};
		int refDistance[] = {3, 80};
		int alignStart[] = {2, 5};
		int alignEnd[] = {15, 20};
		int minFragSize = 5;
		int maxFragSize = 35;
		uint revComp[] = {0, 1};

		chooseBestAlignment_paired3_wrap(qryName, qryNameLen, revComp, refName,
				refNameLen, score, refDistance, alignStart, alignEnd, numQrs,
				minFragSize, maxFragSize);

		if (score[0] == -1.0 && score[1] == -1.0)
		{ }
		else
			fail("Incorrect behavior when there are 2 queries, no mate.\n");
	}

	/* Three queries of the same pair. */
	{
		int numQrs = 3;
		int qryNameLen = MAX_QRY_NAME_LENGTH;
		char *qryName = (char *) malloc(numQrs * qryNameLen * sizeof(char));
		strcpy(qryName, "qry1/1");
		strcpy(qryName + qryNameLen, "qry1/1*");
		strcpy(qryName + (2 * qryNameLen), "qry1/2*");
		int refNameLen = MAX_REF_NAME_LENGTH;
		char *refName = (char *) malloc(numQrs * refNameLen * sizeof(char));
		strcpy(refName, "ref1");
		strcpy(refName + MAX_REF_NAME_LENGTH, "ref1");
		strcpy(refName + (2 * MAX_REF_NAME_LENGTH), "ref1");
		float score[] = {5.0, 12.0, 9.0};
		int refDistance[] = {3, 80, 25};
		int alignStart[] = {2, 5, 0};
		int alignEnd[] = {15, 20, 12};
		int minFragSize = 5;
		int maxFragSize = 35;
		uint revComp[] = {0, 1, 1};

		chooseBestAlignment_paired3_wrap(qryName, qryNameLen, revComp, refName,
				refNameLen, score, refDistance, alignStart, alignEnd, numQrs,
				minFragSize, maxFragSize);

		if (score[0] == 5.0 && score[1] == -1.0 && score[2] == 9.0)
		{ }
		else
			fail("Incorrect behavior when there are 3 queries of the same "
					"pair.\n");
	}

	/* Four queries, 3 belonging to the same pair. */
	{
		int numQrs = 4;
		int qryNameLen = MAX_QRY_NAME_LENGTH;
		char *qryName = (char *) malloc(numQrs * qryNameLen * sizeof(char));
		strcpy(qryName, "qry1/1");
		strcpy(qryName + qryNameLen, "qry1/1*");
		strcpy(qryName + (2 * qryNameLen), "qry1/2*");
		strcpy(qryName + (3 * qryNameLen), "qry2/1");
		int refNameLen = MAX_REF_NAME_LENGTH;
		char *refName = (char *) malloc(numQrs * refNameLen * sizeof(char));
		strcpy(refName, "ref1");
		strcpy(refName + MAX_REF_NAME_LENGTH, "ref1");
		strcpy(refName + (2 * MAX_REF_NAME_LENGTH), "ref1");
		strcpy(refName + (2 * MAX_REF_NAME_LENGTH), "ref1");
		float score[] = {5.0, 12.0, 9.0, 10.0};
		int refDistance[] = {3, 80, 25, 13};
		int alignStart[] = {2, 5, 0, 8};
		int alignEnd[] = {15, 20, 12, 27};
		int minFragSize = 5;
		int maxFragSize = 35;
		uint revComp[] = {0, 1, 1, 0};

		chooseBestAlignment_paired3_wrap(qryName, qryNameLen, revComp, refName,
				refNameLen, score, refDistance, alignStart, alignEnd, numQrs,
				minFragSize, maxFragSize);

		if (score[0] == 5.0 && score[1] == -1.0 && score[2] == 9.0
				&& score[3] == -1.0)
		{ }
		else
			fail("Incorrect behavior when there are 4 queries, 3 belonging to "
					" the same pair.\n");
	}

	/* Four queries, all belonging to the same pair. */
	{
		int numQrs = 4;
		int qryNameLen = MAX_QRY_NAME_LENGTH;
		char *qryName = (char *) malloc(numQrs * qryNameLen * sizeof(char));
		strcpy(qryName, "qry1/1");
		strcpy(qryName + qryNameLen, "qry1/1*");
		strcpy(qryName + (2 * qryNameLen), "qry1/2");
		strcpy(qryName + (3 * qryNameLen), "qry1/2*");
		int refNameLen = MAX_REF_NAME_LENGTH;
		char *refName = (char *) malloc(numQrs * refNameLen * sizeof(char));
		strcpy(refName, "ref1");
		strcpy(refName + MAX_REF_NAME_LENGTH, "ref1");
		strcpy(refName + (2 * MAX_REF_NAME_LENGTH), "ref1");
		strcpy(refName + (3 * MAX_REF_NAME_LENGTH), "ref1");
		float score[] = {5.0, 12.0, 9.0, 10.0};
		int refDistance[] = {3, 80, 25, 13};
		int alignStart[] = {2, 5, 0, 8};
		int alignEnd[] = {15, 20, 12, 27};
		int minFragSize = 5;
		int maxFragSize = 35;
		uint revComp[] = {0, 1, 0, 1};

		chooseBestAlignment_paired3_wrap(qryName, qryNameLen, revComp, refName,
				refNameLen, score, refDistance, alignStart, alignEnd, numQrs,
				minFragSize, maxFragSize);

		if (score[0] == 5.0 && score[1] == -1.0 && score[2] == -1.0
				&& score[3] == 10.0)
		{ }
		else
			fail("Incorrect behavior when there are 4 queries, all belonging to"
					" the same pair.\n");
	}

	/* Five queries, last one has no mate. */
	{
		int numQrs = 5;
		int qryNameLen = MAX_QRY_NAME_LENGTH;
		char *qryName = (char *) malloc(numQrs * qryNameLen * sizeof(char));
		strcpy(qryName, "qry1/1");
		strcpy(qryName + qryNameLen, "qry1/1*");
		strcpy(qryName + (2 * qryNameLen), "qry1/2");
		strcpy(qryName + (3 * qryNameLen), "qry1/2*");
		strcpy(qryName + (4 * qryNameLen), "qry2/2");
		int refNameLen = MAX_REF_NAME_LENGTH;
		char *refName = (char *) malloc(numQrs * refNameLen * sizeof(char));
		strcpy(refName, "ref1");
		strcpy(refName + MAX_REF_NAME_LENGTH, "ref1");
		strcpy(refName + (2 * MAX_REF_NAME_LENGTH), "ref1");
		strcpy(refName + (3 * MAX_REF_NAME_LENGTH), "ref1");
		strcpy(refName + (4 * MAX_REF_NAME_LENGTH), "ref1");
		float score[] = {5.0, 12.0, 9.0, 10.0, 14.0};
		int refDistance[] = {3, 80, 25, 13, 7};
		int alignStart[] = {2, 5, 0, 8, 4};
		int alignEnd[] = {15, 20, 12, 27, 16};
		int minFragSize = 5;
		int maxFragSize = 35;
		uint revComp[] = {0, 1, 0, 1, 0};

		chooseBestAlignment_paired3_wrap(qryName, qryNameLen, revComp, refName,
				refNameLen, score, refDistance, alignStart, alignEnd, numQrs,
				minFragSize, maxFragSize);

		if (score[0] == 5.0 && score[1] == -1.0 && score[2] == -1.0
				&& score[3] == 10.0 && score[4] == -1.0)
		{ }
		else
			fail("Incorrect behavior when there are 5 queries, last one has "
					"no mate.\n");
	}

	/* Three queries, all belonging to the same pair, and have same scores. */
	{
		int numQrs = 3;
		int qryNameLen = MAX_QRY_NAME_LENGTH;
		char *qryName = (char *) malloc(numQrs * qryNameLen * sizeof(char));
		strcpy(qryName, "qry1/1");
		strcpy(qryName + qryNameLen, "qry1/1");
		strcpy(qryName + (2 * qryNameLen), "qry1/2*");
		int refNameLen = MAX_REF_NAME_LENGTH;
		char *refName = (char *) malloc(numQrs * refNameLen * sizeof(char));
		strcpy(refName, "ref1");
		strcpy(refName + MAX_REF_NAME_LENGTH, "ref1");
		strcpy(refName + (2 * MAX_REF_NAME_LENGTH), "ref1");
		float score[] = {5.0, 5.0, 9.0};
		int refDistance[] = {3, 8, 25};
		int alignStart[] = {2, 5, 0};
		int alignEnd[] = {15, 20, 12};
		int minFragSize = 5;
		int maxFragSize = 35;
		uint revComp[] = {0, 0, 1};

		chooseBestAlignment_paired3_wrap(qryName, qryNameLen, revComp, refName,
				refNameLen, score, refDistance, alignStart, alignEnd, numQrs,
				minFragSize, maxFragSize);

		if (((score[0] == 5.0 && score[1] == -1.0)
				|| (score[0] == -1.0 && score[1] == 5.0)) && score[2] == 9.0)
		{ }
		else
			fail("Incorrect behavior when there are 3 queries of the same "
					"pair and have the same scores.\n");
	}

	/* Three queries, all belonging to the same pair, and map
	 * to different references. */
	{
		int numQrs = 3;
		int qryNameLen = MAX_QRY_NAME_LENGTH;
		char *qryName = (char *) malloc(numQrs * qryNameLen * sizeof(char));
		strcpy(qryName, "qry1/1");
		strcpy(qryName + qryNameLen, "qry1/1");
		strcpy(qryName + (2 * qryNameLen), "qry1/2*");
		int refNameLen = MAX_REF_NAME_LENGTH;
		char *refName = (char *) malloc(numQrs * refNameLen * sizeof(char));
		strcpy(refName, "ref1");
		strcpy(refName + MAX_REF_NAME_LENGTH, "ref2");
		strcpy(refName + (2 * MAX_REF_NAME_LENGTH), "ref1");
		float score[] = {5.0, 8.0, 9.0};
		int refDistance[] = {3, 8, 25};
		int alignStart[] = {2, 5, 0};
		int alignEnd[] = {15, 20, 12};
		int minFragSize = 5;
		int maxFragSize = 35;
		uint revComp[] = {0, 0, 1};

		chooseBestAlignment_paired3_wrap(qryName, qryNameLen, revComp, refName,
				refNameLen, score, refDistance, alignStart, alignEnd, numQrs,
				minFragSize, maxFragSize);

		if (score[0] == -1.0 && score[1] == 8.0 && score[2] == 9.0)
		{ }
		else
			fail("Incorrect behavior when there are 3 queries of the same "
					"pair and map to different references.\n");
	}

	/* Three queries, first one has no mate. */
	{
		int numQrs = 3;
		int qryNameLen = MAX_QRY_NAME_LENGTH;
		char *qryName = (char *) malloc(numQrs * qryNameLen * sizeof(char));
		strcpy(qryName, "qry0/1");
		strcpy(qryName + qryNameLen, "qry1/1");
		strcpy(qryName + (2 * qryNameLen), "qry1/2*");
		int refNameLen = MAX_REF_NAME_LENGTH;
		char *refName = (char *) malloc(numQrs * refNameLen * sizeof(char));
		strcpy(refName, "ref1");
		strcpy(refName + MAX_REF_NAME_LENGTH, "ref2");
		strcpy(refName + (2 * MAX_REF_NAME_LENGTH), "ref1");
		float score[] = {5.0, 8.0, 9.0};
		int refDistance[] = {3, 8, 25};
		int alignStart[] = {2, 5, 0};
		int alignEnd[] = {15, 20, 12};
		int minFragSize = 5;
		int maxFragSize = 35;
		uint revComp[] = {0, 0, 1};

		chooseBestAlignment_paired3_wrap(qryName, qryNameLen, revComp, refName,
				refNameLen, score, refDistance, alignStart, alignEnd, numQrs,
				minFragSize, maxFragSize);

		if (score[0] == -1.0 && score[1] == 8.0 && score[2] == 9.0)
		{ }
		else
			fail("Incorrect behavior when there are 3 queries, first one "
					"has no mate.\n");
	}
}
END_TEST


/**
 * Tests @a getNumQrsPerIter_paired function.
 */
START_TEST(getNumQrsPerIter_paired)
{
	/* Max queries per iteration is 5 and there are 4 hits belonging to the
	 * same query. */
	{
		int maxNumQrsPerIter = 5;
		char mapFile[MAX_FILE_NAME_LENGTH];
		sprintf(mapFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_map.txt");
		FILE *mapFilePtr = fopen(mapFile, "w");
		char map1[] = "1\tq1/1\n";
		char map2[] = "1\tq1/1\n";
		char map3[] = "1\tq1/1*\n";
		char map4[] = "1\tq1/1*\n";
		fputs(map1, mapFilePtr);
		fputs(map2, mapFilePtr);
		fputs(map3, mapFilePtr);
		fputs(map4, mapFilePtr);
		fclose(mapFilePtr);

		int *numQrsPerIter = NULL;
		int arrSize = 0;

		getNumQrsPerIter_paired_wrap(mapFile, maxNumQrsPerIter, &numQrsPerIter,
				&arrSize);
		if (arrSize != 1 || numQrsPerIter[0] != 4)
			fail("Incorrect behavior when the max queries per iteration is 5 "
					"and there are 4 hits belonging to the same query.\n");

		remove(mapFile);
	}

	/* Max queries per iteration is 5 and there are 5 hits belonging to the
	 * same query. */
	{
		int maxNumQrsPerIter = 5;
		char mapFile[MAX_FILE_NAME_LENGTH];
		sprintf(mapFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_map.txt");
		FILE *mapFilePtr = fopen(mapFile, "w");
		char map1[] = "1\tq1/1\n";
		char map2[] = "1\tq1/1\n";
		char map3[] = "1\tq1/1*\n";
		char map4[] = "1\tq1/1*\n";
		char map5[] = "1\tq1/2\n";
		fputs(map1, mapFilePtr);
		fputs(map2, mapFilePtr);
		fputs(map3, mapFilePtr);
		fputs(map4, mapFilePtr);
		fputs(map5, mapFilePtr);
		fclose(mapFilePtr);

		int *numQrsPerIter = NULL;
		int arrSize = 0;

		getNumQrsPerIter_paired_wrap(mapFile, maxNumQrsPerIter, &numQrsPerIter,
				&arrSize);
		if (arrSize != 1 || numQrsPerIter[0] != 5)
			fail("Incorrect behavior when the max queries per iteration is 5 "
					"and there are 5 hits belonging to the same query.\n");

		remove(mapFile);
	}

	/* Max queries per iteration is 5 and there are 5 hits belonging to 2
	 * queries. */
	{
		int maxNumQrsPerIter = 5;
		char mapFile[MAX_FILE_NAME_LENGTH];
		sprintf(mapFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_map.txt");
		FILE *mapFilePtr = fopen(mapFile, "w");
		char map1[] = "1\tq1/1\n";
		char map2[] = "1\tq1/1\n";
		char map3[] = "1\tq1/1*\n";
		char map4[] = "2\tq2/1*\n";
		char map5[] = "2\tq2/2\n";
		fputs(map1, mapFilePtr);
		fputs(map2, mapFilePtr);
		fputs(map3, mapFilePtr);
		fputs(map4, mapFilePtr);
		fputs(map5, mapFilePtr);
		fclose(mapFilePtr);

		int *numQrsPerIter = NULL;
		int arrSize = 0;

		getNumQrsPerIter_paired_wrap(mapFile, maxNumQrsPerIter, &numQrsPerIter,
				&arrSize);
		if (arrSize != 1 || numQrsPerIter[0] != 5)
			fail("Incorrect behavior when the max queries per iteration is 5 "
					"and there are 5 hits belonging to 2 queries.\n");

		remove(mapFile);
	}

	/* Max queries per iteration is 5 and there are 6 hits belonging to
	 * 2 queries. */
	{
		int maxNumQrsPerIter = 5;
		char mapFile[MAX_FILE_NAME_LENGTH];
		sprintf(mapFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_map.txt");
		FILE *mapFilePtr = fopen(mapFile, "w");
		char map1[] = "1\n";
		char map2[] = "1\n";
		char map3[] = "1\n";
		char map4[] = "2\n";
		char map5[] = "2\n";
		char map6[] = "2\n";
		fputs(map1, mapFilePtr);
		fputs(map2, mapFilePtr);
		fputs(map3, mapFilePtr);
		fputs(map4, mapFilePtr);
		fputs(map5, mapFilePtr);
		fputs(map6, mapFilePtr);
		fclose(mapFilePtr);

		int *numQrsPerIter = NULL;
		int arrSize = 0;

		getNumQrsPerIter_paired_wrap(mapFile, maxNumQrsPerIter, &numQrsPerIter,
				&arrSize);
		if (arrSize != 2 || numQrsPerIter[0] != 3 || numQrsPerIter[1] != 3)
			fail("Incorrect behavior when the max queries per iteration is 5 "
					"and there are 6 hits belonging to 2 queries.\n");

		remove(mapFile);
	}

	/* Max queries per iteration is 5 and there are 8 hits belonging to
	 * 3 queries. */
	{
		int maxNumQrsPerIter = 5;
		char mapFile[MAX_FILE_NAME_LENGTH];
		sprintf(mapFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_map.txt");
		FILE *mapFilePtr = fopen(mapFile, "w");
		char map1[] = "1\n";
		char map2[] = "1\n";
		char map3[] = "1\n";
		char map4[] = "2\n";
		char map5[] = "2\n";
		char map6[] = "2\n";
		char map7[] = "3\n";
		char map8[] = "3\n";
		fputs(map1, mapFilePtr);
		fputs(map2, mapFilePtr);
		fputs(map3, mapFilePtr);
		fputs(map4, mapFilePtr);
		fputs(map5, mapFilePtr);
		fputs(map6, mapFilePtr);
		fputs(map7, mapFilePtr);
		fputs(map8, mapFilePtr);
		fclose(mapFilePtr);

		int *numQrsPerIter = NULL;
		int arrSize = 0;

		getNumQrsPerIter_paired_wrap(mapFile, maxNumQrsPerIter, &numQrsPerIter,
				&arrSize);
		if (arrSize != 2 || numQrsPerIter[0] != 3 || numQrsPerIter[1] != 5)
			fail("Incorrect behavior when the max queries per iteration is 5 "
					"and there are 8 hits belonging to 3 queries.\n");

		remove(mapFile);
	}

	/* Max queries per iteration is 5 and there are 9 hits belonging to
	 * 3 queries. */
	{
		int maxNumQrsPerIter = 5;
		char mapFile[MAX_FILE_NAME_LENGTH];
		sprintf(mapFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_map.txt");
		FILE *mapFilePtr = fopen(mapFile, "w");
		char map1[] = "1\n";
		char map2[] = "1\n";
		char map3[] = "1\n";
		char map4[] = "2\n";
		char map5[] = "2\n";
		char map6[] = "2\n";
		char map7[] = "3\n";
		char map8[] = "3\n";
		char map9[] = "3\n";
		fputs(map1, mapFilePtr);
		fputs(map2, mapFilePtr);
		fputs(map3, mapFilePtr);
		fputs(map4, mapFilePtr);
		fputs(map5, mapFilePtr);
		fputs(map6, mapFilePtr);
		fputs(map7, mapFilePtr);
		fputs(map8, mapFilePtr);
		fputs(map9, mapFilePtr);
		fclose(mapFilePtr);

		int *numQrsPerIter = NULL;
		int arrSize = 0;

		getNumQrsPerIter_paired_wrap(mapFile, maxNumQrsPerIter, &numQrsPerIter,
				&arrSize);
		if (arrSize != 3 || numQrsPerIter[0] != 3 || numQrsPerIter[1] != 3
				|| numQrsPerIter[2] != 3)
			fail("Incorrect behavior when the max queries per iteration is 5 "
					"and there are 9 hits belonging to 3 queries.\n");

		remove(mapFile);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *alignSuite(void)
{
	Suite *s = suite_create("align");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, alignSequences);
	tcase_add_test(testCaseCore, chooseBestAlignment);
	tcase_add_test(testCaseCore, chooseBestAlignment_paired);
	tcase_add_test(testCaseCore, chooseBestAlignment_paired2);
	tcase_add_test(testCaseCore, chooseBestAlignment_paired3);
	tcase_add_test(testCaseCore, qryNameCmp);
	tcase_add_test(testCaseCore, getMatePairName);
	tcase_add_test(testCaseCore, getNumQrsPerIter_paired);
	suite_add_tcase (s, testCaseCore);

	return s;
}
