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
#include "testLookupTable7.h"
#include "../src/lookupTable7.h"
#include "../src/refMap.h"
#include "../src/common.h"
#include "../src/preprocess.h"
#include <limits.h>


/**
 * Tests @a createClusters_gpu function.
 */
START_TEST(createClusters_gpu)
{
	/* Case 1. */
	{
		int shift[] = {9, 8, 3, 100, 3};
		int refPos[] = {37, 13, 22, 65, 10};
		short start = 0;
		short end = 4;
		int size = 5;

		int *shift_d;
		cudaMalloc(&shift_d, size * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shift_d, shift, size * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *refPos_d;
		cudaMalloc(&refPos_d, size * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos_d, refPos, size * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		short *clusterSize_d;
		cudaMalloc(&clusterSize_d, size * sizeof(short));
		PRINT_CUDA_ERROR()

		int *distinctShift_d;
		cudaMalloc(&distinctShift_d, size * sizeof(int));
		PRINT_CUDA_ERROR()

		createClusters_gpu7_2_wrap<<<1, 1>>>(shift_d, refPos_d, start, end,
				clusterSize_d, distinctShift_d, size);
		PRINT_CUDA_ERROR()

		short clusterSize[size];
		cudaMemcpy(clusterSize, clusterSize_d, size * sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int distinctShift[size];
		cudaMemcpy(distinctShift, distinctShift_d, size * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shift, shift_d, size * sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (clusterSize[0] != 1 || clusterSize[1] != 1 || clusterSize[2] != 2
				|| clusterSize[3] != 1 || clusterSize[4] != 0
				|| distinctShift[0] != 9 || distinctShift[1] != 8
				|| distinctShift[2] != 3 || distinctShift[3] != 100
				|| distinctShift[4] != -1 || shift[0] != 0 || shift[1] != 1
				|| shift[2] != 4 || shift[3] != 3 || shift[4] != 3)
			fail("Incorrect behavior (case 1).\n");

		cudaFree(shift_d);
		cudaFree(refPos_d);
		cudaFree(clusterSize_d);
		cudaFree(distinctShift_d);
	}

	/* Case 2. */
	{
		int shift[] = {9, 8, 3, 100, 3};
		int refPos[] = {37, 13, 22, 65, 10};
		short start = 1;
		short end = 4;
		int size = 5;

		int *shift_d;
		cudaMalloc(&shift_d, size * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shift_d, shift, size * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *refPos_d;
		cudaMalloc(&refPos_d, size * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos_d, refPos, size * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		short *clusterSize_d;
		cudaMalloc(&clusterSize_d, size * sizeof(short));
		PRINT_CUDA_ERROR()

		int *distinctShift_d;
		cudaMalloc(&distinctShift_d, size * sizeof(int));
		PRINT_CUDA_ERROR()

		createClusters_gpu7_2_wrap<<<1, 1>>>(shift_d, refPos_d, start, end,
				clusterSize_d, distinctShift_d, size);
		PRINT_CUDA_ERROR()

		short clusterSize[size];
		cudaMemcpy(clusterSize, clusterSize_d, size * sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int distinctShift[size];
		cudaMemcpy(distinctShift, distinctShift_d, size * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shift, shift_d, size * sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (clusterSize[0] != 0 || clusterSize[1] != 1 || clusterSize[2] != 2
				|| clusterSize[3] != 1 || clusterSize[4] != 0
				|| distinctShift[0] != -1 || distinctShift[1] != 8
				|| distinctShift[2] != 3 || distinctShift[3] != 100
				|| distinctShift[4] != -1 || shift[0] != 9 || shift[1] != 1
				|| shift[2] != 4 || shift[3] != 3 || shift[4] != 3)
			fail("Incorrect behavior (case 2).\n");

		cudaFree(shift_d);
		cudaFree(refPos_d);
		cudaFree(clusterSize_d);
		cudaFree(distinctShift_d);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *lookupTable7Suite(void)
{
	Suite *s = suite_create("lookupTable7");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, createClusters_gpu);
	suite_add_tcase (s, testCaseCore);

	return s;
}
