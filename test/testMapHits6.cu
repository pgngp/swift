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
#include "../src/mapHits6.h"
#include "testMapHits6.h"


/**
 * Tests @a swap function.
 */
START_TEST(swap)
{
	int maxNumHits = 3;
	char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
	int *shift = (int *) calloc(maxNumHits, sizeof(int));
	int *refPos = (int *) calloc(maxNumHits, sizeof(int));
	refIdx[0] = 1;
	shift[0] = 11;
	refPos[0] = 111;
	refIdx[1] = 2;
	shift[1] = 22;
	refPos[1] = 222;
	refIdx[2] = 3;
	shift[2] = 33;
	refPos[2] = 333;
	swap_wrap6(1, 2, refIdx, shift, refPos);
	if (refIdx[1] != 3 || refIdx[2] != 2 || shift[1] != 33 || shift[2] != 22
			|| refPos[1] != 333 || refPos[2] != 222)
		fail("Incorrect behavior.\n");
	free(refIdx);
	free(shift);
	free(refPos);
}
END_TEST


/**
 * Tests @a quickSort function.
 */
START_TEST(quickSort)
{
	/* Case 1. */
	{
		int maxNumHits = 3;
		char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
		int *shift = (int *) calloc(maxNumHits, sizeof(int));
		int *refPos = (int *) calloc(maxNumHits, sizeof(int));
		refIdx[0] = 2;
		shift[0] = 22;
		refPos[0] = 222;
		refIdx[1] = 3;
		shift[1] = 33;
		refPos[1] = 333;
		refIdx[2] = 1;
		shift[2] = 11;
		refPos[2] = 111;
		quickSort_wrap6(refIdx, shift, refPos, maxNumHits);
		if (refIdx[0] != 1 || shift[0] != 11 || refPos[0] != 111
				|| refIdx[1] != 2 || shift[1] != 22 || refPos[1] != 222
				|| refIdx[2] != 3 || shift[2] != 33 || refPos[2] != 333)
			fail("Incorrect behavior.\n");
		free(refIdx);
		free(shift);
		free(refPos);
	}
}
END_TEST


/**
 * Tests @a createClusters function.
 */
START_TEST(createClusters)
{
	/* Case 1. */
	{
		int maxNumHits = 3;
		char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
		int *shift = (int *) calloc(maxNumHits, sizeof(int));
		int *refPos = (int *) calloc(maxNumHits, sizeof(int));
		refIdx[0] = 1;
		shift[0] = 5;
		refPos[0] = 20;
		refIdx[1] = 2;
		shift[1] = -1;
		refPos[1] = 10;
		refIdx[2] = 2;
		shift[2] = -1;
		refPos[2] = 30;
		int *clusters = (int *) calloc(maxNumHits, sizeof(int));
		int numClusters = 0;
		createClusters_wrap6(refIdx, shift, refPos, maxNumHits, clusters,
				&numClusters);
		if (numClusters != 2 || clusters[0] != 0 || clusters[1] != 1
				|| clusters[2] != 0)
			fail("Incorrect behavior.\n");
		free(refIdx);
		free(shift);
		free(refPos);
		free(clusters);
	}
}
END_TEST


/**
 * Tests @a findBiggestClusterSize function.
 */
START_TEST(findBiggestClusterSize)
{
	/* Case 1. */
	{
		int maxNumHits = 3;
		char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
		int *shift = (int *) calloc(maxNumHits, sizeof(int));
		int *refPos = (int *) calloc(maxNumHits, sizeof(int));
		refIdx[0] = 1;
		shift[0] = 5;
		refPos[0] = 20;
		refIdx[1] = 2;
		shift[1] = -1;
		refPos[1] = 10;
		refIdx[2] = 2;
		shift[2] = -1;
		refPos[2] = 30;
		int *clusters = (int *) calloc(maxNumHits, sizeof(int));
		int numClusters = 0;
		int *clusterSize = (int *) calloc(maxNumHits, sizeof(int));
		int biggestClusterSize = 0;
		findBiggestClusterSize_wrap6(refIdx, shift, refPos, maxNumHits,
				clusters, &numClusters, clusterSize, &biggestClusterSize);
		if (clusterSize[0] != 1 || clusterSize[1] != 2 || clusterSize[2] != 0
				|| biggestClusterSize != 2)
			fail("Incorrect behavior.\n");
		free(refIdx);
		free(shift);
		free(refPos);
		free(clusters);
		free(clusterSize);
	}
}
END_TEST


/**
 * Tests @a findBiggestClusters function.
 */
START_TEST(findBiggestClusters)
{
	/* Case 1. */
	{
		int maxNumHits = 3;
		char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
		int *shift = (int *) calloc(maxNumHits, sizeof(int));
		int *refPos = (int *) calloc(maxNumHits, sizeof(int));
		refIdx[0] = 1;
		shift[0] = 5;
		refPos[0] = 20;
		refIdx[1] = 2;
		shift[1] = -1;
		refPos[1] = 10;
		refIdx[2] = 2;
		shift[2] = -1;
		refPos[2] = 30;

		char *refIdx_bestHits = (char *) calloc(maxNumHits, sizeof(char));
		int *shift_bestHits = (int *) calloc(maxNumHits, sizeof(int));
		int *refPos_bestHits = (int *) calloc(maxNumHits, sizeof(int));

		int numHits = findBiggestClusters_wrap6(maxNumHits, refIdx_bestHits,
				shift_bestHits, refPos_bestHits, refIdx, shift, refPos,
				maxNumHits);
		if (numHits != 1 || refIdx_bestHits[0] != 2 || shift_bestHits[0] != -1
				|| refPos_bestHits[0] != 10)
			fail("Incorrect behavior.\n");

		free(refIdx);
		free(shift);
		free(refPos);
		free(refIdx_bestHits);
		free(shift_bestHits);
		free(refPos_bestHits);
	}
}
END_TEST


/**
 * Tests @a mapHits6GetBestHits function.
 */
START_TEST(mapHits6GetBestHits)
{
	/* When max number of hits desired is 3 and actual number of best hits
	 * is 1. */
	{
		int maxNumHits = 3;
		int size = 3;
		char *refIdx = (char *) calloc(size, sizeof(char));
		int *shift = (int *) calloc(size, sizeof(int));
		int *refPos = (int *) calloc(size, sizeof(int));
		refIdx[0] = 1;
		shift[0] = 5;
		refPos[0] = 20;
		refIdx[1] = 2;
		shift[1] = -1;
		refPos[1] = 10;
		refIdx[2] = 2;
		shift[2] = -1;
		refPos[2] = 30;
		char *refIdx_bestHits = (char *) calloc(size, sizeof(char));
		int *shift_bestHits = (int *) calloc(size, sizeof(int));
		int *refPos_bestHits = (int *) calloc(size, sizeof(int));
		int numHits = mapHits6GetBestHits(maxNumHits, refIdx, shift, refPos, size,
				refIdx_bestHits, shift_bestHits, refPos_bestHits);
		if (numHits != 1 || refIdx_bestHits[0] != 2 || shift_bestHits[0] != -1
				|| refPos_bestHits[0] != 10)
			fail("Incorrect behavior when max number of hits desired is 3 and "
					"actual number of best hits is 1.\n");
		free(refIdx);
		free(shift);
		free(refPos);
		free(refIdx_bestHits);
		free(shift_bestHits);
		free(refPos_bestHits);
	}

	/* When max number of hits desired is 3 and actual number of best hits
	 * is 2. */
	{
		int maxNumHits = 3;
		int size = 5;
		char *refIdx = (char *) calloc(size, sizeof(char));
		int *shift = (int *) calloc(size, sizeof(int));
		int *refPos = (int *) calloc(size, sizeof(int));
		refIdx[0] = 1;
		shift[0] = 5;
		refPos[0] = 20;
		refIdx[1] = 2;
		shift[1] = -1;
		refPos[1] = 10;
		refIdx[2] = 2;
		shift[2] = -1;
		refPos[2] = 30;
		refIdx[3] = 3;
		shift[3] = 5;
		refPos[3] = 30;
		refIdx[4] = 3;
		shift[4] = 5;
		refPos[4] = 40;
		char *refIdx_bestHits = (char *) calloc(size, sizeof(char));
		int *shift_bestHits = (int *) calloc(size, sizeof(int));
		int *refPos_bestHits = (int *) calloc(size, sizeof(int));
		int numHits = mapHits6GetBestHits(maxNumHits, refIdx, shift, refPos,
				size, refIdx_bestHits, shift_bestHits, refPos_bestHits);
		if (numHits != 2 || refIdx_bestHits[0] != 2 || shift_bestHits[0] != -1
				|| refPos_bestHits[0] != 10 || refIdx_bestHits[1] != 3
				|| shift_bestHits[1] != 5 || refPos_bestHits[1] != 30)
			fail("Incorrect behavior when max number of hits desired is 3 and "
					"actual number of best hits is 2.\n");
		free(refIdx);
		free(shift);
		free(refPos);
		free(refIdx_bestHits);
		free(shift_bestHits);
		free(refPos_bestHits);
	}

	/* When max number of hits desired is 1 and actual number of best hits
	 * is 2. */
	{
		int maxNumHits = 1;
		int size = 5;
		char *refIdx = (char *) calloc(size, sizeof(char));
		int *shift = (int *) calloc(size, sizeof(int));
		int *refPos = (int *) calloc(size, sizeof(int));
		refIdx[0] = 1;
		shift[0] = 5;
		refPos[0] = 20;
		refIdx[1] = 2;
		shift[1] = -1;
		refPos[1] = 10;
		refIdx[2] = 2;
		shift[2] = -1;
		refPos[2] = 30;
		refIdx[3] = 3;
		shift[3] = 5;
		refPos[3] = 30;
		refIdx[4] = 3;
		shift[4] = 5;
		refPos[4] = 40;
		char *refIdx_bestHits = (char *) calloc(size, sizeof(char));
		int *shift_bestHits = (int *) calloc(size, sizeof(int));
		int *refPos_bestHits = (int *) calloc(size, sizeof(int));
		int numHits = mapHits6GetBestHits(maxNumHits, refIdx, shift, refPos,
				size, refIdx_bestHits, shift_bestHits, refPos_bestHits);
		if (numHits != 1 || (refIdx_bestHits[0] != 2 && refIdx_bestHits[0] != 3)
				|| (shift_bestHits[0] != -1 && shift_bestHits[0] != 5)
				|| (refPos_bestHits[0] != 10 && refPos_bestHits[0] != 30))
			fail("Incorrect behavior when max number of hits desired is 1 and "
					"actual number of best hits is 2.\n");
		free(refIdx);
		free(shift);
		free(refPos);
		free(refIdx_bestHits);
		free(shift_bestHits);
		free(refPos_bestHits);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *mapHits6Suite(void)
{
	Suite *s = suite_create("mapHits6");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, swap);
	tcase_add_test(testCaseCore, quickSort);
	tcase_add_test(testCaseCore, createClusters);
	tcase_add_test(testCaseCore, findBiggestClusterSize);
	tcase_add_test(testCaseCore, findBiggestClusters);
	tcase_add_test(testCaseCore, mapHits6GetBestHits);
	suite_add_tcase (s, testCaseCore);

	return s;
}
