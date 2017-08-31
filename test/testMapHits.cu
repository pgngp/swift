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
#include "../src/mapHits.h"
#include "testMapHits.h"


/**
 * Tests @a mapHitsAddHit function.
 */
START_TEST(mapHitsAddHit)
{
	/* One hit is added. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		char refIdx = 1;
		int shift = 2;
		int refPos = 5;
		char *refIdxArr;
		int *shiftArr, *refPosArr, size;
		mapHitsAddHit_wrap(&refIdxArr, &shiftArr, &refPosArr, &size, refIdx,
				shift, refPos);
		if (refIdxArr[0] != 1 || shiftArr[0] != 2 || refPosArr[0] != 5
				|| size != 1)
			fail("Incorrect behavior when a hit is added.\n");
		mapHitsDelete();
	}
}
END_TEST


/**
 * Tests @a mapHitsGetBestHits function.
 */
START_TEST(mapHitsGetBestHits)
{
	/* There are no hits and number of desired hits is only 1. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHitsGetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 0)
			fail("Incorrect behavior when there are no hits and number of best "
					"hits desired is 1.\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHitsDelete();
	}

	/* Number of biggest clusters is 1, there is only 1 hit, and number
	 * of desired hits is only 1. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 2, 5);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHitsGetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 1 || refIdx_bestMatch[0] != 0 || shift_bestMatch[0] != 2
				|| refPos_bestMatch[0] != 5)
			fail("Incorrect behavior when number of biggest clusters is 1, "
					"there is only 1 hit, and number of best hits desired is 1."
					"\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHitsDelete();
	}

	/* Number of biggest clusters is 1 and number of best hits desired is 1. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 2, 5);
		mapHitsAddHit(0, 2, 10);
		mapHitsAddHit(0, 3, 10);
		mapHitsAddHit(1, 2, 5);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHitsGetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 1 || refIdx_bestMatch[0] != 0 || shift_bestMatch[0] != 2
				|| refPos_bestMatch[0] != 5)
			fail("Incorrect behavior when number of biggest clusters is 1 "
					"and number of best hits desired is 1.\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHitsDelete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 1. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 2, 5);
		mapHitsAddHit(0, 2, 10);
		mapHitsAddHit(1, 3, 2);
		mapHitsAddHit(1, 3, 5);
		mapHitsAddHit(2, 3, 5);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHitsGetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 1 || (refIdx_bestMatch[0] != 0
				&& refIdx_bestMatch[0] != 1) || (shift_bestMatch[0] != 2
						&& shift_bestMatch[0] != 3) || (refPos_bestMatch[0] != 5
								&& refPos_bestMatch[0] != 2))
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and number of best hits desired is 1.\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHitsDelete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 2. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 2, 5);
		mapHitsAddHit(0, 2, 10);
		mapHitsAddHit(1, 3, 2);
		mapHitsAddHit(1, 3, 5);
		mapHitsAddHit(2, 3, 5);
		int numBestHits = 2;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHitsGetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 2 || refIdx_bestMatch[0] != 0
				|| refIdx_bestMatch[1] != 1 || shift_bestMatch[0] != 2
				|| shift_bestMatch[1] != 3 || refPos_bestMatch[0] != 5
				|| refPos_bestMatch[1] != 2)
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and number of best hits desired is 2.\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHitsDelete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 3. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 2, 5);
		mapHitsAddHit(0, 2, 10);
		mapHitsAddHit(1, 3, 2);
		mapHitsAddHit(1, 3, 5);
		mapHitsAddHit(2, 3, 5);
		int numBestHits = 3;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHitsGetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 2 || refIdx_bestMatch[0] != 0
				|| refIdx_bestMatch[1] != 1 || shift_bestMatch[0] != 2
				|| shift_bestMatch[1] != 3 || refPos_bestMatch[0] != 5
				|| refPos_bestMatch[1] != 2)
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and number of best hits desired is 3.\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHitsDelete();
	}
}
END_TEST


/**
 * Tests @a swap function.
 */
START_TEST(swap)
{
	int seedLen = 8;
	mapHitsCreate(seedLen);
	mapHitsAddHit(0, 2, 5);
	mapHitsAddHit(0, 2, 10);
	mapHitsAddHit(0, 3, 10);
	mapHitsAddHit(1, 3, 5);
	char *refIdx;
	int *shift, *refPos;
	swap_wrap(1, 3, &refIdx, &shift, &refPos);
	if (refIdx[1] != 1 || refIdx[3] != 0 || shift[1] != 3 || shift[3] != 2
			|| refPos[1] != 5 || refPos[3] != 10)
		fail("Incorrect behavior.\n");
	mapHitsDelete();
}
END_TEST


/**
 * Tests @a quickSort function.
 */
START_TEST(quickSort)
{
	/* When there is 1 hit. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(3, 0, 10);
		char *refIdx;
		int *shift, *refPos;
		quickSort_wrap(&refIdx, &shift, &refPos);
		if (refIdx[0] != 3 || shift[0] != 0 || refPos[0] != 10)
			fail("Incorrect behavior when there is 1 hit.\n");
		mapHitsDelete();
	}

	/* When there are 4 hits. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(3, 0, 10);
		mapHitsAddHit(1, 3, 5);
		mapHitsAddHit(0, 2, 10);
		mapHitsAddHit(0, 2, 5);
		char *refIdx;
		int *shift, *refPos;
		quickSort_wrap(&refIdx, &shift, &refPos);
		if (refIdx[0] != 0 || refIdx[1] != 0 || refIdx[2] != 1
				|| refIdx[3] != 3 || shift[0] != 2 || shift[1] != 2
				|| shift[2] != 3 || shift[3] != 0 || refPos[0] != 5
				|| refPos[1] != 10 || refPos[2] != 5 || refPos[3] != 10)
			fail("Incorrect behavior when there are 4 hits.\n");
		mapHitsDelete();
	}

	/* When there are 2 hits. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(3, 0, 10);
		mapHitsAddHit(0, 1, 20);
		char *refIdx;
		int *shift, *refPos;
		quickSort_wrap(&refIdx, &shift, &refPos);
		if (refIdx[0] != 0 || refIdx[1] != 3 || shift[0] != 1 || shift[1] != 0
				|| refPos[0] != 20 || refPos[1] != 10)
			fail("Incorrect behavior when there are 2 hits.\n");
		mapHitsDelete();
	}

	/* When there are 3 hits. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(3, 0, 10);
		mapHitsAddHit(0, 1, 20);
		mapHitsAddHit(3, 0, 30);
		char *refIdx;
		int *shift, *refPos;
		quickSort_wrap(&refIdx, &shift, &refPos);
		if (refIdx[0] != 0 || refIdx[1] != 3 || refIdx[2] != 3 || shift[0] != 1
				|| shift[1] != 0 || shift[2] != 0 || refPos[0] != 20
				|| refPos[1] != 10 || refPos[2] != 30)
			fail("Incorrect behavior when there are 3 hits.\n");
		mapHitsDelete();
	}
}
END_TEST


/**
 * Tests @a createClusters function.
 */
START_TEST(createClusters)
{
	/* There is 1 hit and 1 cluster. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		int *cluster;
		createClusters_wrap(&cluster);
		if (cluster[0] != 0)
			fail("Incorrect behavior when there is 1 hit and 1 cluster.\n");
		mapHitsDelete();
	}

	/* There are 2 hits and 1 cluster. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		mapHitsAddHit(0, 1, 25);
		int *cluster;
		createClusters_wrap(&cluster);
		if (cluster[0] != 0)
			fail("Incorrect behavior when there are 2 hits and 1 cluster.\n");
		mapHitsDelete();
	}

	/* There are 3 hits and 2 clusters. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		mapHitsAddHit(3, 0, 10);
		mapHitsAddHit(3, 0, 30);
		int *cluster;
		createClusters_wrap(&cluster);
		if (cluster[0] != 0 || cluster[1] != 1)
			fail("Incorrect behavior when there are 3 hits and 2 clusters.\n");
		mapHitsDelete();
	}

	/* There are 5 hits and 3 clusters. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		mapHitsAddHit(3, 0, 10);
		mapHitsAddHit(3, 0, 30);
		mapHitsAddHit(5, -1, 0);
		mapHitsAddHit(5, -1, 8);
		int *cluster;
		createClusters_wrap(&cluster);
		if (cluster[0] != 0 || cluster[1] != 1 || cluster[2] != 3)
			fail("Incorrect behavior when there are 5 hits and 3 clusters.\n");
		mapHitsDelete();
	}
}
END_TEST


/**
 * Tests @a findBiggestClusterSize function.
 */
START_TEST(findBiggestClusterSize)
{
	/* Biggest cluster size is 1 and there is 1 such cluster. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		int biggestClusterSize = findBiggestClusterSize_wrap();
		if (biggestClusterSize != 1)
			fail("Incorrect behavior when the biggest cluster size is 1 "
					"and there is 1 such cluster.\n");
		mapHitsDelete();
	}

	/* Biggest cluster size is 2 and there is 1 such cluster. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		mapHitsAddHit(3, 0, 10);
		mapHitsAddHit(3, 0, 30);
		int biggestClusterSize = findBiggestClusterSize_wrap();
		if (biggestClusterSize != 2)
			fail("Incorrect behavior when the biggest cluster size is 2 "
					"and there is 1 such cluster.\n");
		mapHitsDelete();
	}

	/* Biggest cluster size is 2 and there are 2 such clusters. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		mapHitsAddHit(3, 0, 10);
		mapHitsAddHit(3, 0, 30);
		mapHitsAddHit(5, -1, 0);
		mapHitsAddHit(5, -1, 8);
		int biggestClusterSize = findBiggestClusterSize_wrap();
		if (biggestClusterSize != 2)
			fail("Incorrect behavior when the biggest cluster size is 2 "
					"and there are 2 such clusters.\n");
		mapHitsDelete();
	}

	/* Biggest cluster size is 3 and there are 2 such clusters. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		mapHitsAddHit(3, 0, 10);
		mapHitsAddHit(3, 0, 30);
		mapHitsAddHit(3, 0, 40);
		mapHitsAddHit(3, 2, 40);
		mapHitsAddHit(5, -1, 0);
		mapHitsAddHit(5, -1, 8);
		mapHitsAddHit(5, -1, 15);
		mapHitsAddHit(5, 0, 0);
		int biggestClusterSize = findBiggestClusterSize_wrap();
		if (biggestClusterSize != 3)
			fail("Incorrect behavior when the biggest cluster size is 3 "
					"and there are 2 such clusters.\n");
		mapHitsDelete();
	}
}
END_TEST


/**
 * Tests @a findBiggestClusters function.
 */
START_TEST(findBiggestClusters)
{
	/* Number of biggest clusters is 2 and the number of desired clusters
	 * is 1. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != numBestHits || refIdx_bestMatch[0] != 0
				|| shift_bestMatch[0] != 1 || refPos_bestMatch[0] != 20)
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and the number of desired clusters is 1.\n");
		mapHitsDelete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
	}

	/* Number of biggest clusters is 2 and the number of desired clusters
	 * is 1. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		mapHitsAddHit(3, 0, 10);
		mapHitsAddHit(3, 0, 30);
		mapHitsAddHit(5, -1, 0);
		mapHitsAddHit(5, -1, 8);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != numBestHits || (refIdx_bestMatch[0] != 3
				&& refIdx_bestMatch[0] != 5) || (shift_bestMatch[0] != 0
				&& shift_bestMatch[0] != -1) || (refPos_bestMatch[0] != 10
				&& refPos_bestMatch[0] != 0))
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and the number of desired clusters is 1.\n");
		mapHitsDelete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
	}

	/* Number of biggest clusters is 2 and the number of desired clusters
	 * is 2. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		mapHitsAddHit(3, 0, 10);
		mapHitsAddHit(3, 0, 30);
		mapHitsAddHit(5, -1, 0);
		mapHitsAddHit(5, -1, 8);
		int numBestHits = 2;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != numBestHits || refIdx_bestMatch[0] != 3
				|| refIdx_bestMatch[1] != 5 || shift_bestMatch[0] != 0
				|| shift_bestMatch[1] != -1 || refPos_bestMatch[0] != 10
				|| refPos_bestMatch[1] != 0)
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and the number of desired clusters is 2.\n");
		mapHitsDelete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
	}

	/* Number of biggest clusters is 1 and the number of desired clusters
	 * is 3. */
	{
		int seedLen = 8;
		mapHitsCreate(seedLen);
		mapHitsAddHit(0, 1, 20);
		mapHitsAddHit(3, -1, 8);
		mapHitsAddHit(3, 0, 10);
		mapHitsAddHit(3, 0, 30);
		mapHitsAddHit(3, 0, 40);
		int numBestHits = 3;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != 1 || refIdx_bestMatch[0] != 3
				|| shift_bestMatch[0] != 0 || refPos_bestMatch[0] != 10)
			fail("Incorrect behavior when number of biggest clusters is 1 "
					"and the number of desired clusters is 3.\n");
		mapHitsDelete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *mapHitsSuite(void)
{
	Suite *s = suite_create("mapHits");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, mapHitsAddHit);
	tcase_add_test(testCaseCore, swap);
	tcase_add_test(testCaseCore, quickSort);
	tcase_add_test(testCaseCore, createClusters);
	tcase_add_test(testCaseCore, findBiggestClusterSize);
	tcase_add_test(testCaseCore, findBiggestClusters);
	tcase_add_test(testCaseCore, mapHitsGetBestHits);
	suite_add_tcase (s, testCaseCore);

	return s;
}
