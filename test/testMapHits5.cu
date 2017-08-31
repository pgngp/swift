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
#include "../src/mapHits5.h"
#include "testMapHits5.h"
#include "../src/common.h"


/**
 * Tests @a mapHits5AddHit function.
 */
START_TEST(mapHits5AddHit)
{
	/* One hit is added. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		char refIdx = 1;
		int refPos = 5;
		char *refIdxArr;
		int *refPosArr, size;
		mapHits5AddHit_wrap(&refIdxArr, &refPosArr, &size, refIdx, refPos);
		if (refIdxArr[0] != 1 || refPosArr[0] != 5 || size != 1)
			fail("Incorrect behavior when a hit is added.\n");
		mapHits5Delete();
	}
}
END_TEST


/**
 * Tests @a swap function.
 */
START_TEST(swap)
{
	int seedLen = 8;
	mapHits5Create(seedLen);
	mapHits5AddHit(0, 5);
	mapHits5AddHit(0, 10);
	mapHits5AddHit(0, 10);
	mapHits5AddHit(1, 5);
	char *refIdx;
	int *refPos;
	swap_wrap5(1, 3, &refIdx, &refPos);
	if (refIdx[1] != 1 || refIdx[3] != 0 || refPos[1] != 5 || refPos[3] != 10)
		fail("Incorrect behavior.\n");
	mapHits5Delete();
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
		mapHits5Create(seedLen);
		mapHits5AddHit(3, 10);
		char *refIdx;
		int *refPos;
		quickSort_wrap5(&refIdx, &refPos);
		if (refIdx[0] != 3 || refPos[0] != 10)
			fail("Incorrect behavior when there is 1 hit.\n");
		mapHits5Delete();
	}

	/* When there are 4 hits. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(3, 10);
		mapHits5AddHit(1, 5);
		mapHits5AddHit(0, 10);
		mapHits5AddHit(0, 5);
		char *refIdx;
		int *refPos;
		quickSort_wrap5(&refIdx, &refPos);
		if (refIdx[0] != 0 || refIdx[1] != 0 || refIdx[2] != 1
				|| refIdx[3] != 3 || refPos[0] != 5 || refPos[1] != 10
				|| refPos[2] != 5 || refPos[3] != 10)
			fail("Incorrect behavior when there are 4 hits.\n");
		mapHits5Delete();
	}

	/* When there are 2 hits. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(3, 10);
		mapHits5AddHit(0, 20);
		char *refIdx;
		int *refPos;
		quickSort_wrap5(&refIdx, &refPos);
		if (refIdx[0] != 0 || refIdx[1] != 3 || refPos[0] != 20
				|| refPos[1] != 10)
			fail("Incorrect behavior when there are 2 hits.\n");
		mapHits5Delete();
	}

	/* When there are 3 hits. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(3, 10);
		mapHits5AddHit(0, 20);
		mapHits5AddHit(3, 30);
		char *refIdx;
		int *refPos;
		quickSort_wrap5(&refIdx, &refPos);
		if (refIdx[0] != 0 || refIdx[1] != 3 || refIdx[2] != 3
				|| refPos[0] != 20 || refPos[1] != 10 || refPos[2] != 30)
			fail("Incorrect behavior when there are 3 hits.\n");
		mapHits5Delete();
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
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		int *cluster;
		createClusters_wrap5(&cluster);
		if (cluster[0] != 0)
			fail("Incorrect behavior when there is 1 hit and 1 cluster.\n");
		mapHits5Delete();
	}

	/* There are 2 hits and 1 cluster. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		mapHits5AddHit(0, 25);
		int *cluster;
		createClusters_wrap5(&cluster);
		if (cluster[0] != 0)
			fail("Incorrect behavior when there are 2 hits and 1 cluster.\n");
		mapHits5Delete();
	}

	/* There are 3 hits and 2 clusters. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		mapHits5AddHit(3, 10);
		mapHits5AddHit(3, 30);
		int *cluster;
		createClusters_wrap5(&cluster);
		if (cluster[0] != 0 || cluster[1] != 1 || cluster[2] != 0)
			fail("Incorrect behavior when there are 3 hits and 2 clusters.\n");
		mapHits5Delete();
	}

	/* There are 5 hits and 3 clusters. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		mapHits5AddHit(3, 10);
		mapHits5AddHit(3, 30);
		mapHits5AddHit(5, 0);
		mapHits5AddHit(5, 8);
		int *cluster;
		createClusters_wrap5(&cluster);
		if (cluster[0] != 0 || cluster[1] != 1 || cluster[2] != 3
				|| cluster[3] != 0)
			fail("Incorrect behavior when there are 5 hits and 3 clusters.\n");
		mapHits5Delete();
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
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		int biggestClusterSize = findBiggestClusterSize_wrap5();
		if (biggestClusterSize != 1)
			fail("Incorrect behavior when the biggest cluster size is 1 "
					"and there is 1 such cluster.\n");
		mapHits5Delete();
	}

	/* Biggest cluster size is 2 and there is 1 such cluster. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		mapHits5AddHit(3, 10);
		mapHits5AddHit(3, 30);
		int biggestClusterSize = findBiggestClusterSize_wrap5();
		if (biggestClusterSize != 2)
			fail("Incorrect behavior when the biggest cluster size is 2 "
					"and there is 1 such cluster.\n");
		mapHits5Delete();
	}

	/* Biggest cluster size is 2 and there are 2 such clusters. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		mapHits5AddHit(3, 10);
		mapHits5AddHit(3, 30);
		mapHits5AddHit(5, 0);
		mapHits5AddHit(5, 8);
		int biggestClusterSize = findBiggestClusterSize_wrap5();
		if (biggestClusterSize != 2)
			fail("Incorrect behavior when the biggest cluster size is 2 "
					"and there are 2 such clusters.\n");
		mapHits5Delete();
	}

	/* Biggest cluster size is 4 and there are 2 such clusters. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		mapHits5AddHit(3, 10);
		mapHits5AddHit(3, 30);
		mapHits5AddHit(3, 40);
		mapHits5AddHit(3, 40);
		mapHits5AddHit(5, 0);
		mapHits5AddHit(5, 8);
		mapHits5AddHit(5, 15);
		mapHits5AddHit(5, 0);
		int biggestClusterSize = findBiggestClusterSize_wrap5();
		if (biggestClusterSize != 4)
			fail("Incorrect behavior when the biggest cluster size is 4 "
					"and there are 2 such clusters.\n");
		mapHits5Delete();
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
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		int numBestHits = 1;
		char *refIdx = (char *) calloc(numBestHits, sizeof(char));
		int *refPos = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap5(numBestHits,
				refIdx, refPos);
		if (numBiggestClusters != numBestHits || refIdx[0] != 0
				|| refPos[0] != 20)
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and the number of desired clusters is 1.\n");
		mapHits5Delete();
		free(refIdx);
		free(refPos);
	}

	/* Number of biggest clusters is 2 and the number of desired clusters
	 * is 1. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		mapHits5AddHit(3, 10);
		mapHits5AddHit(3, 30);
		mapHits5AddHit(5, 0);
		mapHits5AddHit(5, 8);
		int numBestHits = 1;
		char *refIdx = (char *) calloc(numBestHits, sizeof(char));
		int *refPos = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap5(numBestHits,
				refIdx, refPos);
		if (numBiggestClusters != numBestHits || (refIdx[0] != 3
				&& refIdx[0] != 5) || (refPos[0] != 10 && refPos[0] != 0))
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and the number of desired clusters is 1.\n");
		mapHits5Delete();
		free(refIdx);
		free(refPos);
	}

	/* Number of biggest clusters is 2 and the number of desired clusters
	 * is 2. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		mapHits5AddHit(3, 10);
		mapHits5AddHit(3, 30);
		mapHits5AddHit(5, 0);
		mapHits5AddHit(5, 8);
		int numBestHits = 2;
		char *refIdx = (char *) calloc(numBestHits, sizeof(char));
		int *refPos = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap5(numBestHits,
				refIdx, refPos);
		if (numBiggestClusters != numBestHits || refIdx[0] != 3
				|| refIdx[1] != 5 || refPos[0] != 10 || refPos[1] != 0)
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and the number of desired clusters is 2.\n");
		mapHits5Delete();
		free(refIdx);
		free(refPos);
	}

	/* Number of biggest clusters is 1 and the number of desired clusters
	 * is 3. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 20);
		mapHits5AddHit(2, 8);
		mapHits5AddHit(3, 10);
		mapHits5AddHit(3, 30);
		mapHits5AddHit(3, 40);
		int numBestHits = 3;
		char *refIdx = (char *) calloc(numBestHits, sizeof(char));
		int *refPos = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap5(numBestHits,
				refIdx, refPos);
		if (numBiggestClusters != 1 || refIdx[0] != 3 || refPos[0] != 10)
			fail("Incorrect behavior when number of biggest clusters is 1 "
					"and the number of desired clusters is 3.\n");
		mapHits5Delete();
		free(refIdx);
		free(refPos);
	}
}
END_TEST


/**
 * Tests @a mapHitsGetBestHits function.
 */
START_TEST(mapHits5GetBestHits)
{
	/* There are no hits and number of desired hits is only 1. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		int numBestHits = 1;
		char *refIdx = (char *) calloc(numBestHits, sizeof(char));
		int *refPos = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits5GetBestHits(numBestHits, refIdx, refPos);
		if (numHits != 0)
			fail("Incorrect behavior when there are no hits and number of best "
					"hits desired is 1.\n");
		free(refIdx);
		free(refPos);
		mapHits5Delete();
	}

	/* Number of biggest clusters is 1, there is only 1 hit, and number
	 * of desired hits is only 1. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 5);
		int numBestHits = 1;
		char *refIdx = (char *) calloc(numBestHits, sizeof(char));
		int *refPos = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits5GetBestHits(numBestHits, refIdx, refPos);
		if (numHits != 1 || refIdx[0] != 0 || refPos[0] != 5)
			fail("Incorrect behavior when number of biggest clusters is 1, "
					"there is only 1 hit, and number of best hits desired is 1."
					"\n");
		free(refIdx);
		free(refPos);
		mapHits5Delete();
	}

	/* Number of biggest clusters is 1 and number of best hits desired is 1. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 5);
		mapHits5AddHit(0, 10);
		mapHits5AddHit(0, 10);
		mapHits5AddHit(1, 5);
		int numBestHits = 1;
		char *refIdx = (char *) calloc(numBestHits, sizeof(char));
		int *refPos = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits5GetBestHits(numBestHits, refIdx, refPos);
		if (numHits != 1 || refIdx[0] != 0 || refPos[0] != 5)
			fail("Incorrect behavior when number of biggest clusters is 1 "
					"and number of best hits desired is 1.\n");
		free(refIdx);
		free(refPos);
		mapHits5Delete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 1. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 5);
		mapHits5AddHit(0, 10);
		mapHits5AddHit(1, 2);
		mapHits5AddHit(1, 5);
		mapHits5AddHit(2, 5);
		int numBestHits = 1;
		char *refIdx = (char *) calloc(numBestHits, sizeof(char));
		int *refPos = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits5GetBestHits(numBestHits, refIdx, refPos);
		if (numHits != 1 || (refIdx[0] != 0 && refIdx[0] != 1)
				|| (refPos[0] != 5 && refPos[0] != 2))
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and number of best hits desired is 1.\n");
		free(refIdx);
		free(refPos);
		mapHits5Delete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 2. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 5);
		mapHits5AddHit(0, 10);
		mapHits5AddHit(1, 2);
		mapHits5AddHit(1, 5);
		mapHits5AddHit(2, 5);
		int numBestHits = 2;
		char *refIdx = (char *) calloc(numBestHits, sizeof(char));
		int *refPos = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits5GetBestHits(numBestHits, refIdx, refPos);
		if (numHits != 2 || refIdx[0] != 0 || refIdx[1] != 1 || refPos[0] != 5
				|| refPos[1] != 2)
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and number of best hits desired is 2.\n");
		free(refIdx);
		free(refPos);
		mapHits5Delete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 3. */
	{
		int seedLen = 8;
		mapHits5Create(seedLen);
		mapHits5AddHit(0, 5);
		mapHits5AddHit(0, 10);
		mapHits5AddHit(1, 2);
		mapHits5AddHit(1, 5);
		mapHits5AddHit(2, 5);
		int numBestHits = 3;
		char *refIdx = (char *) calloc(numBestHits, sizeof(char));
		int *refPos = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits5GetBestHits(numBestHits, refIdx, refPos);
		if (numHits != 2 || refIdx[0] != 0 || refIdx[1] != 1 || refPos[0] != 5
				|| refPos[1] != 2)
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and number of best hits desired is 3.\n");
		free(refIdx);
		free(refPos);
		mapHits5Delete();
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *mapHits5Suite(void)
{
	Suite *s = suite_create("mapHits5");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, mapHits5AddHit);
	tcase_add_test(testCaseCore, swap);
	tcase_add_test(testCaseCore, quickSort);
	tcase_add_test(testCaseCore, createClusters);
	tcase_add_test(testCaseCore, findBiggestClusterSize);
	tcase_add_test(testCaseCore, findBiggestClusters);
	tcase_add_test(testCaseCore, mapHits5GetBestHits);
	suite_add_tcase (s, testCaseCore);

	return s;
}
