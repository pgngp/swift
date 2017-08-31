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
#include "../src/mapHits4.h"
#include "testMapHits4.h"


/**
 * Tests @a mapHits4AddHit function.
 */
START_TEST(mapHits4AddHit)
{
	/* One hit is added. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		char refIdx = 1;
		int shift = 2;
		int refPos = 5;
		ulong *hits;
		int size;
		mapHits4AddHit_wrap(&hits, &size, refIdx, shift, refPos);
		if (hits[0] != 72057594574798853 || size != 1)
			fail("Incorrect behavior when a hit is added.\n");
		mapHits4Delete();
	}
}
END_TEST


/**
 * Tests @a mapHits4GetBestHits function.
 */
START_TEST(mapHits4GetBestHits)
{
	/* There are no hits and number of desired hits is only 1. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits4GetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 0)
			fail("Incorrect behavior when there are no hits and number of best "
					"hits desired is 1.\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHits4Delete();
	}

	/* Number of biggest clusters is 1, there is only 1 hit, and number
	 * of desired hits is only 1. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 2, 5);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits4GetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 1 || refIdx_bestMatch[0] != 0 || shift_bestMatch[0] != 2
				|| refPos_bestMatch[0] != 5)
			fail("Incorrect behavior when number of biggest clusters is 1, "
					"there is only 1 hit, and number of best hits desired is 1."
					"\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHits4Delete();
	}

	/* Number of biggest clusters is 1 and number of best hits desired is 1. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 2, 5);
		mapHits4AddHit(0, 2, 10);
		mapHits4AddHit(0, 3, 10);
		mapHits4AddHit(1, 2, 5);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits4GetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 1 || refIdx_bestMatch[0] != 0 || shift_bestMatch[0] != 2
				|| refPos_bestMatch[0] != 5)
			fail("Incorrect behavior when number of biggest clusters is 1 "
					"and number of best hits desired is 1.\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHits4Delete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 1. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 2, 5);
		mapHits4AddHit(0, 2, 10);
		mapHits4AddHit(1, 3, 2);
		mapHits4AddHit(1, 3, 5);
		mapHits4AddHit(2, 3, 5);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits4GetBestHits(numBestHits, refIdx_bestMatch,
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
		mapHits4Delete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 2. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 2, 5);
		mapHits4AddHit(0, 2, 10);
		mapHits4AddHit(1, 3, 2);
		mapHits4AddHit(1, 3, 5);
		mapHits4AddHit(2, 3, 5);
		int numBestHits = 2;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits4GetBestHits(numBestHits, refIdx_bestMatch,
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
		mapHits4Delete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 3. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 2, 5);
		mapHits4AddHit(0, 2, 10);
		mapHits4AddHit(1, 3, 2);
		mapHits4AddHit(1, 3, 5);
		mapHits4AddHit(2, 3, 5);
		int numBestHits = 3;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits4GetBestHits(numBestHits, refIdx_bestMatch,
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
		mapHits4Delete();
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
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		int *cluster;
		createClusters_wrap4(&cluster);
		if (cluster[0] != 0)
			fail("Incorrect behavior when there is 1 hit and 1 cluster.\n");
		mapHits4Delete();
	}

	/* There are 2 hits and 1 cluster. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		mapHits4AddHit(0, 1, 25);
		int *cluster;
		createClusters_wrap4(&cluster);
		if (cluster[0] != 0)
			fail("Incorrect behavior when there are 2 hits and 1 cluster.\n");
		mapHits4Delete();
	}

	/* There are 3 hits and 2 clusters. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		mapHits4AddHit(3, 0, 10);
		mapHits4AddHit(3, 0, 30);
		int *cluster;
		createClusters_wrap4(&cluster);
		if (cluster[0] != 0 || cluster[1] != 1)
			fail("Incorrect behavior when there are 3 hits and 2 clusters.\n");
		mapHits4Delete();
	}

	/* There are 5 hits and 3 clusters. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		mapHits4AddHit(3, 0, 10);
		mapHits4AddHit(3, 0, 30);
		mapHits4AddHit(5, -1, 0);
		mapHits4AddHit(5, -1, 8);
		int *cluster;
		createClusters_wrap4(&cluster);
		if (cluster[0] != 0 || cluster[1] != 1 || cluster[2] != 3)
			fail("Incorrect behavior when there are 5 hits and 3 clusters.\n");
		mapHits4Delete();
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
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		int biggestClusterSize = findBiggestClusterSize_wrap4();
		if (biggestClusterSize != 1)
			fail("Incorrect behavior when the biggest cluster size is 1 "
					"and there is 1 such cluster.\n");
		mapHits4Delete();
	}

	/* Biggest cluster size is 2 and there is 1 such cluster. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		mapHits4AddHit(3, 0, 10);
		mapHits4AddHit(3, 0, 30);
		int biggestClusterSize = findBiggestClusterSize_wrap4();
		if (biggestClusterSize != 2)
			fail("Incorrect behavior when the biggest cluster size is 2 "
					"and there is 1 such cluster.\n");
		mapHits4Delete();
	}

	/* Biggest cluster size is 2 and there are 2 such clusters. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		mapHits4AddHit(3, 0, 10);
		mapHits4AddHit(3, 0, 30);
		mapHits4AddHit(5, -1, 0);
		mapHits4AddHit(5, -1, 8);
		int biggestClusterSize = findBiggestClusterSize_wrap4();
		if (biggestClusterSize != 2)
			fail("Incorrect behavior when the biggest cluster size is 2 "
					"and there are 2 such clusters.\n");
		mapHits4Delete();
	}

	/* Biggest cluster size is 3 and there are 2 such clusters. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		mapHits4AddHit(3, 0, 10);
		mapHits4AddHit(3, 0, 30);
		mapHits4AddHit(3, 0, 40);
		mapHits4AddHit(3, 2, 40);
		mapHits4AddHit(5, -1, 0);
		mapHits4AddHit(5, -1, 8);
		mapHits4AddHit(5, -1, 15);
		mapHits4AddHit(5, 0, 0);
		int biggestClusterSize = findBiggestClusterSize_wrap4();
		if (biggestClusterSize != 3)
			fail("Incorrect behavior when the biggest cluster size is 3 "
					"and there are 2 such clusters.\n");
		mapHits4Delete();
	}
}
END_TEST


/**
 * Tests @a findBiggestClusters function.
 */
START_TEST(findBiggestClusters)
{
	/* Number of biggest clusters is 1 and the number of desired clusters
	 * is 1. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap4(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != numBestHits || refIdx_bestMatch[0] != 0
				|| shift_bestMatch[0] != 1 || refPos_bestMatch[0] != 20)
			fail("Incorrect behavior when number of biggest clusters is 1 "
					"and the number of desired clusters is 1.\n");
		mapHits4Delete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
	}

	/* Number of biggest clusters is 2 and the number of desired clusters
	 * is 1. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		mapHits4AddHit(3, 0, 10);
		mapHits4AddHit(3, 0, 30);
		mapHits4AddHit(5, -1, 0);
		mapHits4AddHit(5, -1, 8);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap4(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != numBestHits || (refIdx_bestMatch[0] != 3
				&& refIdx_bestMatch[0] != 5) || (shift_bestMatch[0] != 0
				&& shift_bestMatch[0] != -1) || (refPos_bestMatch[0] != 10
				&& refPos_bestMatch[0] != 0))
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and the number of desired clusters is 1.\n");
		mapHits4Delete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
	}

	/* Number of biggest clusters is 2 and the number of desired clusters
	 * is 2. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		mapHits4AddHit(3, 0, 10);
		mapHits4AddHit(3, 0, 30);
		mapHits4AddHit(5, -1, 0);
		mapHits4AddHit(5, -1, 8);
		int numBestHits = 2;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap4(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != numBestHits || refIdx_bestMatch[0] != 3
				|| refIdx_bestMatch[1] != 5 || shift_bestMatch[0] != 0
				|| shift_bestMatch[1] != -1 || refPos_bestMatch[0] != 10
				|| refPos_bestMatch[1] != 0)
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and the number of desired clusters is 2.\n");
		mapHits4Delete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
	}

	/* Number of biggest clusters is 1 and the number of desired clusters
	 * is 3. */
	{
		int seedLen = 8;
		mapHits4Create(seedLen);
		mapHits4AddHit(0, 1, 20);
		mapHits4AddHit(3, -1, 8);
		mapHits4AddHit(3, 0, 10);
		mapHits4AddHit(3, 0, 30);
		mapHits4AddHit(3, 0, 40);
		int numBestHits = 3;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap4(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != 1 || refIdx_bestMatch[0] != 3
				|| shift_bestMatch[0] != 0 || refPos_bestMatch[0] != 10)
			fail("Incorrect behavior when number of biggest clusters is 1 "
					"and the number of desired clusters is 3.\n");
		mapHits4Delete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *mapHits4Suite(void)
{
	Suite *s = suite_create("mapHits4");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, mapHits4AddHit);
	tcase_add_test(testCaseCore, mapHits4GetBestHits);
	tcase_add_test(testCaseCore, createClusters);
	tcase_add_test(testCaseCore, findBiggestClusterSize);
	tcase_add_test(testCaseCore, findBiggestClusters);
	suite_add_tcase (s, testCaseCore);

	return s;
}
