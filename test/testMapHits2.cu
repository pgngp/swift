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
#include "../src/mapHits2.h"
#include "testMapHits2.h"


/**
 * Tests @a search function.
 */
START_TEST(search)
{
	/* Reference index and 'shift' are not already present in the list
	 * (case 1). */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		char refIdx = 1;
		int shift = 0;
		int nextHighestIdx = 0;
		int index = search_wrap(refIdx, shift, &nextHighestIdx);

		if (index != -1 || nextHighestIdx != 0)
			fail("Incorrect behavior when reference index and 'shift' are "
					"not already present in the list (case 1).\n");
		mapHits2Delete();
	}

	/* Reference index and 'shift' are not already present in the list
	 * (case 2). */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		char refIdx = 1;
		int shift = 0;
		int refPos = 5;
		mapHits2AddHit(refIdx, shift, refPos);
		refIdx = 1;
		shift = -1;
		int nextHighestIdx = 0;
		int index = search_wrap(refIdx, shift, &nextHighestIdx);

		if (index != -1 || nextHighestIdx != 0)
			fail("Incorrect behavior when reference index and 'shift' are "
					"not already present in the list (case 2).\n");
		mapHits2Delete();
	}

	/* Reference index and 'shift' are already present in the list
	 * (case 1). */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		char refIdx = 1;
		int shift = 0;
		int refPos = 5;
		mapHits2AddHit(refIdx, shift, refPos);
		int nextHighestIdx = 0;
		int index = search_wrap(refIdx, shift, &nextHighestIdx);
		if (index != 0)
			fail("Incorrect behavior when reference index and 'shift' are "
					"already present in the list (case 1).\n");
		mapHits2Delete();
	}

	/* Reference index and 'shift' are already present in the list
	 * (case 2). */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		char refIdx = 1;
		int shift = 0;
		int refPos = 5;
		mapHits2AddHit(refIdx, shift, refPos);
		refIdx = 1, shift = -2, refPos = 5;
		mapHits2AddHit(refIdx, shift, refPos);
		int nextHighestIdx = 0;
		int index = search_wrap(refIdx, shift, &nextHighestIdx);
		if (index != 1)
			fail("Incorrect behavior when reference index and 'shift' are "
					"already present in the list (case 2).\n");
		mapHits2Delete();
	}
}
END_TEST


/**
 * Tests @a mapHits2AddHit function.
 */
START_TEST(mapHits2AddHit)
{
	/* The 'hit' to be added is not already in there. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		char refIdx = 1;
		int shift = 2;
		int refPos = 5;
		char *refIdxArr;
		int *shiftArr, *refPosArr, *numClusterHits, size;
		mapHits2AddHit_wrap(&refIdxArr, &shiftArr, &refPosArr, &numClusterHits,
				&size, refIdx, shift, refPos);
		if (refIdxArr[0] != 1 || shiftArr[0] != 2 || refPosArr[0] != 5
				|| numClusterHits[0] != 1 || size != 1)
			fail("Incorrect behavior when the 'hit' to be added is not "
					"already in there.\n");
		mapHits2Delete();
	}

	/* The 'hit' to be added is already in there; reference position is
	 * greater than the reference position of the incumbent. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		char refIdx = 1;
		int shift = 2;
		int refPos = 5;
		char *refIdxArr;
		int *shiftArr, *refPosArr, *numClusterHits, size;
		mapHits2AddHit(refIdx, shift, refPos);
		refPos = 6;
		mapHits2AddHit_wrap(&refIdxArr, &shiftArr, &refPosArr, &numClusterHits,
				&size, refIdx, shift, refPos);
		if (refIdxArr[0] != 1 || shiftArr[0] != 2 || refPosArr[0] != 5
				|| numClusterHits[0] != 2 || size != 1)
			fail("Incorrect behavior when the 'hit' to be added is already "
					"in there; reference position is greater than the "
					"reference position of the incumbent.\n");
		mapHits2Delete();
	}

	/* The 'hit' to be added is already in there; reference position is
	 * less than the reference position of the incumbent. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		char refIdx = 1;
		int shift = 2;
		int refPos = 5;
		char *refIdxArr;
		int *shiftArr, *refPosArr, *numClusterHits, size;
		mapHits2AddHit(refIdx, shift, refPos);
		refPos = 4;
		mapHits2AddHit_wrap(&refIdxArr, &shiftArr, &refPosArr, &numClusterHits,
				&size, refIdx, shift, refPos);
		if (refIdxArr[0] != 1 || shiftArr[0] != 2 || refPosArr[0] != 4
				|| numClusterHits[0] != 2 || size != 1)
			fail("Incorrect behavior when the 'hit' to be added is already "
					"in there; reference position is less than the "
					"reference position of the incumbent.\n");
		mapHits2Delete();
	}
}
END_TEST


/**
 * Tests @a mapHits2AddHit function.
 */
START_TEST(mapHits2AddHit2)
{
	/* The 'hit' to be added is not already in there. */
	{
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp1.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);

		int seedLen = 8;
		mapHits2Create2(refFile, seedLen);
		char refIdx = 0;
		int qryPos = 2;
		int refPos = 8;
		char *ref;
		uchar *hits;
		int *pos, biggestClusterSize;
		mapHits2AddHit2_wrap(&ref, &hits, &pos, &biggestClusterSize, refIdx,
				qryPos, refPos);
		int index = (refPos + 1) / seedLen;
		index = index - ((qryPos + 1) / seedLen) + MAX_QRY_SEQ_LENGTH;
		fprintf(stderr, "index = %d, ", index);
		fprintf(stderr, "index = %d, ref[%d] = %d, hits[%d] = %d, "
				"pos[%d] = %d, biggestClusterSize = %d\n", index, index,
				ref[index], index, hits[index], index, pos[index],
				biggestClusterSize);
		if (ref[index] != 0 || hits[index] != 1 || pos[index] != 8
				|| biggestClusterSize != 1)
			fail("Incorrect behavior when the 'hit' to be added is not "
					"already in there.\n");
		mapHits2Delete2();
		remove(refFile);
	}

	/* The 'hit' to be added is already in there; reference position is
	 * greater than the reference position of the incumbent. */
	{
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp1.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);

		int seedLen = 8;
		mapHits2Create2(refFile, seedLen);
		char refIdx = 0;
		int qryPos = 2;
		int refPos = 0;
		mapHits2AddHit2(refIdx, qryPos, refPos);
		qryPos = 10;
		refPos = 8;
		char *ref;
		uchar *hits;
		int *pos, biggestClusterSize;
		mapHits2AddHit2_wrap(&ref, &hits, &pos, &biggestClusterSize, refIdx,
				qryPos, refPos);
		int index = (refPos + 1) / seedLen;
		index = index - qryPos + MAX_QRY_SEQ_LENGTH;
		fprintf(stderr, "index = %d, ref[199] = %d, hits[199] = %d, "
				"pos[199] = %d, biggestClusterSize = %d\n", index, ref[199],
				hits[199], pos[199], biggestClusterSize);
		if (ref[index] != 0 || hits[index] != 2 || pos[index] != 0
				|| biggestClusterSize != 2)
			fail("Incorrect behavior when the 'hit' to be added is already "
					"in there; reference position is greater than the "
					"reference position of the incumbent.\n");
		mapHits2Delete2();
		remove(refFile);
	}

//	/* The 'hit' to be added is already in there; reference position is
//	 * less than the reference position of the incumbent. */
//	{
//		int seedLen = 8;
//		mapHits2Create(seedLen);
//		char refIdx = 1;
//		int shift = 2;
//		int refPos = 5;
//		char *refIdxArr;
//		int *shiftArr, *refPosArr, *numClusterHits, size;
//		mapHits2AddHit(refIdx, shift, refPos);
//		refPos = 4;
//		mapHits2AddHit_wrap(&refIdxArr, &shiftArr, &refPosArr, &numClusterHits,
//				&size, refIdx, shift, refPos);
//		if (refIdxArr[0] != 1 || shiftArr[0] != 2 || refPosArr[0] != 4
//				|| numClusterHits[0] != 2 || size != 1)
//			fail("Incorrect behavior when the 'hit' to be added is already "
//					"in there; reference position is less than the "
//					"reference position of the incumbent.\n");
//		mapHits2Delete();
//	}
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
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 1, 20);
		int biggestClusterSize = findBiggestClusterSize_wrap2();
		if (biggestClusterSize != 1)
			fail("Incorrect behavior when the biggest cluster size is 1 "
					"and there is 1 such cluster.\n");
		mapHits2Delete();
	}

	/* Biggest cluster size is 2 and there is 1 such cluster. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 1, 20);
		mapHits2AddHit(3, 0, 10);
		mapHits2AddHit(3, 0, 30);
		int biggestClusterSize = findBiggestClusterSize_wrap2();
		if (biggestClusterSize != 2)
			fail("Incorrect behavior when the biggest cluster size is 2 "
					"and there is 1 such cluster.\n");
		mapHits2Delete();
	}

	/* Biggest cluster size is 2 and there are 2 such clusters. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 1, 20);
		mapHits2AddHit(3, 0, 10);
		mapHits2AddHit(3, 0, 30);
		mapHits2AddHit(5, -1, 0);
		mapHits2AddHit(5, -1, 8);
		int biggestClusterSize = findBiggestClusterSize_wrap2();
		if (biggestClusterSize != 2)
			fail("Incorrect behavior when the biggest cluster size is 2 "
					"and there are 2 such clusters.\n");
		mapHits2Delete();
	}

	/* Biggest cluster size is 3 and there are 2 such clusters. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 1, 20);
		mapHits2AddHit(3, 0, 10);
		mapHits2AddHit(3, 0, 30);
		mapHits2AddHit(3, 0, 40);
		mapHits2AddHit(3, 2, 40);
		mapHits2AddHit(5, -1, 0);
		mapHits2AddHit(5, -1, 8);
		mapHits2AddHit(5, -1, 15);
		mapHits2AddHit(5, 0, 0);
		int biggestClusterSize = findBiggestClusterSize_wrap2();
		if (biggestClusterSize != 3)
			fail("Incorrect behavior when the biggest cluster size is 3 "
					"and there are 2 such clusters.\n");
		mapHits2Delete();
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
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 1, 20);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap2(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != numBestHits || refIdx_bestMatch[0] != 0
				|| shift_bestMatch[0] != 1 || refPos_bestMatch[0] != 20)
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and the number of desired clusters is 1.\n");
		mapHits2Delete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
	}

	/* Number of biggest clusters is 2 and the number of desired clusters
	 * is 1. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 1, 20);
		mapHits2AddHit(3, 0, 10);
		mapHits2AddHit(3, 0, 30);
		mapHits2AddHit(5, -1, 0);
		mapHits2AddHit(5, -1, 8);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap2(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != numBestHits || (refIdx_bestMatch[0] != 3
				&& refIdx_bestMatch[0] != 5) || (shift_bestMatch[0] != 0
				&& shift_bestMatch[0] != -1) || (refPos_bestMatch[0] != 10
				&& refPos_bestMatch[0] != 0))
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and the number of desired clusters is 1.\n");
		mapHits2Delete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
	}

	/* Number of biggest clusters is 2 and the number of desired clusters
	 * is 2. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 1, 20);
		mapHits2AddHit(3, 0, 10);
		mapHits2AddHit(3, 0, 30);
		mapHits2AddHit(5, -1, 0);
		mapHits2AddHit(5, -1, 8);
		int numBestHits = 2;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap2(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != numBestHits || refIdx_bestMatch[0] != 3
				|| refIdx_bestMatch[1] != 5 || shift_bestMatch[0] != 0
				|| shift_bestMatch[1] != -1 || refPos_bestMatch[0] != 10
				|| refPos_bestMatch[1] != 0)
			fail("Incorrect behavior when number of biggest clusters is 2 "
					"and the number of desired clusters is 2.\n");
		mapHits2Delete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
	}

	/* Number of biggest clusters is 1 and the number of desired clusters
	 * is 3. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 1, 20);
		mapHits2AddHit(3, -1, 8);
		mapHits2AddHit(3, 0, 10);
		mapHits2AddHit(3, 0, 30);
		mapHits2AddHit(3, 0, 40);
		int numBestHits = 3;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numBiggestClusters = findBiggestClusters_wrap2(numBestHits,
				refIdx_bestMatch, shift_bestMatch, refPos_bestMatch);
		if (numBiggestClusters != 1 || refIdx_bestMatch[0] != 3
				|| shift_bestMatch[0] != 0 || refPos_bestMatch[0] != 10)
			fail("Incorrect behavior when number of biggest clusters is 1 "
					"and the number of desired clusters is 3.\n");
		mapHits2Delete();
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
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
		mapHits2Create(seedLen);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits2GetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 0)
			fail("Incorrect behavior when there are no hits and number of best "
					"hits desired is 1.\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHits2Delete();
	}

	/* Number of biggest clusters is 1, there is only 1 hit, and number
	 * of desired hits is only 1. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 2, 5);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits2GetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 1 || refIdx_bestMatch[0] != 0 || shift_bestMatch[0] != 2
				|| refPos_bestMatch[0] != 5)
			fail("Incorrect behavior when number of biggest clusters is 1, "
					"there is only 1 hit, and number of best hits desired is 1."
					"\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHits2Delete();
	}

	/* Number of biggest clusters is 1 and number of best hits desired is 1. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 2, 5);
		mapHits2AddHit(0, 2, 10);
		mapHits2AddHit(0, 3, 10);
		mapHits2AddHit(1, 2, 5);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits2GetBestHits(numBestHits, refIdx_bestMatch,
				shift_bestMatch, refPos_bestMatch);
		if (numHits != 1 || refIdx_bestMatch[0] != 0 || shift_bestMatch[0] != 2
				|| refPos_bestMatch[0] != 5)
			fail("Incorrect behavior when number of biggest clusters is 1 "
					"and number of best hits desired is 1.\n");
		free(refIdx_bestMatch);
		free(shift_bestMatch);
		free(refPos_bestMatch);
		mapHits2Delete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 1. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 2, 5);
		mapHits2AddHit(0, 2, 10);
		mapHits2AddHit(1, 3, 2);
		mapHits2AddHit(1, 3, 5);
		mapHits2AddHit(2, 3, 5);
		int numBestHits = 1;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits2GetBestHits(numBestHits, refIdx_bestMatch,
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
		mapHits2Delete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 2. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 2, 5);
		mapHits2AddHit(0, 2, 10);
		mapHits2AddHit(1, 3, 2);
		mapHits2AddHit(1, 3, 5);
		mapHits2AddHit(2, 3, 5);
		int numBestHits = 2;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits2GetBestHits(numBestHits, refIdx_bestMatch,
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
		mapHits2Delete();
	}

	/* Number of biggest clusters is 2 and number of best hits desired is 3. */
	{
		int seedLen = 8;
		mapHits2Create(seedLen);
		mapHits2AddHit(0, 2, 5);
		mapHits2AddHit(0, 2, 10);
		mapHits2AddHit(1, 3, 2);
		mapHits2AddHit(1, 3, 5);
		mapHits2AddHit(2, 3, 5);
		int numBestHits = 3;
		char *refIdx_bestMatch = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestMatch = (int *) calloc(numBestHits, sizeof(int));
		int numHits = mapHits2GetBestHits(numBestHits, refIdx_bestMatch,
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
		mapHits2Delete();
	}
}
END_TEST


/**
 * Tests @a createRefMap function.
 */
START_TEST(createRefMap)
{
	/* 1 reference. */
	{
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp1.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);
		int seedLen = 8;
		mapHits2Create2(refFile, seedLen);
		int *lastPos;
		createRefMap_wrap(&lastPos);
		if (lastPos[0] != 0 || lastPos[1] != 2)
			fail("Incorrect behavior when there is only 1 reference.\n");
		mapHits2Delete2();
		remove(refFile);
	}

	/* 2 references (case 1). */
	{
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp1.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n>ref2\nACGTACGTACGT\n", filePtr);
		fclose(filePtr);
		int seedLen = 8;
		mapHits2Create2(refFile, seedLen);
		int *lastPos;
		createRefMap_wrap(&lastPos);
		if (lastPos[0] != 0 || lastPos[1] != 2 || lastPos[2] != 3)
			fail("Incorrect behavior when there are 2 references (case 1).\n");
		mapHits2Delete2();
		remove(refFile);
	}

	/* 2 references (case 2). */
	{
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp1.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n>ref2\nACGTACGT\nACGTACGT\n", filePtr);
		fclose(filePtr);
		int seedLen = 8;
		mapHits2Create2(refFile, seedLen);
		int *lastPos;
		createRefMap_wrap(&lastPos);
		if (lastPos[0] != 0 || lastPos[1] != 2 || lastPos[2] != 4)
			fail("Incorrect behavior when there are 2 references (case 2).\n");
		mapHits2Delete2();
		remove(refFile);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *mapHits2Suite(void)
{
	Suite *s = suite_create("mapHits2");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_set_timeout(testCaseCore, 0);
	tcase_add_test(testCaseCore, search);
	tcase_add_test(testCaseCore, mapHits2AddHit);
	tcase_add_test(testCaseCore, mapHits2AddHit2);
	tcase_add_test(testCaseCore, findBiggestClusterSize);
	tcase_add_test(testCaseCore, createRefMap);
//	tcase_add_test(testCaseCore, findBiggestClusters);
//	tcase_add_test(testCaseCore, mapHitsGetBestHits);
	suite_add_tcase (s, testCaseCore);

	return s;
}
