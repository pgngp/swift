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
#include "testRefPosMap.h"
#include "../src/refPosMap.h"
#include <signal.h>
#include "../src/common.h"
#include "../src/reference.h"


/**
 * Handles OS signals.
 */
static void signalHandler(int param)
{
	if (param == SIGSEGV)
		fail("Received segmentation fault.\n");
}


/**
 * Tests @a refMapCreate function.
 */
START_TEST(refMapCreate)
{
	/* Hashtable is NULL. */
	{
		(void) signal(SIGSEGV, signalHandler);

		RefPosList **hashtable = NULL;
		int numKeys = 0;
		int numVals = 0;
		int seedLen = 8;
		uint *localPosArr = NULL;
		uint *globalPosArr = NULL;
		int arrSize;
		refMapCreate_wrap(hashtable, numKeys, numVals, seedLen, &localPosArr,
				&globalPosArr, &arrSize);
		if (localPosArr != NULL || globalPosArr != NULL || arrSize != 0)
			fail("Incorrect behavior when hash table is NULL.\n");
		refMapDelete();
	}

	/* There is 1 reference, 5 reference tuples and 5 positions. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT", filePtr);
		fclose(filePtr);
		int seedLen = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLen);
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLen);
		int numVals = 5;
		uint *localPosArr = NULL;
		uint *globalPosArr = NULL;
		int arrSize;
		refMapCreate_wrap(hashTable, numKeys, numVals, seedLen, &localPosArr,
				&globalPosArr, &arrSize);
		if (arrSize != numVals || localPosArr[0] != 0 || globalPosArr[0] != 6
				|| localPosArr[1] != 1 || globalPosArr[1] != 14
				|| localPosArr[2] != 2 || globalPosArr[2] != 22
				|| localPosArr[3] != 3 || globalPosArr[3] != 30
				|| localPosArr[4] != 4 || globalPosArr[4] != 38)
			fail("Incorrect behavior when there is 1 reference, 5 reference "
					"tuples and 5 positions.\n");
		refMapDelete();
		refDeleteHashTable(hashTable, numKeys);
		remove(refFile);
	}

	/* There are 2 references, 9 reference tuples and 9 positions. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT\n"
				">ref2\nCCGTACGTGATAACGTTCGTATTTGGGTACAC", filePtr);
		fclose(filePtr);
		int seedLen = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLen);
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLen);
		int numVals = 9;
		uint *localPosArr = NULL;
		uint *globalPosArr = NULL;
		int arrSize;
		refMapCreate_wrap(hashTable, numKeys, numVals, seedLen, &localPosArr,
				&globalPosArr, &arrSize);
		if (arrSize != numVals || localPosArr[0] != 0 || globalPosArr[0] != 6
				|| localPosArr[1] != 1 || globalPosArr[1] != 14
				|| localPosArr[2] != 2 || globalPosArr[2] != 22
				|| localPosArr[3] != 3 || globalPosArr[3] != 30
				|| localPosArr[4] != 4 || globalPosArr[4] != 38
				|| localPosArr[5] != 134217728 || globalPosArr[5] != 53
				|| localPosArr[6] != 134217729 || globalPosArr[6] != 61
				|| localPosArr[7] != 134217730 || globalPosArr[7] != 69
				|| localPosArr[8] != 134217731 || globalPosArr[8] != 77)
			fail("Incorrect behavior when there are 2 references, 9 reference "
					"tuples and 9 positions.\n");
		refMapDelete();
		refDeleteHashTable(hashTable, numKeys);
		remove(refFile);
	}

	/* There are 2 references, 4 reference tuples and 8 positions. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGT\n"
				">ref2\nACGTACGTAATAACGTACGTATTTACGTACACCCGT", filePtr);
		fclose(filePtr);
		int seedLen = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLen);
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLen);
		int numVals = 8;
		uint *localPosArr = NULL;
		uint *globalPosArr = NULL;
		int arrSize;
		refMapCreate_wrap(hashTable, numKeys, numVals, seedLen, &localPosArr,
				&globalPosArr, &arrSize);
		if (arrSize != numVals || localPosArr[0] != 0 || globalPosArr[0] != 6
				|| localPosArr[1] != 1 || globalPosArr[1] != 14
				|| localPosArr[2] != 2 || globalPosArr[2] != 22
				|| localPosArr[3] != 3 || globalPosArr[3] != 30
				|| localPosArr[4] != 134217728 || globalPosArr[4] != 49
				|| localPosArr[5] != 134217729 || globalPosArr[5] != 57
				|| localPosArr[6] != 134217730 || globalPosArr[6] != 65
				|| localPosArr[7] != 134217731 || globalPosArr[7] != 73)
			fail("Incorrect behavior when there are 2 references, 4 reference "
					"tuples and 8 positions.\n");
		refMapDelete();
		refDeleteHashTable(hashTable, numKeys);
		remove(refFile);
	}
}
END_TEST


/**
 * Tests @a refMapGetGlobalPos function.
 */
START_TEST(refMapGetGlobalPosFunc)
{
	/* Case 1. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT", filePtr);
		fclose(filePtr);
		int seedLen = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLen);
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLen);
		int numVals = 5;
		refMapCreate(hashTable, numKeys, numVals, seedLen);

		uint globalPos = refMapGetGlobalPos(0);
		if (globalPos != 6)
			fail ("Incorrect behavior (case 1.1).\n");

		globalPos = refMapGetGlobalPos(1);
		if (globalPos != 14)
			fail ("Incorrect behavior (case 1.2).\n");

		globalPos = refMapGetGlobalPos(2);
		if (globalPos != 22)
			fail ("Incorrect behavior (case 1.3).\n");

		globalPos = refMapGetGlobalPos(3);
		if (globalPos != 30)
			fail ("Incorrect behavior (case 1.4).\n");

		globalPos = refMapGetGlobalPos(4);
		if (globalPos != 38)
			fail ("Incorrect behavior (case 1.5).\n");

		globalPos = refMapGetGlobalPos(10);
		if (globalPos != UINT_MAX)
			fail ("Incorrect behavior (case 1.7).\n");

		globalPos = refMapGetGlobalPos(5);
		if (globalPos != UINT_MAX)
			fail ("Incorrect behavior (case 1.8).\n");

		refMapDelete();
		refDeleteHashTable(hashTable, numKeys);
		remove(refFile);
	}
}
END_TEST


/**
 * Tests @a binarySearch function.
 */
START_TEST(binarySearch)
{
	/* Element to be searched is the first element in the array. */
	{
		uint arr[] = {8, 10, 31};
		int arrSize = 3;
		int index = binarySearch_wrap(8, arr, arrSize);
		if (index != 0)
			fail("Incorrect behavior when element to be searched is the "
					"first element in the array.\n");
	}

	/* Element to be searched is the last element in the array. */
	{
		uint arr[] = {8, 10, 31};
		int arrSize = 3;
		int index = binarySearch_wrap(31, arr, arrSize);
		if (index != 2)
			fail("Incorrect behavior when element to be searched is the "
					"last element in the array.\n");
	}

	/* Element to be searched is the middle element in the array. */
	{
		uint arr[] = {8, 10, 31};
		int arrSize = 3;
		int index = binarySearch_wrap(10, arr, arrSize);
		if (index != 1)
			fail("Incorrect behavior when element to be searched is the "
					"middle element in the array.\n");
	}

	/* Element to be searched is not the first, middle or last element. */
	{
		uint arr[] = {8, 10, 31, 44};
		int arrSize = 4;
		int index = binarySearch_wrap(31, arr, arrSize);
		if (index != 2)
			fail("Incorrect behavior when element to be searched is not the "
					"first, middle, or last element in the array.\n");
	}

	/* Element to be searched is not in the array. */
	{
		uint arr[] = {8, 10, 31};
		int arrSize = 3;
		int index = binarySearch_wrap(5, arr, arrSize);
		if (index != -1)
			fail("Incorrect behavior when element to be searched is not in the "
					"array.\n");
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *refPosMapSuite(void)
{
	Suite *s = suite_create("refPosMap");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, refMapCreate);
	tcase_add_test(testCaseCore, refMapGetGlobalPosFunc);
	tcase_add_test(testCaseCore, binarySearch);
	suite_add_tcase(s, testCaseCore);

	return s;
}
