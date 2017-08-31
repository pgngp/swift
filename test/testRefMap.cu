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
#include "testRefMap.h"
#include "../src/refMap.h"
#include "../src/common.h"


/**
 * Tests @a refMapCreate function.
 */
START_TEST(refMapCreate)
{
	/* 1 references. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGT\n", filePtr);
		fclose(filePtr);

		uint *nameOffset, *seqOffset;
		int *numBasesPerLine, size;
		refMapCreate_wrap(refFile, &nameOffset, &seqOffset, &numBasesPerLine,
				&size);

		if (size != 1 || nameOffset[0] != 0 || seqOffset[0] != 6
				|| numBasesPerLine[0] != 0)
			fail("Incorrect behavior when there is 1 reference.\n");
		refMapFree();
	}

	/* 3 references, all of different lengths. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGT\n\n>ref2\nACGTACGTA\nACGT\n>ref3\nACGTACGTAA\n"
				"ACGTACGTAA\nACGT", filePtr);
		fclose(filePtr);

		uint *nameOffset, *seqOffset;
		int *numBasesPerLine, size;
		refMapCreate_wrap(refFile, &nameOffset, &seqOffset, &numBasesPerLine,
				&size);

		if (size != 3 || nameOffset[0] != 0 || nameOffset[1] != 16
				|| nameOffset[2] != 37 || seqOffset[0] != 6
				|| seqOffset[1] != 22 || seqOffset[2] != 43
				|| numBasesPerLine[0] != 0 || numBasesPerLine[1] != 9
				|| numBasesPerLine[2] != 10)
			fail("Incorrect behavior when there are 3 references, all of "
					"different lengths.\n");
		refMapFree();
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
		int index = binarySearch_wrapper(8, arr, arrSize);
		if (index != 0)
			fail("Incorrect behavior when element to be searched is the "
					"first element in the array.\n");
	}

	/* Element to be searched is the last element in the array. */
	{
		uint arr[] = {8, 10, 31};
		int arrSize = 3;
		int index = binarySearch_wrapper(31, arr, arrSize);
		if (index != 2)
			fail("Incorrect behavior when element to be searched is the "
					"last element in the array.\n");
	}

	/* Element to be searched is the middle element in the array. */
	{
		uint arr[] = {8, 10, 31};
		int arrSize = 3;
		int index = binarySearch_wrapper(10, arr, arrSize);
		if (index != 1)
			fail("Incorrect behavior when element to be searched is the "
					"middle element in the array.\n");
	}

	/* Element to be searched is not the first, middle or last element. */
	{
		uint arr[] = {8, 10, 31, 44};
		int arrSize = 4;
		int index = binarySearch_wrapper(31, arr, arrSize);
		if (index != 2)
			fail("Incorrect behavior when element to be searched is not the "
					"first, middle, or last element in the array.\n");
	}

	/* Element to be searched is not in the array. */
	{
		uint arr[] = {8, 10, 31};
		int arrSize = 3;
		int index = binarySearch_wrapper(5, arr, arrSize);
		if (index != -1)
			fail("Incorrect behavior when element to be searched is not in the "
					"array.\n");
	}
}
END_TEST


/**
 * Tests @a refMapGetFileOffset function.
 */
START_TEST(refMapGetFileOffset)
{
	/* 1 reference. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGT\n", filePtr);
		fclose(filePtr);

		refMapCreate(refFile);

		uint fileOffset = refMapGetFileOffset(0, 0);
		if (fileOffset != 6)
			fail("Incorrect behavior when there is 1 reference (case 1).\n");

		fileOffset = refMapGetFileOffset(0, 3);
		if (fileOffset != 9)
			fail("Incorrect behavior when there is 1 reference (case 2).\n");

		fileOffset = refMapGetFileOffset(0, 7);
		if (fileOffset != 13)
			fail("Incorrect behavior when there is 1 reference (case 3).\n");

		refMapFree();
	}

	/* 3 references, all of different length. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGT\n\n>ref2\nACGTACGTA\nACGT\n>ref3\nACGTACGTAA\n"
				"ACGTACGTAA\nACGT", filePtr);
		fclose(filePtr);

		refMapCreate(refFile);

		uint fileOffset = refMapGetFileOffset(0, 0);
		if (fileOffset != 6)
			fail("Incorrect behavior when there are 3 references, all"
					"of different length (case 1).\n");

		fileOffset = refMapGetFileOffset(0, 3);
		if (fileOffset != 9)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 2).\n");

		fileOffset = refMapGetFileOffset(0, 7);
		if (fileOffset != 13)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 3).\n");

		fileOffset = refMapGetFileOffset(16, 0);
		if (fileOffset != 22)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 4).\n");

		fileOffset = refMapGetFileOffset(16, 8);
		if (fileOffset != 30)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 5).\n");

		fileOffset = refMapGetFileOffset(16, 9);
		if (fileOffset != 32)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 6).\n");

		fileOffset = refMapGetFileOffset(16, 12);
		if (fileOffset != 35)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 7).\n");

		fileOffset = refMapGetFileOffset(37, 0);
		if (fileOffset != 43)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 8).\n");

		fileOffset = refMapGetFileOffset(37, 9);
		if (fileOffset != 52)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 9).\n");

		fileOffset = refMapGetFileOffset(37, 10);
		if (fileOffset != 54)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 10).\n");

		fileOffset = refMapGetFileOffset(37, 19);
		if (fileOffset != 63)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 11).\n");

		fileOffset = refMapGetFileOffset(37, 20);
		if (fileOffset != 65)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 12).\n");

		fileOffset = refMapGetFileOffset(37, 23);
		if (fileOffset != 68)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 13).\n");

		refMapFree();
	}
}
END_TEST


/**
 * Tests @a refMapGetOffsets function.
 */
START_TEST(refMapGetOffsets)
{
	/* 1 reference. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGT\n", filePtr);
		fclose(filePtr);

		refMapCreate(refFile);
		uint nameOffset, seqOffset;
		refMapGetOffsets(0, 0, &nameOffset, &seqOffset);
		if (nameOffset != 0 || seqOffset != 6)
			fail("Incorrect behavior when there is 1 reference (case 1).\n");

		refMapGetOffsets(0, 4, &nameOffset, &seqOffset);
		if (nameOffset != 0 || seqOffset != 10)
			fail("Incorrect behavior when there is 1 reference (case 2).\n");

		refMapGetOffsets(0, 7, &nameOffset, &seqOffset);
		if (nameOffset != 0 || seqOffset != 13)
			fail("Incorrect behavior when there is 1 reference (case 3).\n");

		refMapFree();
	}


	/* 3 references, all of different length. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGT\n\n>ref2\nACGTACGTA\nACGT\n>ref3\nACGTACGTAA\n"
				"ACGTACGTAA\nACGT", filePtr);
		fclose(filePtr);

		refMapCreate(refFile);

		uint nameOffset, seqOffset;
		refMapGetOffsets(0, 0, &nameOffset, &seqOffset);
		if (nameOffset != 0 || seqOffset != 6)
			fail("Incorrect behavior when there are 3 references, all"
					"of different length (case 1).\n");

		refMapGetOffsets(0, 3, &nameOffset, &seqOffset);
		if (nameOffset != 0 || seqOffset != 9)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 2).\n");

		refMapGetOffsets(0, 7, &nameOffset, &seqOffset);
		if (nameOffset != 0 || seqOffset != 13)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 3).\n");

		refMapGetOffsets(1, 0, &nameOffset, &seqOffset);
		if (nameOffset != 16 || seqOffset != 22)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 4).\n");

		refMapGetOffsets(1, 8, &nameOffset, &seqOffset);
		if (nameOffset != 16 || seqOffset != 30)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 5).\n");

		refMapGetOffsets(1, 9, &nameOffset, &seqOffset);
		if (nameOffset != 16 || seqOffset != 32)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 6).\n");

		refMapGetOffsets(1, 12, &nameOffset, &seqOffset);
		if (nameOffset != 16 || seqOffset != 35)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 7).\n");

		refMapGetOffsets(2, 0, &nameOffset, &seqOffset);
		if (nameOffset != 37 || seqOffset != 43)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 8).\n");

		refMapGetOffsets(2, 9, &nameOffset, &seqOffset);
		if (nameOffset != 37 || seqOffset != 52)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 9).\n");

		refMapGetOffsets(2, 10, &nameOffset, &seqOffset);
		if (nameOffset!= 37 || seqOffset != 54)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 10).\n");

		refMapGetOffsets(2, 19, &nameOffset, &seqOffset);
		if (nameOffset != 37 || seqOffset != 63)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 11).\n");

		refMapGetOffsets(2, 20, &nameOffset, &seqOffset);
		if (nameOffset != 37 || seqOffset != 65)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 12).\n");

		refMapGetOffsets(2, 23, &nameOffset, &seqOffset);
		if (nameOffset != 37 || seqOffset != 68)
			fail("Incorrect behavior when there are 3 references, all "
					"of different length (case 13).\n");

		refMapFree();
	}
}
END_TEST


/**
 * Tests @a refMapGetIndexAndPos function.
 */
START_TEST(refMapGetIndexAndPos)
{
	/* 1 reference. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGT\n", filePtr);
		fclose(filePtr);

		refMapCreate(refFile);
		uint seqOffset;
		int refIdx, pos;

		seqOffset = 6;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 0 || pos != 0)
			fail("Incorrect behavior when there is 1 reference (case 1).\n");

		seqOffset = 9;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 0 || pos != 3)
			fail("Incorrect behavior when there is 1 reference (case 2).\n");

		seqOffset = 13;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 0 || pos != 7)
			fail("Incorrect behavior when there is 1 reference (case 3).\n");

		refMapFree();
	}

	/* 3 references, all of different length. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGT\n\n>ref2\nACGTACGTA\nACGT\n>ref3\nACGTACGTAA\n"
				"ACGTACGTAA\nACGT", filePtr);
		fclose(filePtr);

		refMapCreate(refFile);
		uint seqOffset;
		int refIdx, pos;

		seqOffset = 6;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 0 || pos != 0)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 1).\n");

		seqOffset = 10;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 0 || pos != 4)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 2).\n");

		seqOffset = 13;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 0 || pos != 7)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 3).\n");

		seqOffset = 22;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 1 || pos != 0)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 4).\n");

		seqOffset = 26;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 1 || pos != 4)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 5).\n");

		seqOffset = 30;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 1 || pos != 8)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 6).\n");

		seqOffset = 32;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 1 || pos != 9)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 7).\n");

		seqOffset = 35;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 1 || pos != 12)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 8).\n");

		seqOffset = 43;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 2 || pos != 0)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 9).\n");

		seqOffset = 47;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 2 || pos != 4)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 10).\n");

		seqOffset = 52;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 2 || pos != 9)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 11).\n");

		seqOffset = 54;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 2 || pos != 10)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 12).\n");

		seqOffset = 58;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 2 || pos != 14)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 13).\n");

		seqOffset = 63;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 2 || pos != 19)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 14).\n");

		seqOffset = 65;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 2 || pos != 20)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 15).\n");

		seqOffset = 68;
		refMapGetIndexAndPos(seqOffset, &refIdx, &pos);
		if (refIdx != 2 || pos != 23)
			fail("Incorrect behavior when there are 3 references, all of "
					"different length (case 16).\n");

		refMapFree();
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *refMapSuite(void)
{
	Suite *s = suite_create("refMap");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, refMapCreate);
	tcase_add_test(testCaseCore, binarySearch);
	tcase_add_test(testCaseCore, refMapGetFileOffset);
	tcase_add_test(testCaseCore, refMapGetOffsets);
	tcase_add_test(testCaseCore, refMapGetIndexAndPos);
	suite_add_tcase (s, testCaseCore);

	return s;
}
