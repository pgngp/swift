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
#include "../src/refNameMap.h"
#include "testRefNameMap.h"
#include <signal.h>
#include <stdio.h>
#include "../src/common.h"


/**
 * Handles OS signals.
 */
static void signalHandler(int param)
{
	if (param == SIGSEGV)
		fail("Received segmentation fault.\n");
}


/**
 * Tests @a refNameMapCreate function.
 */
START_TEST(refNameMapCreate)
{
	/* Reference file name is NULL. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char *refFile = NULL;
		int *namePos = NULL;
		int *seqStartPos = NULL;
		int *seqEndPos = NULL;
		int size = 0;
		refNameMapCreate_wrap(refFile, &namePos, &seqStartPos, &seqEndPos, &size);
		if (size != 0 || namePos != NULL || seqStartPos != NULL
				|| seqEndPos != NULL)
			fail("Incorrect behavior when reference file name is NULL.\n");
		refNameMapDelete();
		remove(refFile);
	}

	/* 1 reference. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT", filePtr);
		fclose(filePtr);

		int *namePos = NULL;
		int *seqStartPos = NULL;
		int *seqEndPos = NULL;
		int size = 0;
		refNameMapCreate_wrap(refFile, &namePos, &seqStartPos, &seqEndPos, &size);
		if (size != 1 || namePos[0] != 0 || seqStartPos[0] != 6
				|| seqEndPos[0] != 45)
			fail("Incorrect behavior when there is only 1 reference.\n");
		refNameMapDelete();
		remove(refFile);
	}

	/* 3 references. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAA\n>ref2\nACGTACGTAA\n>ref3\nACGTACGTAA",
				filePtr);
		fclose(filePtr);

		int *namePos = NULL;
		int *seqStartPos = NULL;
		int *seqEndPos = NULL;
		int size = 0;
		refNameMapCreate_wrap(refFile, &namePos, &seqStartPos, &seqEndPos, &size);
		if (size != 3 || namePos[0] != 0 || seqStartPos[0] != 6
				|| seqEndPos[0] != 16 || namePos[1] != 17
				|| seqStartPos[1] != 23 || seqEndPos[1] != 33
				|| namePos[2] != 34 || seqStartPos[2] != 40
				|| seqEndPos[2] != 49)
			fail("Incorrect behavior when there are 3 references.\n");
		refNameMapDelete();
		remove(refFile);
	}
}
END_TEST


/**
 * Tests @a refNameMapGetNameOffset function.
 */
START_TEST(refNameMapGetNameOffsetFunc)
{
	/* 1 reference. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT", filePtr);
		fclose(filePtr);

		refNameMapCreate(refFile);
		int nameOffset = refNameMapGetNameOffset(16);
		if (nameOffset != 0)
			fail("Incorrect behavior when there is only 1 reference (case 1)."
					"\n");

		nameOffset = refNameMapGetNameOffset(6);
		if (nameOffset != 0)
			fail("Incorrect behavior when there is only 1 reference (case 2)."
					"\n");

		nameOffset = refNameMapGetNameOffset(45);
		if (nameOffset != 0)
			fail("Incorrect behavior when there is only 1 reference (case 3)."
					"\n");

		nameOffset = refNameMapGetNameOffset(47);
		if (nameOffset != -1)
			fail("Incorrect behavior when there is only 1 reference (case 4)."
					"\n");

		nameOffset = refNameMapGetNameOffset(UINT_MAX);
		if (nameOffset != -1)
			fail("Incorrect behavior when there is only 1 reference (case 5)."
					"\n");

		nameOffset = refNameMapGetNameOffset(100);
		if (nameOffset != -1)
			fail("Incorrect behavior when there is only 1 reference (case 6)."
					"\n");

		refNameMapDelete();
		remove(refFile);
	}

	/* 3 references. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_test.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAA\n>ref2\nACGTACGTAA\n>ref3\nACGTACGTAA",
				filePtr);
		fclose(filePtr);

		refNameMapCreate(refFile);
		int nameOffset = refNameMapGetNameOffset(16);
		if (nameOffset != 0)
			fail("Incorrect behavior when there are 3 references (case 1)."
					"\n");

		nameOffset = refNameMapGetNameOffset(6);
		if (nameOffset != 0)
			fail("Incorrect behavior when there are 3 references (case 2)."
					"\n");

		nameOffset = refNameMapGetNameOffset(5);
		if (nameOffset != -1)
			fail("Incorrect behavior when there are 3 references (case 3)."
					"\n");

		nameOffset = refNameMapGetNameOffset(2);
		if (nameOffset != -1)
			fail("Incorrect behavior when there are 3 references (case 4)."
					"\n");

		nameOffset = refNameMapGetNameOffset(0);
		if (nameOffset != -1)
			fail("Incorrect behavior when there are 3 references (case 5)."
					"\n");

		nameOffset = refNameMapGetNameOffset(23);
		if (nameOffset != 17)
			fail("Incorrect behavior when there are 3 references (case 6)."
					"\n");

		nameOffset = refNameMapGetNameOffset(33);
		if (nameOffset != 17)
			fail("Incorrect behavior when there are 3 references (case 7)."
					"\n");

		nameOffset = refNameMapGetNameOffset(27);
		if (nameOffset != 17)
			fail("Incorrect behavior when there are 3 references (case 8)."
					"\n");

		nameOffset = refNameMapGetNameOffset(22);
		if (nameOffset != -1)
			fail("Incorrect behavior when there are 3 references (case 9)."
					"\n");

		nameOffset = refNameMapGetNameOffset(20);
		if (nameOffset != -1)
			fail("Incorrect behavior when there are 3 references (case 10)."
					"\n");

		nameOffset = refNameMapGetNameOffset(40);
		if (nameOffset != 34)
			fail("Incorrect behavior when there are 3 references (case 11)."
					"\n");

		nameOffset = refNameMapGetNameOffset(49);
		if (nameOffset != 34)
			fail("Incorrect behavior when there are 3 references (case 12)."
					"\n");

		nameOffset = refNameMapGetNameOffset(50);
		if (nameOffset != -1)
			fail("Incorrect behavior when there are 3 references (case 13)."
					"\n");

		nameOffset = refNameMapGetNameOffset(34);
		if (nameOffset != -1)
			fail("Incorrect behavior when there are 3 references (case 14)."
					"\n");

		nameOffset = refNameMapGetNameOffset(39);
		if (nameOffset != -1)
			fail("Incorrect behavior when there are 3 references (case 15)."
					"\n");

		nameOffset = refNameMapGetNameOffset(35);
		if (nameOffset != -1)
			fail("Incorrect behavior when there are 3 references (case 16)."
					"\n");

		nameOffset = refNameMapGetNameOffset(110);
		if (nameOffset != -1)
			fail("Incorrect behavior when there are 3 references (case 17)."
					"\n");

		nameOffset = refNameMapGetNameOffset(UINT_MAX);
		if (nameOffset != -1)
			fail("Incorrect behavior when there are 3 references (case 18)."
					"\n");

		refNameMapDelete();
		remove(refFile);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *refNameMapSuite(void)
{
	Suite *s = suite_create("refNameMap");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, refNameMapCreate);
	tcase_add_test(testCaseCore, refNameMapGetNameOffsetFunc);
	suite_add_tcase(s, testCaseCore);

	return s;
}
