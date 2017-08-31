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

#include <stdio.h>
#include <check.h>
#include "testInput.h"
#include "../src/input.h"
#include "../src/common.h"


/**
 * Tests getOffsets function.
 */
START_TEST(getOffsets)
{
	/* Test whether the output is correct when a valid file is passed. */
	{
		char file[MAX_FILE_NAME_LENGTH];
		sprintf(file, "%s%stestInput.fa", TMP_DIR, PATH_SEPARATOR);
		FILE *filePtr = fopen(file, "w+");
		fputs("q0\t1\t1\t8\t8\t106\t112\n", filePtr);
		rewind(filePtr);
		uint queryNameOffset = 0;
		uint querySeqOffset = 0;
		uint refNameOffset = 0;
		uint refSeqOffset = 0;
		uint refSeqFragOffset = 0;
		int refDistance = 0;
		uint isReverseComplement = 0;
		int refIdx = 0;
		getOffsets(filePtr, &queryNameOffset, &querySeqOffset,
				&refNameOffset, &refSeqOffset, &refSeqFragOffset, &refDistance,
				&isReverseComplement, &refIdx);
		fclose(filePtr);
		remove(file);
		if (queryNameOffset != 106 || querySeqOffset != 112 || refDistance != 8
				|| isReverseComplement != 1 || refIdx != 1)
			fail("Incorrect behavior when a valid file is parsed.\n");
	}
}
END_TEST


/**
 * Tests @a getOffsets_paired function.
 */
START_TEST(getOffsets_paired)
{
	/* Test whether the output is correct when a valid file is passed. */
	{
		char file[MAX_FILE_NAME_LENGTH];
		sprintf(file, "%s%stestInput.fa", TMP_DIR, PATH_SEPARATOR);
		FILE *filePtr = fopen(file, "w+");
		fputs("1\tq0\t1\t2\t1\t8\t8\t106\t112\n", filePtr);
		rewind(filePtr);
		uint queryNameOffset = 0;
		uint querySeqOffset = 0;
		int mateIdx = 0;
		uint refNameOffset = 0;
		uint refSeqOffset = 0;
		uint refSeqFragOffset = 0;
		int refDistance = 0;
		uint isReverseComplement = 0;
		int refIdx = 0;
		int qryIdx = 0;
		getOffsets_paired2(filePtr, &queryNameOffset, &querySeqOffset, &mateIdx,
				&refNameOffset, &refSeqOffset, &refSeqFragOffset, &refDistance,
				&isReverseComplement, &refIdx, &qryIdx);
		fclose(filePtr);
		remove(file);
		if (queryNameOffset != 106 || querySeqOffset != 112 || mateIdx != 2
					|| refDistance != 8 || isReverseComplement != 1
					|| refIdx != 1 || qryIdx != 1)
			fail("Incorrect behavior when a valid file is parsed.\n");
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *inputSuite(void)
{
	Suite *s = suite_create("input");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, getOffsets);
	tcase_add_test(testCaseCore, getOffsets_paired);
	suite_add_tcase(s, testCaseCore);

	return s;
}

