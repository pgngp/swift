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
#include "../src/search.h"
#include "testSearch.h"
#include <signal.h>
#include "../src/refMap.h"
#include "../src/preprocess.h"


/**
 * Handles OS signals.
 */
static void signalHandler(int param)
{
	if (param == SIGSEGV)
		fail("Received segmentation fault.\n");
}


/**
 * Tests searchQuery function.
 */
START_TEST(searchQuery)
{
	/* Query file name is NULL. */
	{
		char *qryFileName = NULL;
		int numQueries = 100;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
	}

	/* Number of queries is 0. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 0;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query length is 0. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 0;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Reference file name is NULL. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 100;
		char *refFileName = NULL;
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFileName);
	}

	/* Seed length is 0. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Match file is NULL. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char *matchFile = NULL;
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query file has only 1 query and reference file has only 1 reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nGGACGTACGT", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 1;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		fputs(">r1\nAATTCCCAACGTACGT", refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		int numMatches = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		searchQueries(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine2[] = "q1\t1\t0\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 2).\n");

		if (numMatches != 2)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 3).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}

	/* Query file has only 1 query and reference file has more than
	 * 1 reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nGGACGTACGT", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 1;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		char refSeq[] = ">r1\nAATTCCCAACGTACGT\n>r2\nGACGTACGACGTACGT"
				"\n>r3\nAAAAAAAAAAAAAAAA";
		fputs(refSeq, refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		int numMatches = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		searchQueries(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine2[] = "q1\t0\t1\t-1\t0\t0\t4\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine3[] = "q1\t0\t1\t6\t8\t0\t4\n";
		if (strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 3).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine4[] = "q1\t1\t0\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine4) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 4).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine5[] = "q1\t1\t1\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 5).\n");

		if (numMatches != 5)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 6).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}

	/* Query file has more than 1 query and reference file has only 1
	 * reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		char qrySeq[] = ">q1\nGGACGTACGT\n>q2\nAATTCCCACC\n>q3\nGAATTCCCAG\n";
		fputs(qrySeq, qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 3;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		char refSeq[] = ">r1\nAATTCCCAACGTACGT\n";
		fputs(refSeq, refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine2[] = "q1\t1\t0\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine3[] = "q2\t0\t0\t0\t0\t15\t19\n";
		if (strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 3).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine4[] = "q3\t0\t0\t-1\t0\t30\t34\n";
		if (strcmp(line, tmpLine4) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 4).\n");

		if (numMatches != 4)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 5).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}

	/* Query file has more than 1 query and reference file has more than
	 * 1 reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		char qrySeq[] = ">q1\nGGACGTACGT\n>q2\nAATTCCCACC\n>q3\nGAATTCCCAG\n";
		fputs(qrySeq, qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 3;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		char refSeq[] = ">r1\nAATTCCCAACGTACGT\n>r2\nGACGTACGACGTACGT"
				"\n>r3\nAAAAAAAAAAAAAAAA";
		fputs(refSeq, refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		int numMatches = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		searchQueries(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine2[] = "q1\t0\t1\t-1\t0\t0\t4\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine3[] = "q1\t0\t1\t6\t8\t0\t4\n";
		if (strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 3).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine4[] = "q1\t1\t0\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine4) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 4).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine5[] = "q1\t1\t1\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 5).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine6[] = "q2\t0\t0\t0\t0\t15\t19\n";
		if (strcmp(line, tmpLine6) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 6).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine7[] = "q3\t0\t0\t-1\t0\t30\t34\n";
		if (strcmp(line, tmpLine7) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 7).\n");

		if (numMatches != 7)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 8).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}
}
END_TEST


/**
 * Tests @a searchQuery2 function.
 */
START_TEST(searchQuery2)
{
	/* Query file name is NULL. */
	{
		char *qryFileName = NULL;
		int numQueries = 100;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries2(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
	}

	/* Number of queries is 0. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 0;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries2(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query length is 0. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 0;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries2(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Reference file name is NULL. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 100;
		char *refFileName = NULL;
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries2(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFileName);
	}

	/* Seed length is 0. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries2(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Match file is NULL. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char *matchFile = NULL;
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries2(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query file has only 1 query and reference file has only 1 reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nGGACGTACGT", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 1;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		fputs(">r1\nAATTCCCAACGTACGT", refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		int numMatches = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		searchQueries2(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine2[] = "q1\t1\t0\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 2).\n");

		if (numMatches != 2)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 3).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}

	/* Query file has only 1 query and reference file has more than
	 * 1 reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nGGACGTACGT", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 1;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		char refSeq[] = ">r1\nAATTCCCAACGTACGT\n>r2\nGACGTACGACGTACGT"
				"\n>r3\nAAAAAAAAAAAAAAAA";
		fputs(refSeq, refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		int numMatches = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		searchQueries2(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		char tmpLine2[] = "q1\t0\t1\t-1\t0\t0\t4\n";
		char tmpLine3[] = "q1\t0\t1\t6\t8\t0\t4\n";
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 3).\n");

		char tmpLine4[] = "q1\t1\t0\t8\t8\t0\t4\n";
		char tmpLine5[] = "q1\t1\t1\t8\t8\t0\t4\n";
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine4) != 0 && strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 4).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine4) != 0 && strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 5).\n");

		if (numMatches != 5)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 6).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}

	/* Query file has more than 1 query and reference file has only 1
	 * reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		char qrySeq[] = ">q1\nGGACGTACGT\n>q2\nAATTCCCACC\n>q3\nGAATTCCCAG\n";
		fputs(qrySeq, qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 3;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		char refSeq[] = ">r1\nAATTCCCAACGTACGT\n";
		fputs(refSeq, refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries2(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine2[] = "q1\t1\t0\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine3[] = "q2\t0\t0\t0\t0\t15\t19\n";
		if (strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 3).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine4[] = "q3\t0\t0\t-1\t0\t30\t34\n";
		if (strcmp(line, tmpLine4) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 4).\n");

		printf("numMatches = %d\n", numMatches);
		if (numMatches != 4)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 5).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}

	/* Query file has more than 1 query and reference file has more than
	 * 1 reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		char qrySeq[] = ">q1\nGGACGTACGT\n>q2\nAATTCCCACC\n>q3\nGAATTCCCAG\n";
		fputs(qrySeq, qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 3;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		char refSeq[] = ">r1\nAATTCCCAACGTACGT\n>r2\nGACGTACGACGTACGT"
				"\n>r3\nAAAAAAAAAAAAAAAA";
		fputs(refSeq, refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		int numMatches = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		searchQueries2(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		char tmpLine2[] = "q1\t0\t1\t-1\t0\t0\t4\n";
		char tmpLine3[] = "q1\t0\t1\t6\t8\t0\t4\n";
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 3).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine4[] = "q1\t1\t0\t8\t8\t0\t4\n";
		char tmpLine5[] = "q1\t1\t1\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine4) != 0 && strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 4).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine4) != 0 && strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 5).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine6[] = "q2\t0\t0\t0\t0\t15\t19\n";
		if (strcmp(line, tmpLine6) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 6).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine7[] = "q3\t0\t0\t-1\t0\t30\t34\n";
		if (strcmp(line, tmpLine7) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 7).\n");

		if (numMatches != 7)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 8).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}
}
END_TEST


/**
 * Tests @a searchQuery_paired function.
 */
START_TEST(searchQuery_paired)
{
	/* Both first and second query files are NULL. */
	{
		char *qryFile1 = NULL;
		char *qryFile2 = NULL;
		int numQrs = 100;
		int qryLen = 100;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int numRefs = 10;
		int refLen = 1000;
		int seedLen = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		uint minFragSize = 2, maxFragSize = 10;
		searchQueries_paired(qryFile1, qryFile2, numQrs, qryLen, refFile,
					numRefs, refLen, seedLen, matchFile, maxNumHits,
					minFragSize, maxFragSize);
		remove(refFile);
	}

	/* Number of queries is 0. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 0;
		int qryLen = 100;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int numRefs = 10;
		int refLen = 1000;
		int seedLen = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		uint minFragSize = 2, maxFragSize = 10;
		searchQueries_paired(qryFile1, qryFile2, numQrs, qryLen, refFile,
					numRefs, refLen, seedLen, matchFile, maxNumHits,
					minFragSize, maxFragSize);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
	}

	/* Query length is 0. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 10;
		int qryLen = 0;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int numRefs = 10;
		int refLen = 1000;
		int seedLen = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		uint minFragSize = 2, maxFragSize = 10;
		searchQueries_paired(qryFile1, qryFile2, numQrs, qryLen, refFile,
					numRefs, refLen, seedLen, matchFile, maxNumHits,
					minFragSize, maxFragSize);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
	}

	/* Reference file name is NULL. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 10;
		int qryLen = 100;
		char *refFile = NULL;
		int numRefs = 10;
		int refLen = 1000;
		int seedLen = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		uint minFragSize = 2, maxFragSize = 10;
		searchQueries_paired(qryFile1, qryFile2, numQrs, qryLen, refFile,
				numRefs, refLen, seedLen, matchFile, maxNumHits, minFragSize,
				maxFragSize);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
	}

	/* Number of references is 0. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 10;
		int qryLen = 100;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int numRefs = 0;
		int refLen = 1000;
		int seedLen = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		uint minFragSize = 2, maxFragSize = 10;
		searchQueries_paired(qryFile1, qryFile2, numQrs, qryLen, refFile,
				numRefs, refLen, seedLen, matchFile, maxNumHits, minFragSize,
				maxFragSize);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
	}

	/* Reference length is 0. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 10;
		int qryLen = 100;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int numRefs = 100;
		int refLen = 0;
		int seedLen = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		uint minFragSize = 2, maxFragSize = 10;
		searchQueries_paired(qryFile1, qryFile2, numQrs, qryLen, refFile,
				numRefs, refLen, seedLen, matchFile, maxNumHits, minFragSize,
				maxFragSize);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
	}

	/* Seed length is 0. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 10;
		int qryLen = 100;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int numRefs = 100;
		int refLen = 1000;
		int seedLen = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		uint minFragSize = 2, maxFragSize = 10;
		searchQueries_paired(qryFile1, qryFile2, numQrs, qryLen, refFile,
				numRefs, refLen, seedLen, matchFile, maxNumHits, minFragSize,
				maxFragSize);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
	}

	/* Match file is NULL. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 10;
		int qryLen = 100;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int numRefs = 100;
		int refLen = 1000;
		int seedLen = 8;
		char *matchFile = NULL;
		uint maxNumHits = 10;
		uint minFragSize = 2, maxFragSize = 10;
		searchQueries_paired(qryFile1, qryFile2, numQrs, qryLen, refFile,
				numRefs, refLen, seedLen, matchFile, maxNumHits, minFragSize,
				maxFragSize);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
	}

	/* 1 query pair and 1 reference. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		FILE *qryFilePtr1 = fopen(qryFile1, "w");
		fputs(">q1_1\nACGTACGTAA", qryFilePtr1);
		fclose(qryFilePtr1);

		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		FILE *qryFilePtr2 = fopen(qryFile2, "w");
		fputs(">q1_2\nACCGGAAGGA", qryFilePtr2);
		fclose(qryFilePtr2);

		int numQrs = 1;
		int qryLen = 10;

		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFile, "w");
		/* AATTCCCA ACGTACGT AAGTATTA CCGGAAGG AATC */
		fputs(">r1\nAATTCCCAACGTACGTAAGTATTACCGGAAGGAATC", refFilePtr);
		fclose(refFilePtr);
		int numRefs = 1;
		int refLen = 36;
		int seedLen = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		uint minFragSize = 8, maxFragSize = 30;
		searchQueries_paired(qryFile1, qryFile2, numQrs, qryLen, refFile,
					numRefs, refLen, seedLen, matchFile, maxNumHits,
					minFragSize, maxFragSize);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
//		char tmpLine1[] = "q1_1\t0\t0\t8\t8\t0\t6\t0\t4\t10\n";
		char tmpLine1[] = "q1_1\t0\t0\t8\t8\t0\t6\t0\t12\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when there is 1 query pair and 1 "
					"reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
//		char tmpLine2[] = "q1_2\t0\t0\t23\t24\t0\t6\t0\t4\t25\n";
		char tmpLine2[] = "q1_2\t0\t0\t23\t24\t0\t6\t0\t28\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when there is 1 query pair and 1 "
					"reference (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
//		char tmpLine3[] = "q1_1\t1\t0\t6\t8\t0\t6\t0\t4\t8\n";
		char tmpLine3[] = "q1_1\t1\t0\t6\t8\t0\t6\t0\t12\n";
		if (strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when there is 1 query pair and 1 "
					"reference (case 3).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
//		char tmpLine4[] = "q1_2\t0\t0\t23\t24\t0\t6\t0\t4\t25\n";
		char tmpLine4[] = "q1_2\t0\t0\t23\t24\t0\t6\t0\t28\n";
		if (strcmp(line, tmpLine4) != 0)
			fail("Incorrect behavior when there is 1 query pair and 1 "
					"reference (case 4).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
	}
}
END_TEST


/**
 * Tests @a searchQuery_paired function.
 */
START_TEST(searchQuery_paired2)
{
	/* Both first and second query files are NULL. */
	{
		char *qryFile1 = NULL;
		char *qryFile2 = NULL;
		int numQrs = 100;
		int qryLen = 100;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLen = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries_paired2(qryFile1, qryFile2, numQrs, qryLen, refFile,
					seedLen, matchFile, maxNumHits, &numMatches);
		remove(refFile);
	}

	/* Number of queries is 0. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 0;
		int qryLen = 100;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFile, "w");
		fputs(">r1\nAATTCCCAACGTACGTAAGTATTACCGGAAGGAATC", refFilePtr);
		fclose(refFilePtr);
		int seedLen = 8;
		int numMatches = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		searchQueries_paired2(qryFile1, qryFile2, numQrs, qryLen, refFile,
					seedLen, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
	}

	/* Query length is 0. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 10;
		int qryLen = 0;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLen = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries_paired2(qryFile1, qryFile2, numQrs, qryLen, refFile,
					seedLen, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
	}

	/* Reference file name is NULL. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 10;
		int qryLen = 100;
		char *refFile = NULL;
		int seedLen = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries_paired2(qryFile1, qryFile2, numQrs, qryLen, refFile,
				seedLen, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
	}

	/* Seed length is 0. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 10;
		int qryLen = 100;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLen = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries_paired2(qryFile1, qryFile2, numQrs, qryLen, refFile,
				seedLen, matchFile, maxNumHits, &numMatches);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
	}

	/* Match file is NULL. */
	{
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		int numQrs = 10;
		int qryLen = 100;
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLen = 8;
		char *matchFile = NULL;
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries_paired2(qryFile1, qryFile2, numQrs, qryLen, refFile,
				seedLen, matchFile, maxNumHits, &numMatches);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
	}

	/* 1 query pair and 1 reference. */
	{
		preprocessCreate();
		char qryFile1[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile1, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry1.fa");
		FILE *qryFilePtr1 = fopen(qryFile1, "w");
		fputs(">q1_1\nACGTACGTAA", qryFilePtr1);
		fclose(qryFilePtr1);

		char qryFile2[MAX_FILE_NAME_LENGTH];
		sprintf(qryFile2, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry2.fa");
		FILE *qryFilePtr2 = fopen(qryFile2, "w");
		fputs(">q1_2\nACCGGAAGGA", qryFilePtr2);
		fclose(qryFilePtr2);

		int numQrs = 1;
		int qryLen = 10;

		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFile, "w");
		/* AATTCCCA ACGTACGT AAGTATTA CCGGAAGG AATC */
		fputs(">r1\nAATTCCCAACGTACGTAAGTATTACCGGAAGGAATC", refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFile);

		int seedLen = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		searchQueries_paired2(qryFile1, qryFile2, numQrs, qryLen, refFile,
					seedLen, matchFile, maxNumHits, &numMatches);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "1\tq1_1\t0\t1\t0\t8\t8\t0\t6\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when there is 1 query pair and 1 "
					"reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine3[] = "1\tq1_1\t1\t1\t0\t6\t8\t0\t6\n";
		if (strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when there is 1 query pair and 1 "
					"reference (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine2[] = "1\tq1_2\t0\t2\t0\t23\t24\t0\t6\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when there is 1 query pair and 1 "
					"reference (case 3).\n");

		if (numMatches != 3)
			fail("Incorrect behavior when there is 1 query pair and 1 "
					"reference (case 4).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFile1);
		remove(qryFile2);
		remove(refFile);
		refMapFree();
		preprocessDelete();
	}
}
END_TEST


/**
 * Tests @a searchQueries_gpu function.
 */
START_TEST(searchQueries_gpu)
{
	/* Query file has only 1 query and reference file has only 1 reference
	 * and max number of hits is 1. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nACGTACGTGG", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 1;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		fputs(">r1\nACGTACGTAATTCCCA", refFilePtr);
		fclose(refFilePtr);
		int numReferences = 1;
		int refLength = 16;
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 1;
		searchQueries_gpu(qryFileName, numQueries, qryLength, refFileName,
				numReferences, refLength, seedLength, matchFile, maxNumHits);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t0\t0\t0\t4\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine2[] = "q1\t1\t0\t-2\t0\t0\t4\t0\t4\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 2).\n");

		if (fgets(line, MAX_LINE_LENGTH, matchFilePtr) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 3).\n");
		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query file has only 1 query and reference file has only 1 reference
	 * and max number of hits is 3. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nACGTACGTGG", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 1;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		fputs(">r1\nACGTACGTAATTCCCA", refFilePtr);
		fclose(refFilePtr);
		int numReferences = 1;
		int refLength = 16;
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 3;
		searchQueries_gpu(qryFileName, numQueries, qryLength, refFileName,
				numReferences, refLength, seedLength, matchFile, maxNumHits);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t0\t0\t0\t4\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference and max number of "
					"hits is 3 (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine4[] = "q1\t1\t0\t-2\t0\t0\t4\t0\t4\n";
		if (strcmp(line, tmpLine4) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference and max number "
					"of hits is 3 (case 2).\n");

		if (fgets(line, MAX_LINE_LENGTH, matchFilePtr) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference and max number "
					"of hits is 3 (case 3).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query file has 2 queries and reference file has only 1 reference
	 * and max number of hits is 1. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nACGTACGTGG\n>q2\nACAATTCCCA", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 2;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		fputs(">r1\nACGTACGTAATTCCCA", refFilePtr);
		fclose(refFilePtr);
		int numReferences = 1;
		int refLength = 16;
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 1;
		searchQueries_gpu(qryFileName, numQueries, qryLength, refFileName,
				numReferences, refLength, seedLength, matchFile, maxNumHits);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t0\t0\t0\t4\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has only 1 reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine2[] = "q1\t1\t0\t-2\t0\t0\t4\t0\t4\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has only 1 reference (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine3[] = "q2\t0\t0\t6\t0\t15\t19\t0\t4\n";
		if (strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has only 1 reference (case 3).\n");

		if (fgets(line, MAX_LINE_LENGTH, matchFilePtr) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has only 1 reference (case 4).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query file has 2 queries and reference file has only 1 reference
	 * and max number of hits is 3. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nACGTACGTGG\n>q2\nACAATTCCCA", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 2;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		fputs(">r1\nACGTACGTAATTCCCA", refFilePtr);
		fclose(refFilePtr);

		int numReferences = 1;
		int refLength = 16;
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 3;
		searchQueries_gpu(qryFileName, numQueries, qryLength, refFileName,
				numReferences, refLength, seedLength, matchFile, maxNumHits);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t0\t0\t0\t4\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has only 1 reference and max number of "
					"hits is 3 (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine4[] = "q1\t1\t0\t-2\t0\t0\t4\t0\t4\n";
		if (strcmp(line, tmpLine4) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has only 1 reference and max number of "
					"hits is 3 (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine7[] = "q2\t0\t0\t6\t0\t15\t19\t0\t4\n";
		if (strcmp(line, tmpLine7) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has only 1 reference and max number of "
					"hits is 3 (case 3).\n");

		if (fgets(line, MAX_LINE_LENGTH, matchFilePtr) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has only 1 reference and max number of "
					"hits is 3 (case 4).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query file has 2 queries and reference file has 2 references
	 * and max number of hits is 1. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nACGTACGTGG\n>q2\nACAATTCCCA", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 2;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		fputs(">r1\nACGTACGTAATTCCCA\n>r2\nAATTCCCAACGTACGT", refFilePtr);
		fclose(refFilePtr);
		int numReferences = 2;
		int refLength = 16;
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 1;
		searchQueries_gpu(qryFileName, numQueries, qryLength, refFileName,
				numReferences, refLength, seedLength, matchFile, maxNumHits);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t0\t0\t0\t4\t0\t4\n";
		char tmpLine2[] = "q1\t0\t1\t8\t0\t0\t4\t21\t25\n";
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has 2 references (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine3[] = "q1\t1\t0\t-2\t0\t0\t4\t0\t4\n";
		char tmpLine4[] = "q1\t1\t1\t6\t0\t0\t4\t21\t25\n";
		if (strcmp(line, tmpLine3) != 0 && strcmp(line, tmpLine4) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has 2 references (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine5[] = "q2\t0\t0\t6\t0\t15\t19\t0\t4\n";
		char tmpLine6[] = "q2\t0\t1\t-2\t0\t15\t19\t21\t25\n";
		if (strcmp(line, tmpLine5) != 0 && strcmp(line, tmpLine6) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has 2 references (case 3).\n");

		if (fgets(line, MAX_LINE_LENGTH, matchFilePtr) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has 2 references (case 4).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query file has 2 queries and reference file has 2 references
	 * and max number of hits is 3. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nACGTACGTGG\n>q2\nACAATTCCCA", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 2;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		fputs(">r1\nACGTACGTAATTCCCA\n>r2\nAATTCCCAACGTACGT", refFilePtr);
		fclose(refFilePtr);
		int numReferences = 2;
		int refLength = 16;
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 3;
		searchQueries_gpu(qryFileName, numQueries, qryLength, refFileName,
				numReferences, refLength, seedLength, matchFile, maxNumHits);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t0\t0\t0\t4\t0\t4\n";
		char tmpLine2[] = "q1\t0\t1\t8\t0\t0\t4\t21\t25\n";
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has 2 references and max number of "
					"hits is 3 (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has 2 references and max number of "
					"hits is 3 (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine4[] = "q1\t1\t0\t-2\t0\t0\t4\t0\t4\n";
		char tmpLine5[] = "q1\t1\t1\t6\t0\t0\t4\t21\t25\n";
		if (strcmp(line, tmpLine4) != 0 && strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has 2 references and max number of "
					"hits is 3 (case 3).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine4) != 0 && strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has 2 references and max number of "
					"hits is 3 (case 4).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine7[] = "q2\t0\t0\t6\t0\t15\t19\t0\t4\n";
		char tmpLine8[] = "q2\t0\t1\t-2\t0\t15\t19\t21\t25\n";
		if (strcmp(line, tmpLine7) != 0 && strcmp(line, tmpLine8) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has 2 references and max number of "
					"hits is 3 (case 5).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine7) != 0 && strcmp(line, tmpLine8) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has 2 references and max number of "
					"hits is 3 (case 6).\n");

		if (fgets(line, MAX_LINE_LENGTH, matchFilePtr) != 0)
			fail("Incorrect behavior when query file has 2 queries "
					"and reference file has 2 references and max number of "
					"hits is 3 (case 7).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query file has 1 query and reference file has 2 references
	 * and max number of hits is 3 and the query has more than max number
	 * of hits. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nAAAAAA\n", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 1;
		int qryLength = 6;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		fputs(">r1\nAAAAACCCCC\n>r2\nGGGGGAAAAA\n", refFilePtr);
		fclose(refFilePtr);
		int numReferences = 2;
		int refLength = 10;
		int seedLength = 5;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 2;
		searchQueries_gpu(qryFileName, numQueries, qryLength, refFileName,
				numReferences, refLength, seedLength, matchFile, maxNumHits);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		char tmpLine1[] = "q1\t0\t0\t0\t0\t0\t4\t0\t4\n";
		char tmpLine2[] = "q1\t0\t1\t5\t0\t0\t4\t15\t19\n";
		char tmpLine3[] = "q1\t0\t0\t-1\t0\t0\t4\t0\t4\n";
		char tmpLine4[] = "q1\t0\t1\t4\t0\t0\t4\t15\t19\n";

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0
				&& strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0
				&& strcmp(line, tmpLine4) != 0)
			fail("Incorrect behavior when query file has 1 query "
					"and reference file has 2 references and max number of "
					"hits is 2 (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0
				&& strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0
				&& strcmp(line, tmpLine4) != 0)
			fail("Incorrect behavior when query file has 1 query "
					"and reference file has 2 references and max number of "
					"hits is 2 (case 2).\n");

		if (fgets(line, MAX_LINE_LENGTH, matchFilePtr) != 0)
			fail("Incorrect behavior when query file has 1 query "
					"and reference file has 2 references and max number of "
					"hits is 2 (case 3).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}
}
END_TEST


/**
 * Tests the @a getNumIters function.
 */
START_TEST(getNumIters)
{
	/* There is memory available on the GPU global memory to store query
	 * sequences and results; maximum hits per query = 1. */
	{
		int numQrs = 24000000;
		long gpuTotalMem = 900000000;
		int gpuMemMargin = 100000000;
		int refLen = 250000000;
		int seedLen = 12;
		int qryLen = 100;
		int maxNumHits = 1;
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, seedLen);
		int numTuples = refLen / seedLen;
		int keysTotalMem = numKeys * sizeof(int);
		int valsTotalMem = numTuples * sizeof(int);
		int qrySeqMem = qryLen * sizeof(char);
		int refIdxMem = sizeof(char);
		int shiftMem = sizeof(int);
		int posMem = sizeof(int);
		int numIter;
		long numQrsPerIter, numQrsLastIter;

		getNumIters_wrap(&numIter, &numQrsPerIter, &numQrsLastIter, numQrs,
					maxNumHits, gpuTotalMem, gpuMemMargin, keysTotalMem,
					valsTotalMem, qrySeqMem, refIdxMem, shiftMem, posMem);
		if (numIter != 367 || numQrsPerIter != 65534 || numQrsLastIter != 14556)
			fail("Incorrect behavior when there is memory available on the "
					"GPU global memory to store query sequences and results "
					"and max hits per query is 1.\n");
	}

	/* There is memory available on the GPU global memory to store query
	 * sequences and results; maximum hits per query = 10. */
	{
		int numQrs = 24000000;
		long gpuTotalMem = 900000000;
		int gpuMemMargin = 100000000;
		int refLen = 250000000;
		int seedLen = 12;
		int qryLen = 100;
		int maxNumHits = 10;
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, seedLen);
		int numTuples = refLen / seedLen;
		int keysTotalMem = numKeys * sizeof(int);
		int valsTotalMem = numTuples * sizeof(int);
		int qrySeqMem = qryLen * sizeof(char);
		int refIdxMem = sizeof(char);
		int shiftMem = sizeof(int);
		int posMem = sizeof(int);
		int numIter;
		long numQrsPerIter, numQrsLastIter;

		getNumIters_wrap(&numIter, &numQrsPerIter, &numQrsLastIter, numQrs,
					maxNumHits, gpuTotalMem, gpuMemMargin, keysTotalMem,
					valsTotalMem, qrySeqMem, refIdxMem, shiftMem, posMem);
		if (numIter != 367 || numQrsPerIter != 65534 || numQrsLastIter != 14556)
			fail("Incorrect behavior when there is memory available on the "
					"GPU global memory to store query sequences and results "
					"and max hits per query = 10.\n");
	}

	/* There is no memory available on the GPU global memory to store
	 * query sequences and results; maximum hits per query = 1. */
	{
		int numQrs = 24000000;
		long gpuTotalMem = 900000000;
		int gpuMemMargin = 100000000;
		uint refLen = 4000000000;
		int seedLen = 12;
		int qryLen = 100;
		int maxNumHits = 1;
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, seedLen);
		int numTuples = refLen / seedLen;
		int keysTotalMem = numKeys * sizeof(int);
		int valsTotalMem = numTuples * sizeof(int);
		int qrySeqMem = qryLen * sizeof(char);
		int refIdxMem = sizeof(char);
		int shiftMem = sizeof(int);
		int posMem = sizeof(int);
		int numIter;
		long numQrsPerIter, numQrsLastIter;

		getNumIters_wrap(&numIter, &numQrsPerIter, &numQrsLastIter, numQrs,
					maxNumHits, gpuTotalMem, gpuMemMargin, keysTotalMem,
					valsTotalMem, qrySeqMem, refIdxMem, shiftMem, posMem);
		if (numIter != -1 || numQrsPerIter != -1 || numQrsLastIter != -1)
			fail("Incorrect behavior when there is no memory available on the "
					"GPU global memory to store query sequences and results and "
					"max number of hits per query is 1."
					"\n");
	}

	/* There is no memory available on the GPU global memory to store
	 * query sequences and results; maximum hits per query = 10. */
	{
		int numQrs = 24000000;
		long gpuTotalMem = 900000000;
		int gpuMemMargin = 100000000;
		uint refLen = 4000000000;
		int seedLen = 12;
		int qryLen = 100;
		int maxNumHits = 10;
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, seedLen);
		int numTuples = refLen / seedLen;
		int keysTotalMem = numKeys * sizeof(int);
		int valsTotalMem = numTuples * sizeof(int);
		int qrySeqMem = qryLen * sizeof(char);
		int refIdxMem = sizeof(char);
		int shiftMem = sizeof(int);
		int posMem = sizeof(int);
		int numIter;
		long numQrsPerIter, numQrsLastIter;

		getNumIters_wrap(&numIter, &numQrsPerIter, &numQrsLastIter, numQrs,
					maxNumHits, gpuTotalMem, gpuMemMargin, keysTotalMem,
					valsTotalMem, qrySeqMem, refIdxMem, shiftMem, posMem);
		if (numIter != -1 || numQrsPerIter != -1 || numQrsLastIter != -1)
			fail("Incorrect behavior when there is no memory available on the "
					"GPU global memory to store query sequences and results and "
					"max number of hits per query is 10.\n");
	}

	/* There is memory available on the GPU global memory to store query
	 * sequences and results; maximum hits per query = 1; number of queries
	 * is small; reference length is small. */
	{
		int numQrs = 2;
		long gpuTotalMem = 900000000;
		int gpuMemMargin = 100000000;
		int refLen = 16;
		int seedLen = 8;
		int qryLen = 10;
		int maxNumHits = 1;
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, seedLen);
		int numTuples = refLen / seedLen;
		int keysTotalMem = numKeys * sizeof(int);
		int valsTotalMem = numTuples * sizeof(int);
		int qrySeqMem = qryLen * sizeof(char);
		int refIdxMem = sizeof(char);
		int shiftMem = sizeof(int);
		int posMem = sizeof(int);
		int numIter;
		long numQrsPerIter, numQrsLastIter;

		getNumIters_wrap(&numIter, &numQrsPerIter, &numQrsLastIter, numQrs,
					maxNumHits, gpuTotalMem, gpuMemMargin, keysTotalMem,
					valsTotalMem, qrySeqMem, refIdxMem, shiftMem, posMem);
		if (numIter != 1 || numQrsPerIter != 2 || numQrsLastIter != 2)
			fail("Incorrect behavior when there is memory available on the "
					"GPU global memory to store query sequences and results "
					"and max hits per query is 1 and number of queries is "
					"small and reference length is small.\n");
	}

	/* There is memory available on the GPU global memory to store query
	 * sequences and results; maximum hits per query = 10; number of queries
	 * is small; reference length is small. */
	{
		int numQrs = 2;
		long gpuTotalMem = 900000000;
		int gpuMemMargin = 100000000;
		int refLen = 16;
		int seedLen = 8;
		int qryLen = 10;
		int maxNumHits = 10;
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, seedLen);
		int numTuples = refLen / seedLen;
		int keysTotalMem = numKeys * sizeof(int);
		int valsTotalMem = numTuples * sizeof(int);
		int qrySeqMem = qryLen * sizeof(char);
		int refIdxMem = sizeof(char);
		int shiftMem = sizeof(int);
		int posMem = sizeof(int);
		int numIter;
		long numQrsPerIter, numQrsLastIter;

		getNumIters_wrap(&numIter, &numQrsPerIter, &numQrsLastIter, numQrs,
					maxNumHits, gpuTotalMem, gpuMemMargin, keysTotalMem,
					valsTotalMem, qrySeqMem, refIdxMem, shiftMem, posMem);
		if (numIter != 1 || numQrsPerIter != 2 || numQrsLastIter != 2)
			fail("Incorrect behavior when there is memory available on the "
					"GPU global memory to store query sequences and results "
					"and max hits per query is 10 and number of queries is "
					"small and reference length is small.\n");
	}
}
END_TEST


/**
 * Tests @a getNumBlocks function.
 */
START_TEST(getNumBlocks)
{
	/* There are 65535 x 65535 blocks and number of queries per iteration
	 * is 100. */
	{
		int blocksX, blocksY;
		long availBlocksX = 65535;
		long availBlocksY = 65535;
		long numQrsPerIter = 100;

		getNumBlocks_wrap(&blocksX, &blocksY, availBlocksX, availBlocksY,
				numQrsPerIter);

		if (blocksX != 10 && blocksY != 10)
			fail("Incorrect behavior when there are 65535 x 65535 blocks "
					"and number of queries per iteration is 100.\n");
	}

	/* There are 100 x 100 blocks and number of queries per iteration
	 * is 100. */
	{
		int blocksX, blocksY;
		long availBlocksX = 100;
		long availBlocksY = 100;
		long numQrsPerIter = 100;

		getNumBlocks_wrap(&blocksX, &blocksY, availBlocksX, availBlocksY,
				numQrsPerIter);

		if (blocksX != 10 && blocksY != 10)
			fail("Incorrect behavior when there are 100 x 100 blocks "
					"and number of queries per iteration is 100.\n");
	}

	/* There are 50 x 50 blocks and number of queries per iteration
	 * is 3000. */
	{
		int blocksX, blocksY;
		long availBlocksX = 50;
		long availBlocksY = 50;
		long numQrsPerIter = 3000;

		getNumBlocks_wrap(&blocksX, &blocksY, availBlocksX, availBlocksY,
				numQrsPerIter);

		if (blocksX != 50 && blocksY != 50)
			fail("Incorrect behavior when there are 50 x 50 blocks "
					"and number of queries per iteration is 3000.\n");
	}
}
END_TEST


/**
 * Tests @a searchQueries3 function.
 */
START_TEST(searchQuery3)
{
	/* Query file name is NULL. */
	{
		char *qryFileName = NULL;
		int numQueries = 100;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		int numThreads = 1;
		searchQueries3(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches, numThreads);
	}

	/* Number of queries is 0. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 0;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		int numThreads = 1;
		searchQueries3(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches, numThreads);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query length is 0. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 0;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		int numThreads = 1;
		searchQueries3(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches, numThreads);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Reference file name is NULL. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 100;
		char *refFileName = NULL;
		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		int numThreads = 1;
		searchQueries3(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches, numThreads);
		remove(matchFile);
		remove(qryFileName);
	}

	/* Seed length is 0. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		int numThreads = 1;
		searchQueries3(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches, numThreads);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Match file is NULL. */
	{
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		int numQueries = 10;
		int qryLength = 100;
		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		int seedLength = 8;
		char *matchFile = NULL;
		uint maxNumHits = 10;
		int numMatches = 0;
		int numThreads = 1;
		searchQueries3(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches, numThreads);
		remove(qryFileName);
		remove(refFileName);
	}

	/* Query file has only 1 query and reference file has only 1 reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nGGACGTACGT", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 1;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		fputs(">r1\nAATTCCCAACGTACGT", refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		int numMatches = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numThreads = 1;
		searchQueries3(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches, numThreads);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine2[] = "q1\t1\t0\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 2).\n");

		if (numMatches != 2)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has only 1 reference (case 3).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}

	/* Query file has only 1 query and reference file has more than
	 * 1 reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		fputs(">q1\nGGACGTACGT", qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 1;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		char refSeq[] = ">r1\nAATTCCCAACGTACGT\n>r2\nGACGTACGACGTACGT"
				"\n>r3\nAAAAAAAAAAAAAAAA";
		fputs(refSeq, refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		int numMatches = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numThreads = 1;
		searchQueries3(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches, numThreads);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		char tmpLine2[] = "q1\t0\t1\t-1\t0\t0\t4\n";
		char tmpLine3[] = "q1\t0\t1\t6\t8\t0\t4\n";
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 3).\n");

		char tmpLine4[] = "q1\t1\t0\t8\t8\t0\t4\n";
		char tmpLine5[] = "q1\t1\t1\t8\t8\t0\t4\n";
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine4) != 0 && strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 4).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine4) != 0 && strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 5).\n");

		if (numMatches != 5)
			fail("Incorrect behavior when query file has only 1 query "
					"and reference file has more than 1 reference "
					"(case 6).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}

	/* Query file has more than 1 query and reference file has only 1
	 * reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		char qrySeq[] = ">q1\nGGACGTACGT\n>q2\nAATTCCCACC\n>q3\nGAATTCCCAG\n";
		fputs(qrySeq, qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 3;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		char refSeq[] = ">r1\nAATTCCCAACGTACGT\n";
		fputs(refSeq, refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numMatches = 0;
		int numThreads = 1;
		searchQueries3(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches, numThreads);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		if (strcmp(line, tmpLine1) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine2[] = "q1\t1\t0\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine2) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine3[] = "q2\t0\t0\t0\t0\t15\t19\n";
		if (strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 3).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine4[] = "q3\t0\t0\t-1\t0\t30\t34\n";
		if (strcmp(line, tmpLine4) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 4).\n");

		printf("numMatches = %d\n", numMatches);
		if (numMatches != 4)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has only 1 reference (case 5).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(qryFileName);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}

	/* Query file has more than 1 query and reference file has more than
	 * 1 reference. */
	{
		preprocessCreate();
		char qryFileName[MAX_FILE_NAME_LENGTH];
		sprintf(qryFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_qry.fa");
		FILE *qryFilePtr = fopen(qryFileName, "w");
		char qrySeq[] = ">q1\nGGACGTACGT\n>q2\nAATTCCCACC\n>q3\nGAATTCCCAG\n";
		fputs(qrySeq, qryFilePtr);
		fclose(qryFilePtr);
		int numQueries = 3;
		int qryLength = 10;

		char refFileName[MAX_FILE_NAME_LENGTH];
		sprintf(refFileName, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *refFilePtr = fopen(refFileName, "w");
		char refSeq[] = ">r1\nAATTCCCAACGTACGT\n>r2\nGACGTACGACGTACGT"
				"\n>r3\nAAAAAAAAAAAAAAAA";
		fputs(refSeq, refFilePtr);
		fclose(refFilePtr);

		refMapCreate(refFileName);

		int seedLength = 8;
		int numMatches = 0;
		char matchFile[MAX_FILE_NAME_LENGTH];
		sprintf(matchFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
				"test_matchFile.fa");
		uint maxNumHits = 10;
		int numThreads = 1;
		searchQueries3(qryFileName, numQueries, qryLength, refFileName,
				seedLength, matchFile, maxNumHits, &numMatches, numThreads);
		FILE *matchFilePtr = fopen(matchFile, "r");
		char line[MAX_LINE_LENGTH];
		char tmpLine1[] = "q1\t0\t0\t6\t8\t0\t4\n";
		char tmpLine2[] = "q1\t0\t1\t-1\t0\t0\t4\n";
		char tmpLine3[] = "q1\t0\t1\t6\t8\t0\t4\n";
		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 1).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 2).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine1) != 0 && strcmp(line, tmpLine2) != 0
				&& strcmp(line, tmpLine3) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 3).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine4[] = "q1\t1\t0\t8\t8\t0\t4\n";
		char tmpLine5[] = "q1\t1\t1\t8\t8\t0\t4\n";
		if (strcmp(line, tmpLine4) != 0 && strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 4).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		if (strcmp(line, tmpLine4) != 0 && strcmp(line, tmpLine5) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 5).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine6[] = "q2\t0\t0\t0\t0\t15\t19\n";
		if (strcmp(line, tmpLine6) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 6).\n");

		fgets(line, MAX_LINE_LENGTH, matchFilePtr);
		char tmpLine7[] = "q3\t0\t0\t-1\t0\t30\t34\n";
		if (strcmp(line, tmpLine7) != 0)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 7).\n");

		if (numMatches != 7)
			fail("Incorrect behavior when query file has more than 1 query "
					"and reference file has more than 1 reference (case 8).\n");

		fclose(matchFilePtr);
		remove(matchFile);
		remove(refFileName);
		refMapFree();
		preprocessDelete();
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *searchSuite(void)
{
	Suite *s = suite_create("search");

	/* Core test case. */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, searchQuery);
	tcase_add_test(testCaseCore, searchQuery2);
	tcase_add_test(testCaseCore, searchQuery3);
//	tcase_add_test(testCaseCore, searchQueries_gpu);
//	tcase_add_test(testCaseCore, searchQuery_paired);
	tcase_add_test(testCaseCore, searchQuery_paired2);
	tcase_add_test(testCaseCore, getNumIters);
	tcase_add_test(testCaseCore, getNumBlocks);
	suite_add_tcase (s, testCaseCore);

	return s;
}
