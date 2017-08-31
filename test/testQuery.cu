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
#include "testQuery.h"
#include "../src/query.h"
#include "../src/common.h"


/**
 * Tests qryCreate function.
 */
START_TEST(create)
{
	/* Test behavior when query name is NULL. */
	{
		char *name = NULL;
		char *seq = "ACGTACGTAA";
		char *revComp = "TTACGTACGT";
		ulong nameOffset = 0;
		ulong seqOffset = 5;
		Query *query = qryCreate(name, seq, revComp, nameOffset, seqOffset);
		if (query != NULL)
			fail("Incorrect behavior when query name is NULL.\n");
	}

	/* Test behavior when query sequence is NULL. */
	{
		char *name = "qry1";
		char *seq = NULL;
		char *revComp = "TTACGTACGT";
		ulong nameOffset = 0;
		ulong seqOffset = 5;
		Query *query = qryCreate(name, seq, revComp, nameOffset, seqOffset);
		if (query != NULL)
			fail("Incorrect behavior when query sequence is NULL.\n");
	}

	/* Test behavior when reverse complement sequence is NULL. */
	{
		char *name = "qry1";
		char *seq = "ACGTACGTAA";
		char *revComp = NULL;
		ulong nameOffset = 0;
		ulong seqOffset = 5;
		Query *query = qryCreate(name, seq, revComp, nameOffset, seqOffset);
		if (query != NULL)
			fail("Incorrect behavior when reverse complement of query "
					"sequence is NULL.\n");
	}

	/* Test behavior when query name offset is greater than
	 * to query sequence offset. */
	{
		char *name = "qry1";
		char *seq = "ACGTACGTAA";
		char *revComp = "TTACGTACGT";
		ulong nameOffset = 5;
		ulong seqOffset = 0;
		Query *query = qryCreate(name, seq, revComp, nameOffset, seqOffset);
		if (query != NULL)
			fail("Incorrect behavior when query name offset is greater "
					"than query sequence offset.\n");
	}

	/* Test behavior when query name offset is equal to query
	 * sequence offset. */
	{
		char *name = "qry1";
		char *seq = "ACGTACGTAA";
		char *revComp = "TTACGTACGT";
		ulong nameOffset = 0;
		ulong seqOffset = 0;
		Query *query = qryCreate(name, seq, revComp, nameOffset, seqOffset);
		if (query != NULL)
			fail("Incorrect behavior when query name offset is equal "
					"to query sequence offset.\n");
	}

	/* Test behavior when parameter values are valid. */
	{
		char *name = "qry1";
		char *seq = "ACGTACGTAA";
		char *revComp = "TTACGTACGT";
		ulong nameOffset = 0;
		ulong seqOffset = 5;
		Query *query = qryCreate(name, seq, revComp, nameOffset, seqOffset);
		if (query == NULL || strcmp(query->name, name) != 0
				|| strcmp(query->seq, seq) != 0
				|| strcmp(query->revCompSeq, revComp) != 0
				|| query->nameOffset != nameOffset
				|| query->seqOffset != seqOffset)
			fail("Incorrect behavior when parameter values are valid.\n");
		qryDelete(query);
	}
}
END_TEST


/**
 * Tests @ qryGetNext function.
 */
START_TEST(getNext)
{
	char file[MAX_FILE_NAME_LENGTH];
	sprintf(file, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_getNext.fa");
	FILE *filePtr = fopen(file, "w+");
	fputs(">q1\n", filePtr);
	fputs("ACGTACGTAA\n", filePtr);
	fputs(">q2\n", filePtr);
	fputs("ACGGACGGCC", filePtr);
	rewind(filePtr);

	Query *qry = qryCreate2();
	qry = qryGetNext2(filePtr, qry);
	if (strcmp(qry->name, "q1") != 0 || strcmp(qry->seq, "ACGTACGTAA") != 0
			|| strcmp(qry->revCompSeq, "TTACGTACGT") != 0
			|| qry->nameOffset != 0 || qry->seqOffset != 4)
		fail("Incorrect behavior when fetching the first query.\n");

	qry = qryGetNext2(filePtr, qry);
	if (strcmp(qry->name, "q2") != 0 || strcmp(qry->seq, "ACGGACGGCC") != 0
			|| strcmp(qry->revCompSeq, "GGCCGTCCGT") != 0
			|| qry->nameOffset != 15 || qry->seqOffset != 19)
		fail("Incorrect behavior when fetching the second query.\n");
	fclose(filePtr);
	remove(file);
}
END_TEST


/**
 * Creates test suite.
 */
Suite *querySuite(void)
{
	Suite *s = suite_create("query");

	/* Core test case. */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, create);
	tcase_add_test(testCaseCore, getNext);
	suite_add_tcase(s, testCaseCore);

	return s;
}
