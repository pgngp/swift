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
#include "testLookupTable.h"
#include "../src/lookupTable.h"
#include "../src/lookupTable2.h"
#include "../src/lookupTable3.h"
#include "../src/lookupTable4.h"
#include "../src/refMap.h"
#include "../src/common.h"
#include "../src/preprocess.h"
#include <limits.h>


/**
 * Tests @a lookupTableMapQry function.
 */
START_TEST(lookupTableMapQry)
{
	/* Test for correct results with a normal query (case 1). */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp1.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);
		int seedLength = 8;
		uint maxNumHits = 10;
		lookupTableCreate(refFile, seedLength, maxNumHits);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		HitList *result = lookupTableMapQry(query, seqLength, isRevComplement);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 6, 8);
		if (result != NULL
				&& hitListCmpNodes(result, node1) == 0
				&& hitListCmpNodes(result->next, node2) == 0)
		{}
		else
			fail("Incorrect behavior when parameters have normal values "
					"(case 1).\n");
		qryDelete(query);
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);
		lookupTableDelete();
		remove(refFile);
		preprocessDelete();
	}

	/* Test for correct results with a normal query (case 2). */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp2.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAACGTACGTATTTACGTA\n"
				">ref3\nCGTACGTAAAATAACGTCGTACGTAAACGTACGTACGTAA", filePtr);
		fclose(filePtr);
		int seedLength = 8;
		uint maxNumHits = 10;
		lookupTableCreate(refFile, seedLength, maxNumHits);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		HitList *result = lookupTableMapQry(query, seqLength, isRevComplement);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(1, 8, 8);
		HitList *node3 = hitListCreateNode(1, 23, 24);
		HitList *node4 = hitListCreateNode(2, -1, 0);
		HitList *node5 = hitListCreateNode(2, 30, 32);
		if (result != NULL && hitListCmpNodes(result, node1) == 0
				&& hitListCmpNodes(result->next, node2) == 0
				&& hitListCmpNodes(result->next->next, node3) == 0
				&& hitListCmpNodes(result->next->next->next, node4) == 0
				&& hitListCmpNodes(result->next->next->next->next, node5) == 0)
		{ }
		else
			fail("Incorrect behavior when parameters have normal values "
					"(case 2).\n");
		qryDelete(query);
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		hitListDelete(node4);
		hitListDelete(node5);
		lookupTableDelete();
		remove(refFile);
		preprocessDelete();
	}

	/* Test for correct results with the reverse complement of a query. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp2.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAACGTACGTATTTACGTA\n"
				">ref3\nCGTACGTAAAATAACGTCGTACGTAAACGTACGTACGTAA", filePtr);
		fclose(filePtr);
		int seedLength = 8;
		uint maxNumHits = 10;
		lookupTableCreate(refFile, seedLength, maxNumHits);
		Query *query = qryCreate("qry1", "TTACGTACGT", "ACGTACGTAA", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 1;
		HitList *result = lookupTableMapQry(query, seqLength, isRevComplement);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(1, 8, 8);
		HitList *node3 = hitListCreateNode(1, 23, 24);
		HitList *node4 = hitListCreateNode(2, -1, 0);
		HitList *node5 = hitListCreateNode(2, 30, 32);
		if (result != NULL && hitListCmpNodes(result, node1) == 0
				&& hitListCmpNodes(result->next, node2) == 0
				&& hitListCmpNodes(result->next->next, node3) == 0
				&& hitListCmpNodes(result->next->next->next, node4) == 0
				&& hitListCmpNodes(result->next->next->next->next, node5) == 0)
		{ }
		else
			fail("Incorrect behavior when sequence is a reverse complement "
					"of the query.\n");
		qryDelete(query);
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		hitListDelete(node4);
		hitListDelete(node5);
		lookupTableDelete();
		remove(refFile);
		preprocessDelete();
	}
}
END_TEST


/**
 * Tests @a lookupTable2Create function.
 */
START_TEST(lookupTable2Create)
{
	/* 1 reference. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		uint *refPos;
		int seedLen = 8, maxHitsPerQry = 2, tupleIgnoreThreshold = 15;

		lookupTable2Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 65536 && lookupTable != NULL
					&& lookupTable[3662] == 1 && lookupTable[58596] == 2
					&& refPos != NULL && refPos[0] == UINT_MAX && refPos[1]
					== 14 && refPos[2] == 6 && numRepeatsPerTuple[3662] == 1
					&& numRepeatsPerTuple[58596] == 1)
		{}
		else
			fail("Incorrect behavior when there is 1 reference.\n");
		lookupTable2Delete();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references, whole sequence of a given reference is in the same line. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAACGTACGTATTTACG\n"
				">ref3\nCGTACGTAAAATAACGTCGTACGTAAACGTACGTACGTAA", filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		uint *refPos;
		int seedLen = 8, maxHitsPerQry = 2, tupleIgnoreThreshold = 15;

		lookupTable2Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 65536 && lookupTable != NULL
				&& lookupTable[3648] == 1
				&& lookupTable[3662] == 2
				&& lookupTable[14649] == 3
				&& lookupTable[17636] == 5
				&& lookupTable[20032] == 6
				&& lookupTable[37056] == 7
				&& lookupTable[37782] == 8
				&& lookupTable[58416] == 9
				&& lookupTable[58596] == 10
				&& lookupTable[58597] == 12
				&& lookupTable[58599] == 13
				&& lookupTable[64740] == 14
				&& refPos != NULL
				&& refPos[0] == UINT_MAX
				&& refPos[1] == 69
				&& refPos[2] == 130
				&& refPos[3] == 98
				&& refPos[4] == 77
				&& refPos[5] == 30
				&& refPos[6] == 122
				&& refPos[7] == 106
				&& refPos[8] == 53
				&& refPos[9] == 14
				&& refPos[10] == 61
				&& refPos[11] == 6
				&& refPos[12] == 38
				&& refPos[13] == 114
				&& refPos[14] == 22
				&& numRepeatsPerTuple[3648] == 1
				&& numRepeatsPerTuple[3662] == 1
				&& numRepeatsPerTuple[14649] == 2
				&& numRepeatsPerTuple[17636] == 1
				&& numRepeatsPerTuple[20032] == 1
				&& numRepeatsPerTuple[37056] == 1
				&& numRepeatsPerTuple[37782] == 1
				&& numRepeatsPerTuple[58416] == 1
				&& numRepeatsPerTuple[58596] == 2
				&& numRepeatsPerTuple[58597] == 1
				&& numRepeatsPerTuple[58599] == 1
				&& numRepeatsPerTuple[64740] == 1
		)
		{}
		else
			fail("Incorrect behavior when there are 3 references, whole "
					"sequence of a given reference is in the same line.\n");
		lookupTable2Delete();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references, each sequence split into different lines. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTT\nACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAA\nCGTACGTATTTACG\n"
				">ref3\nCGTACGTAAA\nATAACGTCGT\nACGTAAACGT\nACGTACGTAA",
				filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		uint *refPos;
		int seedLen = 8, maxHitsPerQry = 2, tupleIgnoreThreshold = 4;

		lookupTable2Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 65536 && lookupTable != NULL
					&& lookupTable[3648] == 1
					&& lookupTable[3662] == 2
					&& lookupTable[14649] == 3
					&& lookupTable[17636] == 5
					&& lookupTable[20032] == 6
					&& lookupTable[37056] == 7
					&& lookupTable[37782] == 8
					&& lookupTable[58416] == 9
					&& lookupTable[58596] == 10
					&& lookupTable[58597] == 12
					&& lookupTable[58599] == 13
					&& lookupTable[64740] == 14
					&& refPos != NULL
					&& refPos[0] == UINT_MAX
					&& refPos[1] == 70
					&& refPos[2] == 135
					&& refPos[3] == 100
					&& refPos[4] == 79
					&& refPos[5] == 31
					&& refPos[6] == 126
					&& refPos[7] == 108
					&& refPos[8] == 54
					&& refPos[9] == 14
					&& refPos[10] == 62
					&& refPos[11] == 6
					&& refPos[12] == 39
					&& refPos[13] == 117
					&& refPos[14] == 22
					&& numRepeatsPerTuple[3648] == 1
					&& numRepeatsPerTuple[3662] == 1
					&& numRepeatsPerTuple[14649] == 2
					&& numRepeatsPerTuple[17636] == 1
					&& numRepeatsPerTuple[20032] == 1
					&& numRepeatsPerTuple[37056] == 1
					&& numRepeatsPerTuple[37782] == 1
					&& numRepeatsPerTuple[58416] == 1
					&& numRepeatsPerTuple[58596] == 2
					&& numRepeatsPerTuple[58597] == 1
					&& numRepeatsPerTuple[58599] == 1
					&& numRepeatsPerTuple[64740] == 1
		)
		{}
		else
			fail("Incorrect behavior when there are 3 references, each "
					"sequence split into different lines.\n");
		lookupTable2Delete();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references and one of the tuples has more number of repeats than
	 * allowed by the tuple-ignore-threshold. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTT\nACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAA\nCGTACGTATTTACG\n"
				">ref3\nCGTACGTAAA\nATAACGTCGT\nACGTAAACGT\nACGTACGTAA",
				filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		uint *refPos;
		int seedLen = 4, maxHitsPerQry = 2, tupleIgnoreThreshold = 7;

		lookupTable2Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 256 && lookupTable != NULL
				&& lookupTable[14] == 1
				&& lookupTable[48] == 3
				&& lookupTable[57] == 4
				&& lookupTable[63] == 8
				&& lookupTable[64] == 9
				&& lookupTable[68] == 11
				&& lookupTable[78] == 12
				&& lookupTable[144] == 14
				&& lookupTable[147] == 15
				&& lookupTable[150] == 16
				&& lookupTable[192] == 17
				&& lookupTable[229] == 18
				&& lookupTable[231] == 19
				&& lookupTable[252] == 20
				&& refPos != NULL
				&& refPos[0] == UINT_MAX
				&& refPos[1] == 139
				&& refPos[2] == 74
				&& refPos[3] == 14
				&& refPos[4] == 104
				&& refPos[5] == 100
				&& refPos[6] == 83
				&& refPos[7] == 79
				&& refPos[8] == 87
				&& refPos[9] == 126
				&& refPos[10] == 70
				&& refPos[11] == 35
				&& refPos[12] == 135
				&& refPos[13] == 130
				&& refPos[14] == 113
				&& refPos[15] == 58
				&& refPos[16] == 54
				&& refPos[17] == 108
				&& refPos[18] == 39
				&& refPos[19] == 117
				&& refPos[20] == 26
				&& numRepeatsPerTuple[14] == 2
				&& numRepeatsPerTuple[48] == 1
				&& numRepeatsPerTuple[57] == 4
				&& numRepeatsPerTuple[63] == 1
				&& numRepeatsPerTuple[64] == 2
				&& numRepeatsPerTuple[68] == 1
				&& numRepeatsPerTuple[78] == 2
				&& numRepeatsPerTuple[144] == 1
				&& numRepeatsPerTuple[147] == 1
				&& numRepeatsPerTuple[150] == 1
				&& numRepeatsPerTuple[192] == 1
				&& numRepeatsPerTuple[229] == 1
				&& numRepeatsPerTuple[231] == 1
				&& numRepeatsPerTuple[252] == 1)
		{}
		else
			fail("Incorrect behavior when there are 3 referencesand one of the "
					"tuples has more number of repeats than allowed by the "
					"tuple-ignore-threshold.\n");

		lookupTable2Delete();
		remove(refFile);
		preprocessDelete();
	}
}
END_TEST


/**
 * Tests @a lookupTable3Create function.
 */
START_TEST(lookupTable3Create)
{
	/* 1 reference. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		int *refPos;
		char *refIdx;
		int seedLen = 8, maxHitsPerQry = 2, tupleIgnoreThreshold = 15;

		lookupTable3Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refIdx, &refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 65536 && lookupTable != NULL
					&& lookupTable[3662] == 1 && lookupTable[58596] == 2
					&& refIdx != NULL && refIdx[0] == -1 && refIdx[1] == 0
					&& refIdx[2] == 0 && refPos != NULL && refPos[0] == -1
					&& refPos[1] == 8 && refPos[2] == 0
					&& numRepeatsPerTuple[3662] == 1
					&& numRepeatsPerTuple[58596] == 1)
		{}
		else
			fail("Incorrect behavior when there is 1 reference.\n");
		lookupTable3Delete();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references, whole sequence of a given reference is in the same line. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAACGTACGTATTTACG\n"
				">ref3\nCGTACGTAAAATAACGTCGTACGTAAACGTACGTACGTAA", filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		int *refPos;
		char *refIdx;
		int seedLen = 8, maxHitsPerQry = 2, tupleIgnoreThreshold = 15;

		lookupTable3Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refIdx, &refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 65536 && lookupTable != NULL
				&& lookupTable[3648] == 1
				&& lookupTable[3662] == 2
				&& lookupTable[14649] == 3
				&& lookupTable[17636] == 5
				&& lookupTable[20032] == 6
				&& lookupTable[37056] == 7
				&& lookupTable[37782] == 8
				&& lookupTable[58416] == 9
				&& lookupTable[58596] == 10
				&& lookupTable[58597] == 12
				&& lookupTable[58599] == 13
				&& lookupTable[64740] == 14
				&& refIdx != NULL
				&& refIdx[0] == -1
				&& refIdx[1] == 1
				&& refIdx[2] == 2
				&& refIdx[3] == 2
				&& refIdx[4] == 1
				&& refIdx[5] == 0
				&& refIdx[6] == 2
				&& refIdx[7] == 2
				&& refIdx[8] == 1
				&& refIdx[9] == 0
				&& refIdx[10] == 1
				&& refIdx[11] == 0
				&& refIdx[12] == 0
				&& refIdx[13] == 2
				&& refIdx[14] == 0
				&& refPos != NULL
				&& refPos[0] == -1
				&& refPos[1] == 16
				&& refPos[2] == 32
				&& refPos[3] == 0
				&& refPos[4] == 24
				&& refPos[5] == 24
				&& refPos[6] == 24
				&& refPos[7] == 8
				&& refPos[8] == 0
				&& refPos[9] == 8
				&& refPos[10] == 8
				&& refPos[11] == 0
				&& refPos[12] == 32
				&& refPos[13] == 16
				&& refPos[14] == 16
				&& numRepeatsPerTuple[3648] == 1
				&& numRepeatsPerTuple[3662] == 1
				&& numRepeatsPerTuple[14649] == 2
				&& numRepeatsPerTuple[17636] == 1
				&& numRepeatsPerTuple[20032] == 1
				&& numRepeatsPerTuple[37056] == 1
				&& numRepeatsPerTuple[37782] == 1
				&& numRepeatsPerTuple[58416] == 1
				&& numRepeatsPerTuple[58596] == 2
				&& numRepeatsPerTuple[58597] == 1
				&& numRepeatsPerTuple[58599] == 1
				&& numRepeatsPerTuple[64740] == 1
		)
		{}
		else
			fail("Incorrect behavior when there are 3 references, whole "
					"sequence of a given reference is in the same line.\n");
		lookupTable3Delete();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references, each sequence split into different lines. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTT\nACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAA\nCGTACGTATTTACG\n"
				">ref3\nCGTACGTAAA\nATAACGTCGT\nACGTAAACGT\nACGTACGTAA",
				filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		int *refPos;
		char *refIdx;
		int seedLen = 8, maxHitsPerQry = 2, tupleIgnoreThreshold = 4;

		lookupTable3Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refIdx, &refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 65536 && lookupTable != NULL
					&& lookupTable[3648] == 1
					&& lookupTable[3662] == 2
					&& lookupTable[14649] == 3
					&& lookupTable[17636] == 5
					&& lookupTable[20032] == 6
					&& lookupTable[37056] == 7
					&& lookupTable[37782] == 8
					&& lookupTable[58416] == 9
					&& lookupTable[58596] == 10
					&& lookupTable[58597] == 12
					&& lookupTable[58599] == 13
					&& lookupTable[64740] == 14
					&& refIdx != NULL
					&& refIdx[0] == -1
					&& refIdx[1] == 1
					&& refIdx[2] == 2
					&& refIdx[3] == 2
					&& refIdx[4] == 1
					&& refIdx[5] == 0
					&& refIdx[6] == 2
					&& refIdx[7] == 2
					&& refIdx[8] == 1
					&& refIdx[9] == 0
					&& refIdx[10] == 1
					&& refIdx[11] == 0
					&& refIdx[12] == 0
					&& refIdx[13] == 2
					&& refIdx[14] == 0
					&& refPos != NULL
					&& refPos[0] == -1
					&& refPos[1] == 16
					&& refPos[2] == 32
					&& refPos[3] == 0
					&& refPos[4] == 24
					&& refPos[5] == 24
					&& refPos[6] == 24
					&& refPos[7] == 8
					&& refPos[8] == 0
					&& refPos[9] == 8
					&& refPos[10] == 8
					&& refPos[11] == 0
					&& refPos[12] == 32
					&& refPos[13] == 16
					&& refPos[14] == 16
					&& numRepeatsPerTuple[3648] == 1
					&& numRepeatsPerTuple[3662] == 1
					&& numRepeatsPerTuple[14649] == 2
					&& numRepeatsPerTuple[17636] == 1
					&& numRepeatsPerTuple[20032] == 1
					&& numRepeatsPerTuple[37056] == 1
					&& numRepeatsPerTuple[37782] == 1
					&& numRepeatsPerTuple[58416] == 1
					&& numRepeatsPerTuple[58596] == 2
					&& numRepeatsPerTuple[58597] == 1
					&& numRepeatsPerTuple[58599] == 1
					&& numRepeatsPerTuple[64740] == 1
		)
		{}
		else
			fail("Incorrect behavior when there are 3 references, each "
					"sequence split into different lines.\n");
		lookupTable3Delete();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references and one of the tuples has more number of repeats than
	   allowed by the tuple-ignore-threshold. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTT\nACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAA\nCGTACGTATTTACG\n"
				">ref3\nCGTACGTAAA\nATAACGTCGT\nACGTAAACGT\nACGTACGTAA",
				filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		int *refPos;
		char *refIdx;
		int seedLen = 4, maxHitsPerQry = 2, tupleIgnoreThreshold = 7;

		lookupTable3Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refIdx, &refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 256 && lookupTable != NULL
				&& lookupTable[14] == 1
				&& lookupTable[48] == 3
				&& lookupTable[57] == 4
				&& lookupTable[63] == 8
				&& lookupTable[64] == 9
				&& lookupTable[68] == 11
				&& lookupTable[78] == 12
				&& lookupTable[144] == 14
				&& lookupTable[147] == 15
				&& lookupTable[150] == 16
				&& lookupTable[192] == 17
				&& lookupTable[229] == 18
				&& lookupTable[231] == 19
				&& lookupTable[252] == 20
				&& refIdx != NULL
				&& refIdx[0] == -1
				&& refIdx[1] == 2
				&& refIdx[2] == 1
				&& refIdx[3] == 0
				&& refIdx[4] == 2
				&& refIdx[5] == 2
				&& refIdx[6] == 1
				&& refIdx[7] == 1
				&& refIdx[8] == 1
				&& refIdx[9] == 2
				&& refIdx[10] == 1
				&& refIdx[11] == 0
				&& refIdx[12] == 2
				&& refIdx[13] == 2
				&& refIdx[14] == 2
				&& refIdx[15] == 1
				&& refIdx[16] == 1
				&& refIdx[17] == 2
				&& refIdx[18] == 0
				&& refIdx[19] == 2
				&& refIdx[20] == 0
				&& refPos != NULL
				&& refPos[0] == -1
				&& refPos[1] == 36
				&& refPos[2] == 20
				&& refPos[3] == 8
				&& refPos[4] == 4
				&& refPos[5] == 0
				&& refPos[6] == 28
				&& refPos[7] == 24
				&& refPos[8] == 32
				&& refPos[9] == 24
				&& refPos[10] == 16
				&& refPos[11] == 28
				&& refPos[12] == 32
				&& refPos[13] == 28
				&& refPos[14] == 12
				&& refPos[15] == 4
				&& refPos[16] == 0
				&& refPos[17] == 8
				&& refPos[18] == 32
				&& refPos[19] == 16
				&& refPos[20] == 20
				&& numRepeatsPerTuple[14] == 2
				&& numRepeatsPerTuple[48] == 1
				&& numRepeatsPerTuple[57] == 4
				&& numRepeatsPerTuple[63] == 1
				&& numRepeatsPerTuple[64] == 2
				&& numRepeatsPerTuple[68] == 1
				&& numRepeatsPerTuple[78] == 2
				&& numRepeatsPerTuple[144] == 1
				&& numRepeatsPerTuple[147] == 1
				&& numRepeatsPerTuple[150] == 1
				&& numRepeatsPerTuple[192] == 1
				&& numRepeatsPerTuple[229] == 1
				&& numRepeatsPerTuple[231] == 1
				&& numRepeatsPerTuple[252] == 1)
		{}
		else
			fail("Incorrect behavior when there are 3 referencesand one of the "
					"tuples has more number of repeats than allowed by the "
					"tuple-ignore-threshold.\n");

		lookupTable3Delete();
		remove(refFile);
		preprocessDelete();
	}
}
END_TEST


/**
 * Tests @a lookupTable3Create function.
 */
START_TEST(lookupTable4Create)
{
	/* 1 reference. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		int *refPos;
		char *refIdx;
		int seedLen = 8, maxHitsPerQry = 2, tupleIgnoreThreshold = 15;

		lookupTable4Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refIdx, &refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 65536 && lookupTable != NULL
					&& lookupTable[3662] == 1 && lookupTable[58596] == 2
					&& refIdx != NULL && refIdx[0] == -1 && refIdx[1] == 0
					&& refIdx[2] == 0 && refPos != NULL && refPos[0] == -1
					&& refPos[1] == 8 && refPos[2] == 0
					&& numRepeatsPerTuple[3662] == 1
					&& numRepeatsPerTuple[58596] == 1)
		{}
		else
			fail("Incorrect behavior when there is 1 reference.\n");
		lookupTable4Delete();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references, whole sequence of a given reference is in the same line. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAACGTACGTATTTACG\n"
				">ref3\nCGTACGTAAAATAACGTCGTACGTAAACGTACGTACGTAA", filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		int *refPos;
		char *refIdx;
		int seedLen = 8, maxHitsPerQry = 2, tupleIgnoreThreshold = 15;

		lookupTable4Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refIdx, &refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 65536 && lookupTable != NULL
				&& lookupTable[3648] == 1
				&& lookupTable[3662] == 2
				&& lookupTable[14649] == 3
				&& lookupTable[17636] == 5
				&& lookupTable[20032] == 6
				&& lookupTable[37056] == 7
				&& lookupTable[37782] == 8
				&& lookupTable[58416] == 9
				&& lookupTable[58596] == 10
				&& lookupTable[58597] == 12
				&& lookupTable[58599] == 13
				&& lookupTable[64740] == 14
				&& refIdx != NULL
				&& refIdx[0] == -1
				&& refIdx[1] == 1
				&& refIdx[2] == 2
				&& refIdx[3] == 2
				&& refIdx[4] == 1
				&& refIdx[5] == 0
				&& refIdx[6] == 2
				&& refIdx[7] == 2
				&& refIdx[8] == 1
				&& refIdx[9] == 0
				&& refIdx[10] == 1
				&& refIdx[11] == 0
				&& refIdx[12] == 0
				&& refIdx[13] == 2
				&& refIdx[14] == 0
				&& refPos != NULL
				&& refPos[0] == -1
				&& refPos[1] == 16
				&& refPos[2] == 32
				&& refPos[3] == 0
				&& refPos[4] == 24
				&& refPos[5] == 24
				&& refPos[6] == 24
				&& refPos[7] == 8
				&& refPos[8] == 0
				&& refPos[9] == 8
				&& refPos[10] == 8
				&& refPos[11] == 0
				&& refPos[12] == 32
				&& refPos[13] == 16
				&& refPos[14] == 16
				&& numRepeatsPerTuple[3648] == 1
				&& numRepeatsPerTuple[3662] == 1
				&& numRepeatsPerTuple[14649] == 2
				&& numRepeatsPerTuple[17636] == 1
				&& numRepeatsPerTuple[20032] == 1
				&& numRepeatsPerTuple[37056] == 1
				&& numRepeatsPerTuple[37782] == 1
				&& numRepeatsPerTuple[58416] == 1
				&& numRepeatsPerTuple[58596] == 2
				&& numRepeatsPerTuple[58597] == 1
				&& numRepeatsPerTuple[58599] == 1
				&& numRepeatsPerTuple[64740] == 1
		)
		{}
		else
			fail("Incorrect behavior when there are 3 references, whole "
					"sequence of a given reference is in the same line.\n");
		lookupTable4Delete();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references, each sequence split into different lines. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTT\nACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAA\nCGTACGTATTTACG\n"
				">ref3\nCGTACGTAAA\nATAACGTCGT\nACGTAAACGT\nACGTACGTAA",
				filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		int *refPos;
		char *refIdx;
		int seedLen = 8, maxHitsPerQry = 2, tupleIgnoreThreshold = 4;

		lookupTable4Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refIdx, &refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 65536 && lookupTable != NULL
					&& lookupTable[3648] == 1
					&& lookupTable[3662] == 2
					&& lookupTable[14649] == 3
					&& lookupTable[17636] == 5
					&& lookupTable[20032] == 6
					&& lookupTable[37056] == 7
					&& lookupTable[37782] == 8
					&& lookupTable[58416] == 9
					&& lookupTable[58596] == 10
					&& lookupTable[58597] == 12
					&& lookupTable[58599] == 13
					&& lookupTable[64740] == 14
					&& refIdx != NULL
					&& refIdx[0] == -1
					&& refIdx[1] == 1
					&& refIdx[2] == 2
					&& refIdx[3] == 2
					&& refIdx[4] == 1
					&& refIdx[5] == 0
					&& refIdx[6] == 2
					&& refIdx[7] == 2
					&& refIdx[8] == 1
					&& refIdx[9] == 0
					&& refIdx[10] == 1
					&& refIdx[11] == 0
					&& refIdx[12] == 0
					&& refIdx[13] == 2
					&& refIdx[14] == 0
					&& refPos != NULL
					&& refPos[0] == -1
					&& refPos[1] == 16
					&& refPos[2] == 32
					&& refPos[3] == 0
					&& refPos[4] == 24
					&& refPos[5] == 24
					&& refPos[6] == 24
					&& refPos[7] == 8
					&& refPos[8] == 0
					&& refPos[9] == 8
					&& refPos[10] == 8
					&& refPos[11] == 0
					&& refPos[12] == 32
					&& refPos[13] == 16
					&& refPos[14] == 16
					&& numRepeatsPerTuple[3648] == 1
					&& numRepeatsPerTuple[3662] == 1
					&& numRepeatsPerTuple[14649] == 2
					&& numRepeatsPerTuple[17636] == 1
					&& numRepeatsPerTuple[20032] == 1
					&& numRepeatsPerTuple[37056] == 1
					&& numRepeatsPerTuple[37782] == 1
					&& numRepeatsPerTuple[58416] == 1
					&& numRepeatsPerTuple[58596] == 2
					&& numRepeatsPerTuple[58597] == 1
					&& numRepeatsPerTuple[58599] == 1
					&& numRepeatsPerTuple[64740] == 1
		)
		{}
		else
			fail("Incorrect behavior when there are 3 references, each "
					"sequence split into different lines.\n");
		lookupTable4Delete();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references and one of the tuples has more number of repeats than
	   allowed by the tuple-ignore-threshold. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTT\nACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAA\nCGTACGTATTTACG\n"
				">ref3\nCGTACGTAAA\nATAACGTCGT\nACGTAAACGT\nACGTACGTAA",
				filePtr);
		fclose(filePtr);

		int *lookupTable, numDistinctTuples, *numRepeatsPerTuple;
		int *refPos;
		char *refIdx;
		int seedLen = 4, maxHitsPerQry = 2, tupleIgnoreThreshold = 7;

		lookupTable4Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
				&refIdx, &refPos, &numDistinctTuples, &numRepeatsPerTuple,
				tupleIgnoreThreshold);

		if (numDistinctTuples == 256 && lookupTable != NULL
				&& lookupTable[14] == 1
				&& lookupTable[48] == 3
				&& lookupTable[57] == 4
				&& lookupTable[63] == 8
				&& lookupTable[64] == 9
				&& lookupTable[68] == 11
				&& lookupTable[78] == 12
				&& lookupTable[144] == 14
				&& lookupTable[147] == 15
				&& lookupTable[150] == 16
				&& lookupTable[192] == 17
				&& lookupTable[229] == 18
				&& lookupTable[231] == 19
				&& lookupTable[252] == 20
				&& refIdx != NULL
				&& refIdx[0] == -1
				&& refIdx[1] == 2
				&& refIdx[2] == 1
				&& refIdx[3] == 0
				&& refIdx[4] == 2
				&& refIdx[5] == 2
				&& refIdx[6] == 1
				&& refIdx[7] == 1
				&& refIdx[8] == 1
				&& refIdx[9] == 2
				&& refIdx[10] == 1
				&& refIdx[11] == 0
				&& refIdx[12] == 2
				&& refIdx[13] == 2
				&& refIdx[14] == 2
				&& refIdx[15] == 1
				&& refIdx[16] == 1
				&& refIdx[17] == 2
				&& refIdx[18] == 0
				&& refIdx[19] == 2
				&& refIdx[20] == 0
				&& refPos != NULL
				&& refPos[0] == -1
				&& refPos[1] == 36
				&& refPos[2] == 20
				&& refPos[3] == 8
				&& refPos[4] == 4
				&& refPos[5] == 0
				&& refPos[6] == 28
				&& refPos[7] == 24
				&& refPos[8] == 32
				&& refPos[9] == 24
				&& refPos[10] == 16
				&& refPos[11] == 28
				&& refPos[12] == 32
				&& refPos[13] == 28
				&& refPos[14] == 12
				&& refPos[15] == 4
				&& refPos[16] == 0
				&& refPos[17] == 8
				&& refPos[18] == 32
				&& refPos[19] == 16
				&& refPos[20] == 20
				&& numRepeatsPerTuple[14] == 2
				&& numRepeatsPerTuple[48] == 1
				&& numRepeatsPerTuple[57] == 4
				&& numRepeatsPerTuple[63] == 1
				&& numRepeatsPerTuple[64] == 2
				&& numRepeatsPerTuple[68] == 1
				&& numRepeatsPerTuple[78] == 2
				&& numRepeatsPerTuple[144] == 1
				&& numRepeatsPerTuple[147] == 1
				&& numRepeatsPerTuple[150] == 1
				&& numRepeatsPerTuple[192] == 1
				&& numRepeatsPerTuple[229] == 1
				&& numRepeatsPerTuple[231] == 1
				&& numRepeatsPerTuple[252] == 1)
		{}
		else
			fail("Incorrect behavior when there are 3 referencesand one of the "
					"tuples has more number of repeats than allowed by the "
					"tuple-ignore-threshold.\n");

		lookupTable4Delete();
		remove(refFile);
		preprocessDelete();
	}
}
END_TEST


/**
 * Tests @a lookupTable2MapQry function.
 */
START_TEST(lookupTable2MapQry)
{
	/* 1 reference. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);
		refMapCreate(refFile);
		int seedLength = 8, tupleIgnoreThreshold = 15;
		uint maxNumHits = 10;
		lookupTable2Create(refFile, seedLength, maxNumHits,
				tupleIgnoreThreshold);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		HitList *result = lookupTable2MapQry(query, seqLength, isRevComplement);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 6, 8);
		if (result != NULL
				&& hitListCmpNodes(result, node1) == 0
				&& hitListCmpNodes(result->next, node2) == 0)
		{}
		else
			fail("Incorrect behavior when there is 1 reference (case 1).\n");
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);

		HitList *node3 = hitListCreateNode(0, -2, 0);
		isRevComplement = 1;
		result = lookupTable2MapQry(query, seqLength, isRevComplement);
		if (result != NULL && hitListCmpNodes(result, node3) == 0)
		{}
		else
			fail("Incorrect behavior when there is 1 reference (case 2).\n");
		hitListDelete(result);
		qryDelete(query);
		lookupTable2Delete();
		hitListDelete(node3);
		refMapFree();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references, sequence of each reference in one line. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTACGTACGTACGT\n>ref2\nGCCGTACGACGTACGT\n"
				">ref3\nCGTACGTAAAATAA", filePtr);
		fclose(filePtr);
		refMapCreate(refFile);

		int seedLength = 8, tupleIgnoreThres = 15;
		uint maxNumHits = 2;
		lookupTable2Create(refFile, seedLength, maxNumHits, tupleIgnoreThres);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		HitList *result = lookupTable2MapQry(query, seqLength, isRevComplement);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 8, 8);
		HitList *node3 = hitListCreateNode(1, 8, 8);
		HitList *node4 = hitListCreateNode(2, -1, 0);
		if (result != NULL && hitListGetSize(result) == 2
				&& ((hitListCmpNodes(result, node1) == 0
						&& hitListCmpNodes(result->next, node2) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node4) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node4) == 0)
						|| (hitListCmpNodes(result, node3) == 0
								&& hitListCmpNodes(result->next, node4) == 0)))
		{ }
		else
			fail("Incorrect behavior when there are 3 references, sequence "
					"of each reference in 1 line (case 1).\n");
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		hitListDelete(node4);

		isRevComplement = 1;
		result = lookupTable2MapQry(query, seqLength, isRevComplement);
		node1 = hitListCreateNode(0, -2, 0);
		node2 = hitListCreateNode(0, 6, 8);
		node3 = hitListCreateNode(1, 6, 8);
		if (result != NULL && hitListGetSize(result) == 2
				&& ((hitListCmpNodes(result, node1) == 0
						&& hitListCmpNodes(result->next, node2) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node3) == 0)))
		{}
		else
			fail("Incorrect behavior when there are 3 references, sequence "
					"of each reference in 1 line (case 2).\n");
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		qryDelete(query);
		lookupTable2Delete();
		refMapFree();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references, sequence of each reference split into different lines. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTACGTACGTACGT\n>ref2\nGCCGTACGACGTACGT\n"
				">ref3\nCGTACGTAAAATAA", filePtr);
		fclose(filePtr);
		refMapCreate(refFile);

		int seedLength = 8, tupleIgnoreThres = 15;
		uint maxNumHits = 2;
		lookupTable2Create(refFile, seedLength, maxNumHits, tupleIgnoreThres);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		HitList *result = lookupTable2MapQry(query, seqLength, isRevComplement);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 8, 8);
		HitList *node3 = hitListCreateNode(1, 8, 8);
		HitList *node4 = hitListCreateNode(2, -1, 0);
		if (result != NULL && hitListGetSize(result) == 2
				&& ((hitListCmpNodes(result, node1) == 0
						&& hitListCmpNodes(result->next, node2) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node4) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node4) == 0)
						|| (hitListCmpNodes(result, node3) == 0
								&& hitListCmpNodes(result->next, node4) == 0)))
		{ }
		else
			fail("Incorrect behavior when there are 3 references, sequence "
					"of each reference in 1 line (case 1).\n");
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		hitListDelete(node4);

		isRevComplement = 1;
		result = lookupTable2MapQry(query, seqLength, isRevComplement);
		node1 = hitListCreateNode(0, -2, 0);
		node2 = hitListCreateNode(0, 6, 8);
		node3 = hitListCreateNode(1, 6, 8);
		if (result != NULL && hitListGetSize(result) == 2
				&& ((hitListCmpNodes(result, node1) == 0
						&& hitListCmpNodes(result->next, node2) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node3) == 0)))
		{}
		else
			fail("Incorrect behavior when there are 3 references, sequence "
					"of each reference in 1 line (case 2).\n");
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		qryDelete(query);
		lookupTable2Delete();
		refMapFree();
		remove(refFile);
		preprocessDelete();
	}
}
END_TEST


/**
 * Tests @a lookupTable3MapQry function.
 */
START_TEST(lookupTable3MapQry)
{
	/* 1 reference. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);
		int seedLength = 8, tupleIgnoreThreshold = 15;
		uint maxNumHits = 10;
		lookupTable3Create(refFile, seedLength, maxNumHits,
				tupleIgnoreThreshold);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		HitList *result = lookupTable3MapQry(query, seqLength, isRevComplement);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 6, 8);
		if (result != NULL
				&& hitListCmpNodes(result, node1) == 0
				&& hitListCmpNodes(result->next, node2) == 0)
		{}
		else
			fail("Incorrect behavior when there is 1 reference (case 1).\n");
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);

		HitList *node3 = hitListCreateNode(0, -2, 0);
		isRevComplement = 1;
		result = lookupTable3MapQry(query, seqLength, isRevComplement);
		if (result != NULL && hitListCmpNodes(result, node3) == 0)
		{}
		else
			fail("Incorrect behavior when there is 1 reference (case 2).\n");
		hitListDelete(result);
		qryDelete(query);
		lookupTable3Delete();
		hitListDelete(node3);
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references, sequence of each reference in one line. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTACGTACGTACGT\n>ref2\nGCCGTACGACGTACGT\n"
				">ref3\nCGTACGTAAAATAA", filePtr);
		fclose(filePtr);

		int seedLength = 8, tupleIgnoreThres = 15;
		uint maxNumHits = 2;
		lookupTable3Create(refFile, seedLength, maxNumHits, tupleIgnoreThres);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		HitList *result = lookupTable3MapQry(query, seqLength, isRevComplement);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 8, 8);
		HitList *node3 = hitListCreateNode(1, 8, 8);
		HitList *node4 = hitListCreateNode(2, -1, 0);
		if (result != NULL && hitListGetSize(result) == 2
				&& ((hitListCmpNodes(result, node1) == 0
						&& hitListCmpNodes(result->next, node2) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node4) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node4) == 0)
						|| (hitListCmpNodes(result, node3) == 0
								&& hitListCmpNodes(result->next, node4) == 0)))
		{ }
		else
			fail("Incorrect behavior when there are 3 references, sequence "
					"of each reference in 1 line (case 1).\n");
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		hitListDelete(node4);

		isRevComplement = 1;
		result = lookupTable3MapQry(query, seqLength, isRevComplement);
		node1 = hitListCreateNode(0, -2, 0);
		node2 = hitListCreateNode(0, 6, 8);
		node3 = hitListCreateNode(1, 6, 8);
		if (result != NULL && hitListGetSize(result) == 2
				&& ((hitListCmpNodes(result, node1) == 0
						&& hitListCmpNodes(result->next, node2) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node3) == 0)))
		{}
		else
			fail("Incorrect behavior when there are 3 references, sequence "
					"of each reference in 1 line (case 2).\n");
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		qryDelete(query);
		lookupTable3Delete();
		remove(refFile);
		preprocessDelete();
	}

	/* 3 references, sequence of each reference split into different lines. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTACGTACGTACGT\n>ref2\nGCCGTACGACGTACGT\n"
				">ref3\nCGTACGTAAAATAA", filePtr);
		fclose(filePtr);

		int seedLength = 8, tupleIgnoreThres = 15;
		uint maxNumHits = 2;
		lookupTable3Create(refFile, seedLength, maxNumHits, tupleIgnoreThres);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		HitList *result = lookupTable3MapQry(query, seqLength, isRevComplement);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 8, 8);
		HitList *node3 = hitListCreateNode(1, 8, 8);
		HitList *node4 = hitListCreateNode(2, -1, 0);
		if (result != NULL && hitListGetSize(result) == 2
				&& ((hitListCmpNodes(result, node1) == 0
						&& hitListCmpNodes(result->next, node2) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node4) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node4) == 0)
						|| (hitListCmpNodes(result, node3) == 0
								&& hitListCmpNodes(result->next, node4) == 0)))
		{ }
		else
			fail("Incorrect behavior when there are 3 references, sequence "
					"of each reference in 1 line (case 1).\n");
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		hitListDelete(node4);

		isRevComplement = 1;
		result = lookupTable3MapQry(query, seqLength, isRevComplement);
		node1 = hitListCreateNode(0, -2, 0);
		node2 = hitListCreateNode(0, 6, 8);
		node3 = hitListCreateNode(1, 6, 8);
		if (result != NULL && hitListGetSize(result) == 2
				&& ((hitListCmpNodes(result, node1) == 0
						&& hitListCmpNodes(result->next, node2) == 0)
						|| (hitListCmpNodes(result, node1) == 0
								&& hitListCmpNodes(result->next, node3) == 0)
						|| (hitListCmpNodes(result, node2) == 0
								&& hitListCmpNodes(result->next, node3) == 0)))
		{}
		else
			fail("Incorrect behavior when there are 3 references, sequence "
					"of each reference in 1 line (case 2).\n");
		hitListDelete(result);
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		qryDelete(query);
		lookupTable3Delete();
		remove(refFile);
		preprocessDelete();
	}
}
END_TEST


/**
 * Tests @a lookupTable3MapQry function.
 */
START_TEST(lookupTable4MapQry)
{
	/* 1 reference. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);
		int seedLength = 8, tupleIgnoreThreshold = 15;
		uint maxNumHits = 10;
		lookupTable4Create(refFile, seedLength, maxNumHits,
				tupleIgnoreThreshold);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
		int *shift = (int *) calloc(maxNumHits, sizeof(int));
		int *refPos = (int *) calloc(maxNumHits, sizeof(int));
		int numHits = lookupTable4MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		if (numHits != 2 || refIdx[0] != 0 || refIdx[1] != 0 || shift[0] != 0
				|| shift[1] != 6 || refPos[0] != 0 || refPos[1] != 8)
			fail("Incorrect behavior when there is 1 reference (forward strand)."
					"\n");

		isRevComplement = 1;
		numHits = lookupTable4MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		if (numHits != 1 || refIdx[0] != 0 || shift[0] != -2 || refPos[0] != 0)
			fail("Incorrect behavior when there is 1 reference (reverse "
					"complement).\n");

		qryDelete(query);
		lookupTable4Delete();
		remove(refFile);
		free(refIdx);
		free(shift);
		free(refPos);
		preprocessDelete();
	}

	/* 3 references, sequence of each reference in one line. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTACGTACGTACGT\n>ref2\nGCCGTACGACGTACGT\n"
				">ref3\nCGTACGTAAAATAA", filePtr);
		fclose(filePtr);

		int seedLength = 8, tupleIgnoreThres = 15;
		uint maxNumHits = 2;
		lookupTable4Create(refFile, seedLength, maxNumHits, tupleIgnoreThres);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
		int *shift = (int *) calloc(maxNumHits, sizeof(int));
		int *refPos = (int *) calloc(maxNumHits, sizeof(int));
		Hit *hit1 = createHit(0, 0, 0);
		Hit *hit2 = createHit(0, 8, 8);
		Hit *hit3 = createHit(1, 8, 8);
		Hit *hit4 = createHit(2, -1, 0);
		int numHits = lookupTable4MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		Hit *result1 = createHit(refIdx[0], shift[0], refPos[0]);
		Hit *result2 = createHit(refIdx[1], shift[1], refPos[1]);
		if (numHits == 2 && ((areHitsEqual(result1, hit1) == 0
						&& areHitsEqual(result2, hit2) == 0)
						|| (areHitsEqual(result1, hit2) == 0
								&& areHitsEqual(result2, hit1) == 0)
						|| (areHitsEqual(result1, hit1) == 0
								&& areHitsEqual(result2, hit3) == 0)
						|| (areHitsEqual(result1, hit3) == 0
								&& areHitsEqual(result2, hit1) == 0)
						|| (areHitsEqual(result1, hit1) == 0
								&& areHitsEqual(result2, hit4) == 0)
						|| (areHitsEqual(result1, hit4) == 0
								&& areHitsEqual(result2, hit1) == 0)
						|| (areHitsEqual(result1, hit2) == 0
								&& areHitsEqual(result2, hit3) == 0)
						|| (areHitsEqual(result1, hit3) == 0
								&& areHitsEqual(result2, hit2) == 0)
						|| (areHitsEqual(result1, hit2) == 0
								&& areHitsEqual(result2, hit4) == 0)
						|| (areHitsEqual(result1, hit4) == 0
								&& areHitsEqual(result2, hit2) == 0)
						|| (areHitsEqual(result1, hit3) == 0
								&& areHitsEqual(result2, hit4) == 0)
						|| (areHitsEqual(result1, hit4) == 0
								&& areHitsEqual(result2, hit3) == 0)))
		{ }
		else
			fail("Incorrect behavior when there are 3 references, sequence "
					"of each reference in 1 line (case 1).\n");
		free(hit1);
		free(hit2);
		free(hit3);
		free(hit4);
		free(result1);
		free(result2);

		isRevComplement = 1;
		numHits = lookupTable4MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		hit1 = createHit(0, -2, 0);
		hit2 = createHit(0, 6, 8);
		hit3 = createHit(1, 6, 8);
		result1 = createHit(refIdx[0], shift[0], refPos[0]);
		result2 = createHit(refIdx[1], shift[1], refPos[1]);
		if (numHits == 2 && ((areHitsEqual(result1, hit1) == 0
						&& areHitsEqual(result2, hit2) == 0)
						|| (areHitsEqual(result1, hit2) == 0
								&& areHitsEqual(result2, hit1) == 0)
						|| (areHitsEqual(result1, hit1) == 0
								&& areHitsEqual(result2, hit3) == 0)
						|| (areHitsEqual(result1, hit3) == 0
								&& areHitsEqual(result2, hit1) == 0)
						|| (areHitsEqual(result1, hit2) == 0
								&& areHitsEqual(result2, hit3) == 0)
						|| (areHitsEqual(result1, hit3) == 0
								&& areHitsEqual(result2, hit2) == 0)))
		{}
		else
			fail("Incorrect behavior when there are 3 references, sequence "
					"of each reference in 1 line (case 2).\n");

		qryDelete(query);
		lookupTable4Delete();
		remove(refFile);
		free(refIdx);
		free(shift);
		free(refPos);
		free(hit1);
		free(hit2);
		free(hit3);
		free(result1);
		free(result2);
		preprocessDelete();
	}

	/* 3 references, sequence of each reference split into different lines. */
	{
		preprocessCreate();
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTACGTACGTACGT\n>ref2\nGCCGTACGACGTACGT\n"
				">ref3\nCGTACGTAAAATAA", filePtr);
		fclose(filePtr);

		int seedLength = 8, tupleIgnoreThres = 15;
		uint maxNumHits = 2;
		lookupTable4Create(refFile, seedLength, maxNumHits, tupleIgnoreThres);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
		int *shift = (int *) calloc(maxNumHits, sizeof(int));
		int *refPos = (int *) calloc(maxNumHits, sizeof(int));
		int numHits = lookupTable4MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		Hit *hit1 = createHit(0, 0, 0);
		Hit *hit2 = createHit(0, 8, 8);
		Hit *hit3 = createHit(1, 8, 8);
		Hit *hit4 = createHit(2, -1, 0);
		Hit *result1 = createHit(refIdx[0], shift[0], refPos[0]);
		Hit *result2 = createHit(refIdx[1], shift[1], refPos[1]);
		if (numHits == 2 && ((areHitsEqual(result1, hit1) == 0
						&& areHitsEqual(result2, hit2) == 0)
						|| (areHitsEqual(result1, hit2) == 0
								&& areHitsEqual(result2, hit1) == 0)
						|| (areHitsEqual(result1, hit1) == 0
								&& areHitsEqual(result2, hit3) == 0)
						|| (areHitsEqual(result1, hit3) == 0
								&& areHitsEqual(result2, hit1) == 0)
						|| (areHitsEqual(result1, hit1) == 0
								&& areHitsEqual(result2, hit4) == 0)
						|| (areHitsEqual(result1, hit4) == 0
								&& areHitsEqual(result2, hit1) == 0)
						|| (areHitsEqual(result1, hit2) == 0
								&& areHitsEqual(result2, hit3) == 0)
						|| (areHitsEqual(result1, hit3) == 0
								&& areHitsEqual(result2, hit2) == 0)
						|| (areHitsEqual(result1, hit2) == 0
								&& areHitsEqual(result2, hit4) == 0)
						|| (areHitsEqual(result1, hit4) == 0
								&& areHitsEqual(result2, hit2) == 0)
						|| (areHitsEqual(result1, hit3) == 0
								&& areHitsEqual(result2, hit4) == 0)
						|| (areHitsEqual(result1, hit4) == 0
								&& areHitsEqual(result2, hit3) == 0)))
		{ }
		else
			fail("Incorrect behavior when there are 3 references, sequence of "
					"each reference split into different lines (case 1).\n");
		free(result1);
		free(result2);
		free(hit1);
		free(hit2);
		free(hit3);
		free(hit4);

		isRevComplement = 1;
		numHits = lookupTable4MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		hit1 = createHit(0, -2, 0);
		hit2 = createHit(0, 6, 8);
		hit3 = createHit(1, 6, 8);
		result1 = createHit(refIdx[0], shift[0], refPos[0]);
		result2 = createHit(refIdx[1], shift[1], refPos[1]);
		if (numHits == 2 && ((areHitsEqual(result1, hit1) == 0
						&& areHitsEqual(result2, hit2) == 0)
						|| (areHitsEqual(result1, hit2) == 0
								&& areHitsEqual(result2, hit1) == 0)
						|| (areHitsEqual(result1, hit1) == 0
								&& areHitsEqual(result2, hit3) == 0)
						|| (areHitsEqual(result1, hit3) == 0
								&& areHitsEqual(result2, hit1) == 0)
						|| (areHitsEqual(result1, hit2) == 0
								&& areHitsEqual(result2, hit3) == 0)
						|| (areHitsEqual(result1, hit3) == 0
								&& areHitsEqual(result2, hit2) == 0)))
		{}
		else
			fail("Incorrect behavior when there are 3 references, sequence of "
					"each reference split into different lines (case 2).\n");
		free(result1);
		free(result2);
		free(hit1);
		free(hit2);
		free(hit3);
		qryDelete(query);
		lookupTable4Delete();
		remove(refFile);
		free(refIdx);
		free(shift);
		free(refPos);
		preprocessDelete();
	}
}
END_TEST


/**
 * Returns a 'Hit' data structure.
 *
 * @param	refIdx	Reference index
 * @param	shift		Shift.
 * @param	refPos	Reference position.
 * @return	Creates a 'hit' data structure.
 */
static Hit *createHit(char refIdx, int shift, int refPos)
{
	Hit *hit = (Hit *) malloc(sizeof(Hit));
	hit->refIdx = refIdx;
	hit->shift = shift;
	hit->refPos = refPos;
	return hit;
}


/**
 * Returns 0 if hits are equal, 1 otherwise.
 *
 * @param	hit1	First hit.
 * @param	hit2	Second hit.
 * @return	0 if hits are equal, 1 otherwise.
 */
static int areHitsEqual(const Hit *hit1, const Hit *hit2)
{
	if (hit1->refIdx == hit2->refIdx && hit1->shift == hit2->shift
			&& hit1->refPos == hit2->refPos)
		return 0;
	else
		return 1;
}


/**
 * Prints the given 'Hit'.
 *
 * @param hit The 'Hit' to be printed.
 */
static void printHit(const Hit *hit)
{
	fprintf(stderr, "refIdx: %d; shift: %d; refPos: %d\n", hit->refIdx,
			hit->shift, hit->refPos);
}


/**
 * Creates test suite.
 */
Suite *lookupTableSuite(void)
{
	Suite *s = suite_create("lookupTable");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, lookupTableMapQry);
	tcase_add_test(testCaseCore, lookupTable2MapQry);
	tcase_add_test(testCaseCore, lookupTable2Create);
	tcase_add_test(testCaseCore, lookupTable3Create);
	tcase_add_test(testCaseCore, lookupTable4Create);
	tcase_add_test(testCaseCore, lookupTable3MapQry);
	tcase_add_test(testCaseCore, lookupTable4MapQry);
	suite_add_tcase (s, testCaseCore);

	return s;
}
