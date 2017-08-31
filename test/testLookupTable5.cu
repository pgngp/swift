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
#include "../src/common.h"
#include "testLookupTable.h"
#include "testLookupTable5.h"
#include "../src/lookupTable5.h"
#include "../src/refMap.h"
#include "../src/common.h"
#include "../src/preprocess.h"
#include <limits.h>


/**
 * Tests @a lookupTable5Create function.
 */
START_TEST(lookupTable5Create)
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

		lookupTable5Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
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
		lookupTable5Delete();
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

		lookupTable5Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
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
		lookupTable5Delete();
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

		lookupTable5Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
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
		lookupTable5Delete();
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

		lookupTable5Create_wrap(refFile, seedLen, maxHitsPerQry, &lookupTable,
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

		lookupTable5Delete();
		remove(refFile);
		preprocessDelete();
	}
}
END_TEST


/**
 * Tests @a lookupTable5MapQry function.
 */
START_TEST(lookupTable5MapQry)
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
		int totalTuples = 0;
		lookupTable5Create(refFile, seedLength, maxNumHits,
				tupleIgnoreThreshold, &totalTuples);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
		int *shift = (int *) calloc(maxNumHits, sizeof(int));
		int *refPos = (int *) calloc(maxNumHits, sizeof(int));
		int numHits = lookupTable5MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		if (numHits != 2 || refIdx[0] != 0 || refIdx[1] != 0 || shift[0] != 0
				|| shift[1] != 6 || refPos[0] != 0 || refPos[1] != 8)
			fail("Incorrect behavior when there is 1 reference (forward strand)."
					"\n");

		isRevComplement = 1;
		numHits = lookupTable5MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		if (numHits != 1 || refIdx[0] != 0 || shift[0] != -2 || refPos[0] != 0)
			fail("Incorrect behavior when there is 1 reference (reverse "
					"complement).\n");

		qryDelete(query);
		lookupTable5Delete();
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
		int totalTuples = 0;
		lookupTable5Create(refFile, seedLength, maxNumHits, tupleIgnoreThres,
				&totalTuples);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
		int *shift = (int *) calloc(maxNumHits, sizeof(int));
		int *refPos = (int *) calloc(maxNumHits, sizeof(int));
		Hit5 *hit1 = createHit5(0, 0, 0);
		Hit5 *hit2 = createHit5(0, 8, 8);
		Hit5 *hit3 = createHit5(1, 8, 8);
		Hit5 *hit4 = createHit5(2, -1, 0);
		int numHits = lookupTable5MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		Hit5 *result1 = createHit5(refIdx[0], shift[0], refPos[0]);
		Hit5 *result2 = createHit5(refIdx[1], shift[1], refPos[1]);
		if (numHits == 2 && ((areHitsEqual5(result1, hit1) == 0
						&& areHitsEqual5(result2, hit2) == 0)
						|| (areHitsEqual5(result1, hit2) == 0
								&& areHitsEqual5(result2, hit1) == 0)
						|| (areHitsEqual5(result1, hit1) == 0
								&& areHitsEqual5(result2, hit3) == 0)
						|| (areHitsEqual5(result1, hit3) == 0
								&& areHitsEqual5(result2, hit1) == 0)
						|| (areHitsEqual5(result1, hit1) == 0
								&& areHitsEqual5(result2, hit4) == 0)
						|| (areHitsEqual5(result1, hit4) == 0
								&& areHitsEqual5(result2, hit1) == 0)
						|| (areHitsEqual5(result1, hit2) == 0
								&& areHitsEqual5(result2, hit3) == 0)
						|| (areHitsEqual5(result1, hit3) == 0
								&& areHitsEqual5(result2, hit2) == 0)
						|| (areHitsEqual5(result1, hit2) == 0
								&& areHitsEqual5(result2, hit4) == 0)
						|| (areHitsEqual5(result1, hit4) == 0
								&& areHitsEqual5(result2, hit2) == 0)
						|| (areHitsEqual5(result1, hit3) == 0
								&& areHitsEqual5(result2, hit4) == 0)
						|| (areHitsEqual5(result1, hit4) == 0
								&& areHitsEqual5(result2, hit3) == 0)))
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
		numHits = lookupTable5MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		hit1 = createHit5(0, -2, 0);
		hit2 = createHit5(0, 6, 8);
		hit3 = createHit5(1, 6, 8);
		result1 = createHit5(refIdx[0], shift[0], refPos[0]);
		result2 = createHit5(refIdx[1], shift[1], refPos[1]);
		if (numHits == 2 && ((areHitsEqual5(result1, hit1) == 0
						&& areHitsEqual5(result2, hit2) == 0)
						|| (areHitsEqual5(result1, hit2) == 0
								&& areHitsEqual5(result2, hit1) == 0)
						|| (areHitsEqual5(result1, hit1) == 0
								&& areHitsEqual5(result2, hit3) == 0)
						|| (areHitsEqual5(result1, hit3) == 0
								&& areHitsEqual5(result2, hit1) == 0)
						|| (areHitsEqual5(result1, hit2) == 0
								&& areHitsEqual5(result2, hit3) == 0)
						|| (areHitsEqual5(result1, hit3) == 0
								&& areHitsEqual5(result2, hit2) == 0)))
		{}
		else
			fail("Incorrect behavior when there are 3 references, sequence "
					"of each reference in 1 line (case 2).\n");

		qryDelete(query);
		lookupTable5Delete();
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
		int totalTuples = 0;
		lookupTable5Create(refFile, seedLength, maxNumHits, tupleIgnoreThres,
				&totalTuples);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		char *refIdx = (char *) calloc(maxNumHits, sizeof(char));
		int *shift = (int *) calloc(maxNumHits, sizeof(int));
		int *refPos = (int *) calloc(maxNumHits, sizeof(int));
		int numHits = lookupTable5MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		Hit5 *hit1 = createHit5(0, 0, 0);
		Hit5 *hit2 = createHit5(0, 8, 8);
		Hit5 *hit3 = createHit5(1, 8, 8);
		Hit5 *hit4 = createHit5(2, -1, 0);
		Hit5 *result1 = createHit5(refIdx[0], shift[0], refPos[0]);
		Hit5 *result2 = createHit5(refIdx[1], shift[1], refPos[1]);
		if (numHits == 2 && ((areHitsEqual5(result1, hit1) == 0
						&& areHitsEqual5(result2, hit2) == 0)
						|| (areHitsEqual5(result1, hit2) == 0
								&& areHitsEqual5(result2, hit1) == 0)
						|| (areHitsEqual5(result1, hit1) == 0
								&& areHitsEqual5(result2, hit3) == 0)
						|| (areHitsEqual5(result1, hit3) == 0
								&& areHitsEqual5(result2, hit1) == 0)
						|| (areHitsEqual5(result1, hit1) == 0
								&& areHitsEqual5(result2, hit4) == 0)
						|| (areHitsEqual5(result1, hit4) == 0
								&& areHitsEqual5(result2, hit1) == 0)
						|| (areHitsEqual5(result1, hit2) == 0
								&& areHitsEqual5(result2, hit3) == 0)
						|| (areHitsEqual5(result1, hit3) == 0
								&& areHitsEqual5(result2, hit2) == 0)
						|| (areHitsEqual5(result1, hit2) == 0
								&& areHitsEqual5(result2, hit4) == 0)
						|| (areHitsEqual5(result1, hit4) == 0
								&& areHitsEqual5(result2, hit2) == 0)
						|| (areHitsEqual5(result1, hit3) == 0
								&& areHitsEqual5(result2, hit4) == 0)
						|| (areHitsEqual5(result1, hit4) == 0
								&& areHitsEqual5(result2, hit3) == 0)))
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
		numHits = lookupTable5MapQry(query, seqLength, isRevComplement,
				refIdx, shift, refPos);
		hit1 = createHit5(0, -2, 0);
		hit2 = createHit5(0, 6, 8);
		hit3 = createHit5(1, 6, 8);
		result1 = createHit5(refIdx[0], shift[0], refPos[0]);
		result2 = createHit5(refIdx[1], shift[1], refPos[1]);
		if (numHits == 2 && ((areHitsEqual5(result1, hit1) == 0
						&& areHitsEqual5(result2, hit2) == 0)
						|| (areHitsEqual5(result1, hit2) == 0
								&& areHitsEqual5(result2, hit1) == 0)
						|| (areHitsEqual5(result1, hit1) == 0
								&& areHitsEqual5(result2, hit3) == 0)
						|| (areHitsEqual5(result1, hit3) == 0
								&& areHitsEqual5(result2, hit1) == 0)
						|| (areHitsEqual5(result1, hit2) == 0
								&& areHitsEqual5(result2, hit3) == 0)
						|| (areHitsEqual5(result1, hit3) == 0
								&& areHitsEqual5(result2, hit2) == 0)))
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
		lookupTable5Delete();
		remove(refFile);
		free(refIdx);
		free(shift);
		free(refPos);
		preprocessDelete();
	}
}
END_TEST


/**
 * Tests the @a getHash_gpu function.
 */
START_TEST(getHash_gpu)
{
	char str[] = "ACGT";
	int len = strlen(str);

	char *str_d;
	cudaMalloc((void **) &str_d, len * sizeof(char));
	PRINT_CUDA_ERROR()
	cudaMemcpy(str_d, str, len * sizeof(char), cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	int *hash_d;
	cudaMalloc((void **) &hash_d, sizeof(int));
	PRINT_CUDA_ERROR()

	lookupTable5CpyConstMemToGPU();

	getHash_gpu_wrap<<<1, 1>>>(str_d, len, hash_d);
	PRINT_CUDA_ERROR()

	int hash;
	cudaMemcpy(&hash, hash_d, sizeof(int), cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	if (hash != 228)
		fail("Incorrect behavior.\n");

	cudaFree(str_d);
	cudaFree(hash_d);
}
END_TEST


/**
 * Tests @a sort_gpu function.
 */
START_TEST(sort_gpu)
{
	/* Array size is 2. */
	{
		char refIdx[] = {2, 0};
		int shift[] = {-1, 5};
		int refPos[] = {5, 10};
		int arrSize = 2;

		char *refIdx_d;
		cudaMalloc(&refIdx_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdx_d, refIdx, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shift_d;
		cudaMalloc(&shift_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shift_d, shift, arrSize * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *refPos_d;
		cudaMalloc(&refPos_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos_d, refPos, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Using 1 thread. */
		sort_gpu_wrap<<<1, 1>>>(refIdx_d, shift_d, refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdx, refIdx_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shift, shift_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx[0] != 0 || shift[0] != 5 || refPos[0] != 10
				|| refIdx[1] != 2 || shift[1] != -1 || refPos[1] != 5)
			fail("Incorrect behavior when the array size is 2 and number of "
					"threads launched is 1.\n");

		/* Using 2 threads. */
		sort_gpu_wrap<<<1, 2>>>(refIdx_d, shift_d, refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdx, refIdx_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shift, shift_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx[0] != 0 || shift[0] != 5 || refPos[0] != 10
				|| refIdx[1] != 2 || shift[1] != -1 || refPos[1] != 5)
			fail("Incorrect behavior when the array size is 2 and number of "
					"threads launched is 2.\n");

		/* Using 3 threads. */
		sort_gpu_wrap<<<1, 3>>>(refIdx_d, shift_d, refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdx, refIdx_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shift, shift_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx[0] != 0 || shift[0] != 5 || refPos[0] != 10
				|| refIdx[1] != 2 || shift[1] != -1 || refPos[1] != 5)
			fail("Incorrect behavior when the array size is 2 and number of "
					"threads launched is 3.\n");

		cudaFree(refIdx_d);
		cudaFree(shift_d);
		cudaFree(refPos_d);
	}

	/* Array size is 5. */
	{
		char refIdx[] = {2, 0, 3, 3, 1};
		int shift[] = {-1, 5, -5, -5, 7};
		int refPos[] = {5, 10, 18, 13, 0};
		int arrSize = 5;

		char *refIdx_d;
		cudaMalloc(&refIdx_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdx_d, refIdx, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shift_d;
		cudaMalloc(&shift_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shift_d, shift, arrSize * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *refPos_d;
		cudaMalloc(&refPos_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos_d, refPos, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Using 1 thread. */
		sort_gpu_wrap<<<1, 1>>>(refIdx_d, shift_d, refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdx, refIdx_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shift, shift_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx[0] != 0 || shift[0] != 5 || refPos[0] != 10
				|| refIdx[1] != 1 || shift[1] != 7 || refPos[1] != 0
				|| refIdx[2] != 2 || shift[2] != -1 || refPos[2] != 5
				|| refIdx[3] != 3 || shift[3] != -5 || refPos[3] != 13
				|| refIdx[4] != 3 || shift[4] != -5 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 1.\n");

		/* Using 2 threads. */
		sort_gpu_wrap<<<1, 2>>>(refIdx_d, shift_d, refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdx, refIdx_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shift, shift_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx[0] != 0 || shift[0] != 5 || refPos[0] != 10
				|| refIdx[1] != 1 || shift[1] != 7 || refPos[1] != 0
				|| refIdx[2] != 2 || shift[2] != -1 || refPos[2] != 5
				|| refIdx[3] != 3 || shift[3] != -5 || refPos[3] != 13
				|| refIdx[4] != 3 || shift[4] != -5 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 2.\n");

		/* Using 3 threads. */
		sort_gpu_wrap<<<1, 3>>>(refIdx_d, shift_d, refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdx, refIdx_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shift, shift_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx[0] != 0 || shift[0] != 5 || refPos[0] != 10
				|| refIdx[1] != 1 || shift[1] != 7 || refPos[1] != 0
				|| refIdx[2] != 2 || shift[2] != -1 || refPos[2] != 5
				|| refIdx[3] != 3 || shift[3] != -5 || refPos[3] != 13
				|| refIdx[4] != 3 || shift[4] != -5 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 3.\n");

		/* Using 4 threads. */
		sort_gpu_wrap<<<1, 4>>>(refIdx_d, shift_d, refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdx, refIdx_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shift, shift_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx[0] != 0 || shift[0] != 5 || refPos[0] != 10
				|| refIdx[1] != 1 || shift[1] != 7 || refPos[1] != 0
				|| refIdx[2] != 2 || shift[2] != -1 || refPos[2] != 5
				|| refIdx[3] != 3 || shift[3] != -5 || refPos[3] != 13
				|| refIdx[4] != 3 || shift[4] != -5 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 4.\n");

		/* Using 5 threads. */
		sort_gpu_wrap<<<1, 5>>>(refIdx_d, shift_d, refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdx, refIdx_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shift, shift_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx[0] != 0 || shift[0] != 5 || refPos[0] != 10
				|| refIdx[1] != 1 || shift[1] != 7 || refPos[1] != 0
				|| refIdx[2] != 2 || shift[2] != -1 || refPos[2] != 5
				|| refIdx[3] != 3 || shift[3] != -5 || refPos[3] != 13
				|| refIdx[4] != 3 || shift[4] != -5 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 5.\n");

		/* Using 6 threads. */
		sort_gpu_wrap<<<1, 6>>>(refIdx_d, shift_d, refPos_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdx, refIdx_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shift, shift_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refPos, refPos_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx[0] != 0 || shift[0] != 5 || refPos[0] != 10
				|| refIdx[1] != 1 || shift[1] != 7 || refPos[1] != 0
				|| refIdx[2] != 2 || shift[2] != -1 || refPos[2] != 5
				|| refIdx[3] != 3 || shift[3] != -5 || refPos[3] != 13
				|| refIdx[4] != 3 || shift[4] != -5 || refPos[4] != 18)
			fail("Incorrect behavior when the array size is 5 and number of "
					"threads launched is 6.\n");

		cudaFree(refIdx_d);
		cudaFree(shift_d);
		cudaFree(refPos_d);
	}
}
END_TEST


/**
 * Tests @a createClusters_gpu function.
 */
START_TEST(createClusters_gpu)
{
	/* Array size is 5, there are 4 clusters, and 1 biggest cluster. */
	{
		char refIdx[] = {0, 1, 3, 3, 9};
		int shift[] = {-1, 5, -5, -5, 7};
		int arrSize = 5;

		char *refIdx_d;
		cudaMalloc(&refIdx_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdx_d, refIdx, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shift_d;
		cudaMalloc(&shift_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shift_d, shift, arrSize * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *clusterSize_d;
		cudaMalloc(&clusterSize_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		short *biggestClusterSize_d;
		cudaMalloc(&biggestClusterSize_d, sizeof(short));
		PRINT_CUDA_ERROR()

		short *numClusters_d;
		cudaMalloc(&numClusters_d, sizeof(short));
		PRINT_CUDA_ERROR()

		createClusters_gpu_wrap<<<1, 1>>>(refIdx_d, shift_d, clusterSize_d,
				arrSize, biggestClusterSize_d, numClusters_d);

		int clusterSize[arrSize];
		cudaMemcpy(clusterSize, clusterSize_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short biggestClusterSize;
		cudaMemcpy(&biggestClusterSize, biggestClusterSize_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short numClusters;
		cudaMemcpy(&numClusters, numClusters_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (numClusters != 4 || biggestClusterSize != 2 || clusterSize[0] != 1
				|| clusterSize[1] != 1 || clusterSize[2] != 2
				|| clusterSize[3] != 1)
			fail("Incorrect behavior when array size is 5, there are 4 clusters, "
					" and 1 biggest cluster.\n");

		cudaFree(refIdx_d);
		cudaFree(shift_d);
		cudaFree(clusterSize_d);
		cudaFree(biggestClusterSize_d);
		cudaFree(numClusters_d);
	}

	/* Array size is 5, there is 1 cluster and 1 biggest cluster. */
	{
		char refIdx[] = {1, 1, 1, 1, 1};
		int shift[] = {-5, -5, -5, -5, -5};
		int arrSize = 5;

		char *refIdx_d;
		cudaMalloc(&refIdx_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdx_d, refIdx, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shift_d;
		cudaMalloc(&shift_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shift_d, shift, arrSize * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *clusterSize_d;
		cudaMalloc(&clusterSize_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		short *biggestClusterSize_d;
		cudaMalloc(&biggestClusterSize_d, sizeof(short));
		PRINT_CUDA_ERROR()

		short *numClusters_d;
		cudaMalloc(&numClusters_d, sizeof(short));
		PRINT_CUDA_ERROR()

		createClusters_gpu_wrap<<<1, 1>>>(refIdx_d, shift_d, clusterSize_d,
				arrSize, biggestClusterSize_d, numClusters_d);

		int clusterSize[arrSize];
		cudaMemcpy(clusterSize, clusterSize_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short biggestClusterSize;
		cudaMemcpy(&biggestClusterSize, biggestClusterSize_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short numClusters;
		cudaMemcpy(&numClusters, numClusters_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (numClusters != 1 || biggestClusterSize != 5 || clusterSize[0] != 5)
			fail("Incorrect behavior when array size is 5, there is 1 cluster "
					"and 1 biggest cluster.\n");

		cudaFree(refIdx_d);
		cudaFree(shift_d);
		cudaFree(clusterSize_d);
		cudaFree(biggestClusterSize_d);
		cudaFree(numClusters_d);
	}
}
END_TEST


/**
 * Tests @a findBiggestClusters_gpu function.
 */
START_TEST(findBiggestClusters_gpu)
{
	/* Array size is 5, there are 4 clusters, and 1 biggest cluster. */
	{
		char refIdx[] = {0, 1, 3, 3, 9};
		int shift[] = {-1, 5, -5, -5, 7};
		int arrSize = 5;

		char *refIdx_d;
		cudaMalloc(&refIdx_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdx_d, refIdx, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shift_d;
		cudaMalloc(&shift_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shift_d, shift, arrSize * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *clusterSize_d;
		cudaMalloc(&clusterSize_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		short *biggestClusterSize_d;
		cudaMalloc(&biggestClusterSize_d, sizeof(short));
		PRINT_CUDA_ERROR()

		short *numClusters_d;
		cudaMalloc(&numClusters_d, sizeof(short));
		PRINT_CUDA_ERROR()

		createClusters_gpu_wrap<<<1, 1>>>(refIdx_d, shift_d, clusterSize_d,
				arrSize, biggestClusterSize_d, numClusters_d);

		int clusterSize[arrSize];
		cudaMemcpy(clusterSize, clusterSize_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short biggestClusterSize;
		cudaMemcpy(&biggestClusterSize, biggestClusterSize_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short numClusters;
		cudaMemcpy(&numClusters, numClusters_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *biggestClusters_d;
		cudaMalloc(&biggestClusters_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		char *numBiggestClusters_d;
		cudaMalloc(&numBiggestClusters_d, sizeof(char));
		PRINT_CUDA_ERROR()

		findBiggestClusters_gpu_wrap<<<1, 1>>>(numClusters, biggestClusterSize,
				clusterSize_d, biggestClusters_d, numBiggestClusters_d);
		PRINT_CUDA_ERROR()

		int biggestClusters[arrSize];
		cudaMemcpy(biggestClusters, biggestClusters_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		char numBiggestClusters;
		cudaMemcpy(&numBiggestClusters, numBiggestClusters_d,
				sizeof(char), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (numBiggestClusters != 1 || biggestClusters[0] != 2)
			fail("Incorrect behavior when array size is 5, there are 4 clusters, "
					" and 1 biggest cluster.\n");

		cudaFree(refIdx_d);
		cudaFree(shift_d);
		cudaFree(clusterSize_d);
		cudaFree(biggestClusterSize_d);
		cudaFree(numClusters_d);
		cudaFree(biggestClusters_d);
		cudaFree(numBiggestClusters_d);
	}

	/* Array size is 5, there are 3 clusters, and 2 biggest clusters. */
	{
		char refIdx[] = {0, 0, 3, 3, 9};
		int shift[] = {5, 5, -5, -5, 7};
		int arrSize = 5;

		char *refIdx_d;
		cudaMalloc(&refIdx_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdx_d, refIdx, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shift_d;
		cudaMalloc(&shift_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shift_d, shift, arrSize * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *clusterSize_d;
		cudaMalloc(&clusterSize_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		short *biggestClusterSize_d;
		cudaMalloc(&biggestClusterSize_d, sizeof(short));
		PRINT_CUDA_ERROR()

		short *numClusters_d;
		cudaMalloc(&numClusters_d, sizeof(short));
		PRINT_CUDA_ERROR()

		createClusters_gpu_wrap<<<1, 1>>>(refIdx_d, shift_d, clusterSize_d,
				arrSize, biggestClusterSize_d, numClusters_d);

		int clusterSize[arrSize];
		cudaMemcpy(clusterSize, clusterSize_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short biggestClusterSize;
		cudaMemcpy(&biggestClusterSize, biggestClusterSize_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		short numClusters;
		cudaMemcpy(&numClusters, numClusters_d, sizeof(short),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *biggestClusters_d;
		cudaMalloc(&biggestClusters_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		char *numBiggestClusters_d;
		cudaMalloc(&numBiggestClusters_d, sizeof(char));
		PRINT_CUDA_ERROR()

		findBiggestClusters_gpu_wrap<<<1, 1>>>(numClusters, biggestClusterSize,
				clusterSize_d, biggestClusters_d, numBiggestClusters_d);
		PRINT_CUDA_ERROR()

		int biggestClusters[arrSize];
		cudaMemcpy(biggestClusters, biggestClusters_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		char numBiggestClusters;
		cudaMemcpy(&numBiggestClusters, numBiggestClusters_d,
				sizeof(char), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (numBiggestClusters != 2 || biggestClusters[0] != 0
				|| biggestClusters[1] != 2)
			fail("Incorrect behavior when array size is 5, there are 3 clusters, "
					" and 2 biggest clusters.\n");

		cudaFree(refIdx_d);
		cudaFree(shift_d);
		cudaFree(clusterSize_d);
		cudaFree(biggestClusterSize_d);
		cudaFree(numClusters_d);
		cudaFree(biggestClusters_d);
		cudaFree(numBiggestClusters_d);
	}
}
END_TEST


/**
 * Tests @a assignResults_gpu function.
 */
START_TEST(assignResults_gpu)
{
	/* Num hits is 1 and max allowed hits is 1. */
	{
		char maxHits = 1;
		char numBgstHits = 1;

		char refIdx[] = {0, 1, 3, 3, 9};
		int refPos[] = {1, 2, 3, 4, 5};
		int arrSize = 5;

		char *refIdx1_d;
		cudaMalloc(&refIdx1_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdx1_d, refIdx, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		char *refIdx2_d;
		cudaMalloc(&refIdx2_d, sizeof(char));
		PRINT_CUDA_ERROR()

		int *refPos1_d;
		cudaMalloc(&refPos1_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos1_d, refPos, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *refPos2_d;
		cudaMalloc(&refPos2_d, sizeof(int));
		PRINT_CUDA_ERROR()

		int *bgstClust_d;
		cudaMalloc(&bgstClust_d, sizeof(int));
		PRINT_CUDA_ERROR()

		int bgstClust = 2;
		cudaMemcpy(bgstClust_d, &bgstClust, numBgstHits * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		srand(time(NULL));
		int randNum = rand();
		assignResults_gpu_wrap<<<1, 1>>>(maxHits, numBgstHits, refIdx1_d,
				refIdx2_d, refPos1_d, refPos2_d, bgstClust_d, randNum);
		PRINT_CUDA_ERROR()

		char refIdx2;
		cudaMemcpy(&refIdx2, refIdx2_d, numBgstHits * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int refPos2;
		cudaMemcpy(&refPos2, refPos2_d, numBgstHits * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx2 != 3 || refPos2 != 3)
			fail("Incorrect behavior when number of hits is 1 and max allowed "
					"hits is 1.\n");

		cudaFree(refIdx1_d);
		cudaFree(refIdx2_d);
		cudaFree(refPos1_d);
		cudaFree(refPos2_d);
		cudaFree(bgstClust_d);
	}

	/* Num hits is 2 and max allowed hits is 1. */
	{
		char maxHits = 1;
		char numBgstHits = 2;

		char refIdx[] = {0, 0, 3, 3, 9};
		int refPos[] = {1, 2, 3, 4, 5};
		int arrSize = 5;

		char *refIdx1_d;
		cudaMalloc(&refIdx1_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdx1_d, refIdx, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		char *refIdx2_d;
		cudaMalloc(&refIdx2_d, maxHits * sizeof(char));
		PRINT_CUDA_ERROR()

		int *refPos1_d;
		cudaMalloc(&refPos1_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos1_d, refPos, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *refPos2_d;
		cudaMalloc(&refPos2_d, maxHits * sizeof(int));
		PRINT_CUDA_ERROR()

		int *bgstClust_d;
		cudaMalloc(&bgstClust_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		int bgstClust[] = {0, 2};
		cudaMemcpy(bgstClust_d, &bgstClust, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		srand(time(NULL));
		int randNum = rand();
		assignResults_gpu_wrap<<<1, 1>>>(maxHits, numBgstHits, refIdx1_d,
				refIdx2_d, refPos1_d, refPos2_d, bgstClust_d, randNum);
		PRINT_CUDA_ERROR()

		char refIdx2;
		cudaMemcpy(&refIdx2, refIdx2_d, maxHits * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int refPos2;
		cudaMemcpy(&refPos2, refPos2_d, maxHits * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if ((refIdx2 != 0 && refIdx2 != 3) || (refPos2 != 1 && refPos2 != 3))
			fail("Incorrect behavior when number of hits is 2 and max allowed "
					"hits is 1.\n");

		cudaFree(refIdx1_d);
		cudaFree(refIdx2_d);
		cudaFree(refPos1_d);
		cudaFree(refPos2_d);
		cudaFree(bgstClust_d);
	}

	/* Num hits is 1 and max allowed hits is 2. */
	{
		char maxHits = 2;
		char numBgstHits = 1;

		char refIdx[] = {0, 1, 3, 3, 9};
		int refPos[] = {1, 2, 3, 4, 5};
		int arrSize = 5;

		char *refIdx1_d;
		cudaMalloc(&refIdx1_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdx1_d, refIdx, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		char *refIdx2_d;
		cudaMalloc(&refIdx2_d, maxHits * sizeof(char));
		PRINT_CUDA_ERROR()

		int *refPos1_d;
		cudaMalloc(&refPos1_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refPos1_d, refPos, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *refPos2_d;
		cudaMalloc(&refPos2_d, maxHits * sizeof(int));
		PRINT_CUDA_ERROR()

		int *bgstClust_d;
		cudaMalloc(&bgstClust_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()

		int bgstClust[] = {2};
		cudaMemcpy(bgstClust_d, &bgstClust, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		srand(time(NULL));
		int randNum = rand();
		assignResults_gpu_wrap<<<1, 1>>>(maxHits, numBgstHits, refIdx1_d,
				refIdx2_d, refPos1_d, refPos2_d, bgstClust_d, randNum);
		PRINT_CUDA_ERROR()

		char refIdx2;
		cudaMemcpy(&refIdx2, refIdx2_d, maxHits * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int refPos2;
		cudaMemcpy(&refPos2, refPos2_d, maxHits * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx2 != 3 || refPos2 != 3)
			fail("Incorrect behavior when number of hits is 1 and max allowed "
					"hits is 2.\n");

		cudaFree(refIdx1_d);
		cudaFree(refIdx2_d);
		cudaFree(refPos1_d);
		cudaFree(refPos2_d);
		cudaFree(bgstClust_d);
	}
}
END_TEST


/**
 * Tests @a lookupTable5MapQry_gpu function.
 */
START_TEST(lookupTable5MapQry_gpu)
{
	/* Case 1. */
	{
		int seedLen = 5;
		int numDistinctTuples = (int) pow((float) DNA_ALPHABET_SIZE,
				(int) seedLen);
		int *keys = (int *) calloc(numDistinctTuples, sizeof(int));
		keys[228] = 1;
		keys[377] = 3;
		keys[916] = 4;

		int numTotalTuples = 4;
		int *vals = (int *) calloc((numTotalTuples + 1), sizeof(int));
		vals[0] = -1;
		vals[1] = 1;
		vals[2] = 3;
		vals[3] = 2;
		vals[4] = 0;

		int *numRptsPerTuple = (int *) calloc(numDistinctTuples, sizeof(int));
		numRptsPerTuple[228] = 2;
		numRptsPerTuple[377] = 1;
		numRptsPerTuple[916] = 1;

		int *keys_d;
		cudaMalloc(&keys_d, numDistinctTuples * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(keys_d, keys, numDistinctTuples * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *vals_d;
		cudaMalloc(&vals_d, (numTotalTuples + 1) * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(vals_d, vals, (numTotalTuples + 1) * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *numRptsPerTuple_d;
		cudaMalloc(&numRptsPerTuple_d, numDistinctTuples * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(numRptsPerTuple_d, numRptsPerTuple,
				numDistinctTuples * sizeof(int), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int maxQrySeqLen = 12;
		int numQrs = 2;
		char *qrs = (char *) calloc(maxQrySeqLen * numQrs, sizeof(char));
		strcpy(qrs, "ACGTACGTCC");
		strcpy(qrs + maxQrySeqLen, "GGACGTACGT");

		uchar *qryLen = (uchar *) calloc(numQrs, sizeof(uchar));
		qryLen[0] = 10;
		qryLen[1] = 10;

		char *qrs_d;
		cudaMalloc(&qrs_d, maxQrySeqLen * numQrs * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(qrs_d, qrs, maxQrySeqLen * numQrs * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		uchar *qryLen_d;
		cudaMalloc(&qryLen_d, numQrs * sizeof(uchar));
		PRINT_CUDA_ERROR()
		cudaMemcpy(qryLen_d, qryLen, numQrs * sizeof(uchar),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		char *refIdx_d;
		cudaMalloc(&refIdx_d, numQrs * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemset(refIdx_d, CHAR_MAX, numQrs * sizeof(char));
		PRINT_CUDA_ERROR()

		int *refPos_d;
		cudaMalloc(&refPos_d, numQrs * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemset(refPos_d, -1, numQrs * sizeof(int));
		PRINT_CUDA_ERROR()

		lookupTable5CpyConstMemToGPU();

		int maxHits = 1;
		srand(time(NULL));
		int randNum = rand();
		int arrSize = 100;
		int shrMem = (4 * arrSize * sizeof(int)) + (maxHits * sizeof(int));
		lookupTable5MapQry2_gpu<<<2, 1, shrMem>>>(keys_d, vals_d,
				numRptsPerTuple_d, qrs_d, qryLen_d, maxQrySeqLen, refIdx_d,
				refPos_d, maxHits, seedLen, randNum, arrSize);
		PRINT_CUDA_ERROR()

		char *refIdx = (char *) malloc(numQrs * sizeof(char));
		cudaMemcpy(refIdx, refIdx_d, numQrs * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *refPos = (int *) malloc(numQrs * sizeof(int));
		cudaMemcpy(refPos, refPos_d, numQrs * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdx[0] != 0 || refPos[0] != 5 || refIdx[1] != 0
				|| (refPos[1] != 5 && refPos[1] != 15))
			fail("Incorrect behavior.\n");

		cudaFree(keys_d);
		cudaFree(vals_d);
		cudaFree(numRptsPerTuple_d);
		cudaFree(qrs_d);
		cudaFree(qryLen_d);
		cudaFree(refIdx_d);
		cudaFree(refPos_d);
		free(keys);
		free(vals);
		free(numRptsPerTuple);
		free(qrs);
		free(qryLen);
		free(refIdx);
		free(refPos);
	}
}
END_TEST


/**
 * Tests @a initializeShrMem_gpu function.
 */
START_TEST(initializeShrMem_gpu)
{
	int arrSize = 5;
	char *refIdx = (char *) calloc(arrSize, sizeof(char));
	char *refIdx_d;
	cudaMalloc(&refIdx_d, arrSize * sizeof(char));
	PRINT_CUDA_ERROR()

	int *shift = (int *) calloc(arrSize, sizeof(int));
	int *shift_d;
	cudaMalloc(&shift_d, arrSize * sizeof(int));
	PRINT_CUDA_ERROR()

	int *refPos = (int *) calloc(arrSize, sizeof(int));
	int *refPos_d;
	cudaMalloc(&refPos_d, arrSize * sizeof(int));
	PRINT_CUDA_ERROR()

	int *clusterSize = (int *) calloc(arrSize, sizeof(int));
	int *clusterSize_d;
	cudaMalloc(&clusterSize_d, arrSize * sizeof(int));
	PRINT_CUDA_ERROR()

	intializeShrMem_gpu_wrap<<<1, 5>>>(refIdx_d, shift_d, refPos_d,
			clusterSize_d, arrSize);

	cudaMemcpy(refIdx, refIdx_d, arrSize * sizeof(char),
			cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	cudaMemcpy(shift, shift_d, arrSize * sizeof(int),
			cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	cudaMemcpy(refPos, refPos_d, arrSize * sizeof(int),
			cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	cudaMemcpy(clusterSize, clusterSize_d, arrSize * sizeof(int),
			cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	int i;
	for (i = 0; i < arrSize; ++i)
	{
		if (refIdx[i] != CHAR_MAX)
			fail("Incorrect behavior: Reference index array is not initialized "
					"properly.\n");
		if (shift[i] != -1 || refPos[i] != -1 || clusterSize[i] != -1)
			fail("Incorrect behavior: Shift array not initialized "
					"properly.\n");
		if (refPos[i] != -1)
			fail("Incorrect behavior: Reference position array is not "
					"initialized properly.\n");
		if (clusterSize[i] != -1)
			fail("Incorrect behavior: Cluster size array is not "
					"initialized properly.\n");
	}

	cudaFree(refIdx_d);
	cudaFree(shift_d);
	cudaFree(refPos_d);
	cudaFree(clusterSize_d);
	free(refIdx);
	free(shift);
	free(refPos);
	free(clusterSize);
}
END_TEST


/**
 * Tests @a cpyHitsFromGlobalToShr_gpu function.
 */
START_TEST(cpyHitsFromGlobalToShr_gpu)
{
	int seedLen = 5;
	int numDistinctTuples = (int) pow((float) DNA_ALPHABET_SIZE,
			(int) seedLen);
	int *keys = (int *) calloc(numDistinctTuples, sizeof(int));
	keys[228] = 1;
	keys[377] = 3;
	keys[916] = 4;

	int numTotalTuples = 4;
	int *vals = (int *) calloc((numTotalTuples + 1), sizeof(int));
	vals[0] = -1;
	vals[1] = 1;
	vals[2] = 33554435;
	vals[3] = 2;
	vals[4] = 0;

	int *numRptsPerTuple = (int *) calloc(numDistinctTuples, sizeof(int));
	numRptsPerTuple[228] = 2;
	numRptsPerTuple[377] = 1;
	numRptsPerTuple[916] = 1;

	int *keys_d;
	cudaMalloc(&keys_d, numDistinctTuples * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(keys_d, keys, numDistinctTuples * sizeof(int),
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	int *vals_d;
	cudaMalloc(&vals_d, (numTotalTuples + 1) * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(vals_d, vals, (numTotalTuples + 1) * sizeof(int),
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	int *numRptsPerTuple_d;
	cudaMalloc(&numRptsPerTuple_d, numDistinctTuples * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(numRptsPerTuple_d, numRptsPerTuple,
			numDistinctTuples * sizeof(int), cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	int arrSize = 100;
	char *refIdx_d;
	cudaMalloc(&refIdx_d, arrSize * sizeof(char));
	PRINT_CUDA_ERROR()

	int *shift_d;
	cudaMalloc(&shift_d, arrSize * sizeof(int));
	PRINT_CUDA_ERROR()

	int *refPos_d;
	cudaMalloc(&refPos_d, arrSize * sizeof(int));
	PRINT_CUDA_ERROR()

	int numQryTuples = 6;
	int *hashes = (int *) calloc(numQryTuples, sizeof(int));
	hashes[0] = 228;
	hashes[1] = 313;
	hashes[2] = 590;
	hashes[3] = 915;
	hashes[4] = 484;
	hashes[5] = 377;

	int *hashes_d;
	cudaMalloc(&hashes_d, numQryTuples * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(hashes_d, hashes, numQryTuples * sizeof(int),
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	cpyHitsFromGlobalToShr_gpu_wrap<<<1, 1>>>(refIdx_d, shift_d, refPos_d,
			hashes_d, arrSize, keys_d, vals_d, numRptsPerTuple_d,
			numQryTuples, seedLen);
	PRINT_CUDA_ERROR()

	char *refIdx = (char *) calloc(arrSize, sizeof(char));
	cudaMemcpy(refIdx, refIdx_d, arrSize * sizeof(char), cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	int *shift = (int *) calloc(arrSize, sizeof(int));
	cudaMemcpy(shift, shift_d, arrSize * sizeof(int), cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	int *refPos = (int *) calloc(arrSize, sizeof(int));
	cudaMemcpy(refPos, refPos_d, arrSize * sizeof(int), cudaMemcpyDeviceToHost);
	PRINT_CUDA_ERROR()

	if (refIdx[0] != 0 || shift[0] != 5 || refPos[0] != 5
			|| refIdx[1] != 1 || shift[1] != 15 || refPos[1] != 15
			|| refIdx[2] != 0 || shift[2] != 5 || refPos[2] != 10)
		fail("Incorrect behavior.\n");

	cudaFree(keys_d);
	cudaFree(vals_d);
	cudaFree(numRptsPerTuple_d);
	cudaFree(refIdx_d);
	cudaFree(shift_d);
	cudaFree(refPos_d);
	cudaFree(hashes_d);
	free(keys);
	free(vals);
	free(numRptsPerTuple);
	free(refIdx);
	free(shift);
	free(refPos);
	free(hashes);
}
END_TEST


/**
 * Creates test suite.
 */
Suite *lookupTable5Suite(void)
{
	Suite *s = suite_create("lookupTable5");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, lookupTable5Create);
	tcase_add_test(testCaseCore, lookupTable5MapQry);
	tcase_add_test(testCaseCore, getHash_gpu);
	tcase_add_test(testCaseCore, sort_gpu);
	tcase_add_test(testCaseCore, createClusters_gpu);
	tcase_add_test(testCaseCore, findBiggestClusters_gpu);
	tcase_add_test(testCaseCore, assignResults_gpu);
	tcase_add_test(testCaseCore, initializeShrMem_gpu);
	tcase_add_test(testCaseCore, cpyHitsFromGlobalToShr_gpu);
//	tcase_add_test(testCaseCore, lookupTable5MapQry_gpu);
	suite_add_tcase (s, testCaseCore);

	return s;
}


/**
 * Returns a 'Hit' data structure.
 *
 * @param	refIdx	Reference index
 * @param	shift		Shift.
 * @param	refPos	Reference position.
 * @return	Creates a 'hit' data structure.
 */
static Hit5 *createHit5(char refIdx, int shift, int refPos)
{
	Hit5 *hit = (Hit5 *) malloc(sizeof(Hit5));
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
static int areHitsEqual5(const Hit5 *hit1, const Hit5 *hit2)
{
	if (hit1->refIdx == hit2->refIdx && hit1->shift == hit2->shift
			&& hit1->refPos == hit2->refPos)
		return 0;
	else
		return 1;
}

