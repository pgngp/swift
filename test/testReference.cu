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
#include "testReference.h"
#include "../src/reference.h"
#include "../src/preprocess.h"
#include <signal.h>
#include <limits.h>


/**
 * Handles OS signals.
 */
static void signalHandler(int param)
{
	if (param == SIGSEGV)
		fail("Received segmentation fault.\n");
}


/**
 * Tests refSearchQuery function.
 */
START_TEST(searchQuery)
{
	/* Test for correct results with a normal query (case 1). */
	{
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp1.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);
		int seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		uint maxNumHits = 10;
		HitList *result = refSearchQuery(hashTable, query, seqLength,
				isRevComplement, seedLength, maxNumHits);
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
		if (hashTable != NULL)
		{
			int i;
			int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE, seedLength);
			for (i = 0; i < numRefTuples; ++i)
				refPosListDelete(hashTable[i]);
			free(hashTable);
		}
		remove(refFile);
	}

	/* Test for correct results with a normal query (case 2). */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp2.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAACGTACGTATTTACGTA\n"
				">ref3\nCGTACGTAAAATAACGTCGTACGTAAACGTACGTACGTAA", filePtr);
		fclose(filePtr);
		int seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);
		Query *query = qryCreate("qry1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 0;
		uint maxNumHits = 10;
		HitList *result = refSearchQuery(hashTable, query, seqLength,
				isRevComplement, seedLength, maxNumHits);
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
		if (hashTable != NULL)
		{
			int i;
			int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE, seedLength);
			for (i = 0; i < numRefTuples; ++i)
				refPosListDelete(hashTable[i]);
			free(hashTable);
		}
		remove(refFile);
	}

	/* Test for correct results with the reverse complement of a query. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp2.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATAACGTACGTATTTACGTACACCCGTACGT\n>ref2\n"
				"GCCGTACGACGTACGTAAACGTAACGTACGTATTTACGTA\n"
				">ref3\nCGTACGTAAAATAACGTCGTACGTAAACGTACGTACGTAA", filePtr);
		fclose(filePtr);
		int seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);
		Query *query = qryCreate("qry1", "TTACGTACGT", "ACGTACGTAA", 0, 5);
		uint seqLength = 10;
		uint isRevComplement = 1;
		uint maxNumHits = 10;
		HitList *result = refSearchQuery(hashTable, query, seqLength,
				isRevComplement, seedLength, maxNumHits);
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
		if (hashTable != NULL)
		{
			int i;
			int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE, seedLength);
			for (i = 0; i < numRefTuples; ++i)
				refPosListDelete(hashTable[i]);
			free(hashTable);
		}
		remove(refFile);
	}
}
END_TEST


/**
 * Tests @a refSearchQuery_paired function.
 */
START_TEST(searchQuery_paired)
{
	/* Test behavior when hash table is NULL. */
	{
		Query *qry1 = qryCreate("q1_1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		Query *qry2 = qryCreate("q1_2", "TTACAGTCGA", "TCGACTGTAA", 0, 5);
		uint seqLen = 0;
		RefPosList **hashTable = NULL;
		int seedLen = 8;
		uint maxNumHits = 10, minFragSize = 8, maxFragSize = 20;
		HitList **result = (HitList **) malloc(8 * sizeof(HitList *));
		result[0] = NULL;
		result[1] = NULL;
		result[2] = NULL;
		result[3] = NULL;
		result[4] = NULL;
		result[5] = NULL;
		result[6] = NULL;
		result[7] = NULL;
		refSearchQuery_paired(hashTable, qry1, qry2, seqLen, seedLen,
					maxNumHits, minFragSize, maxFragSize, result);
		if (result[0] != NULL || result[1] != NULL || result[2] != NULL
				|| result[3] != NULL || result[4] != NULL || result[5] != NULL
				|| result[6] != NULL || result[7] != NULL)
			fail("Incorrect behavior when hash table is NULL.\n");
		qryDelete(qry1);
		qryDelete(qry2);
		free(result);
	}

	/* Test behavior when query object is NULL. */
	{
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp1.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);
		int seedLen = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLen);
		Query *qry1 = NULL;
		Query *qry2 = NULL;
		uint seqLen = 100;
		uint maxNumHits = 10, minFragSize = 8, maxFragSize = 20;
		HitList **result = (HitList **) malloc(8 * sizeof(HitList *));
		result[0] = NULL;
		result[1] = NULL;
		result[2] = NULL;
		result[3] = NULL;
		result[4] = NULL;
		result[5] = NULL;
		result[6] = NULL;
		result[7] = NULL;
		refSearchQuery_paired(hashTable, qry1, qry2, seqLen, seedLen,
					maxNumHits, minFragSize, maxFragSize, result);
		if (result[0] != NULL || result[1] != NULL || result[2] != NULL
				|| result[3] != NULL || result[4] != NULL || result[5] != NULL
				|| result[6] != NULL || result[7] != NULL)
			fail("Incorrect behavior when query object is NULL.\n");
		if (hashTable != NULL)
		{
			int i;
			int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE, seedLen);
			for (i = 0; i < numRefTuples; ++i)
				refPosListDelete(hashTable[i]);
			free(hashTable);
		}
		remove(refFile);
		free(result);
	}

	/* Test behavior when query sequence length is 0. */
	{
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp1.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTGTACGTAACATC\n", filePtr);
		fclose(filePtr);
		int seedLen = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLen);
		Query *qry1 = qryCreate("q1_1", "ACGTACGTAA", "TTACGTACGT", 0, 5);
		Query *qry2 = qryCreate("q1_2", "TTACAGTCGA", "TCGACTGTAA", 0, 5);
		uint seqLen = 0;
		uint maxNumHits = 10, minFragSize = 8, maxFragSize = 20;
		HitList **result = (HitList **) malloc(8 * sizeof(HitList *));
		result[0] = NULL;
		result[1] = NULL;
		result[2] = NULL;
		result[3] = NULL;
		result[4] = NULL;
		result[5] = NULL;
		result[6] = NULL;
		result[7] = NULL;
		refSearchQuery_paired(hashTable, qry1, qry2, seqLen, seedLen,
					maxNumHits, minFragSize, maxFragSize, result);
		if (result[0] != NULL || result[1] != NULL || result[2] != NULL
				|| result[3] != NULL || result[4] != NULL || result[5] != NULL
				|| result[6] != NULL || result[7] != NULL)
			fail("Incorrect behavior when sequence length is 0.\n");
		qryDelete(qry1);
		qryDelete(qry2);
		if (hashTable != NULL)
		{
			int i;
			int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE, seedLen);
			for (i = 0; i < numRefTuples; ++i)
				refPosListDelete(hashTable[i]);
			free(hashTable);
		}
		remove(refFile);
		free(result);
	}

	/* Test for correct results with a normal pair of queries (case 1). */
	{
		char refFile[MAX_FILE_NAME_LENGTH] = "\0";
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "ref_tmp1.fa");
		FILE *filePtr = fopen(refFile, "w");
		/* ACGTACGT GTACGTAA ACAGTCGA CATC */
		fputs(">ref1\nACGTACGTGTACGTAAACAGTCGACATC\n", filePtr);
		fclose(filePtr);
		int seedLen = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLen);
		Query *qry1 = qryCreate("q1_1", "ACGTACGTAA", "TTACGTACGT", 0, 6);
		Query *qry2 = qryCreate("q1_2", "TTACAGTCGA", "TCGACTGTAA", 0, 6);
		uint seqLen = 10;
		uint maxNumHits = 10, minFragSize = 8, maxFragSize = 30;
		HitList **result = (HitList **) malloc(8 * sizeof(HitList *));
		result[0] = NULL;
		result[1] = NULL;
		result[2] = NULL;
		result[3] = NULL;
		result[4] = NULL;
		result[5] = NULL;
		result[6] = NULL;
		result[7] = NULL;
		refSearchQuery_paired(hashTable, qry1, qry2, seqLen, seedLen,
					maxNumHits, minFragSize, maxFragSize, result);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 6, 8);
		HitList *node3 = hitListCreateNode(0, 14, 16);
		HitList *node4 = hitListCreateNode(0, 14, 16);
		HitList *node5 = hitListCreateNode(0, -2, 0);
		HitList *node6 = hitListCreateNode(0, 14, 16);
		if (result != NULL && result[0] != NULL && result[1] != NULL
				&& result[2] != NULL && result[3] != NULL && result[4] == NULL
				&& result[5] == NULL && result[6] == NULL && result[7] == NULL
				&& hitListCmpNodes(result[0], node1) == 0
				&& hitListCmpNodes(result[0]->next, node2) == 0
				&& hitListCmpNodes(result[1], node3) == 0
				&& hitListCmpNodes(result[1]->next, node4) == 0
				&& hitListCmpNodes(result[2], node5) == 0
				&& hitListCmpNodes(result[3], node6) == 0)
		{}
		else
			fail("Incorrect behavior with a normal pair of queries "
					"(case 1).\n");
		qryDelete(qry1);
		qryDelete(qry2);
		hitListDelete(result[0]);
		hitListDelete(result[1]);
		hitListDelete(result[2]);
		hitListDelete(result[3]);
		hitListDelete(result[4]);
		hitListDelete(result[5]);
		hitListDelete(result[6]);
		hitListDelete(result[7]);
		free(result);
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		hitListDelete(node4);
		hitListDelete(node5);
		hitListDelete(node6);
		if (hashTable != NULL)
		{
			int i;
			int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE, seedLen);
			for (i = 0; i < numRefTuples; ++i)
				refPosListDelete(hashTable[i]);
			free(hashTable);
		}
		remove(refFile);
	}
}
END_TEST


/**
 * Tests refPreprocess function.
 */
START_TEST(preprocess)
{
	/* Test behavior when seed length is 0. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref1.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGATCTACTACTA", filePtr);
		fclose(filePtr);
		uint seedLength = 0;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);
		if (hashTable != NULL)
			fail("Incorrect behavior when seed length is 0.\n");
		remove(refFile);
	}

	/* Test behavior when reference size is less than seed length. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref3.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCA\n", filePtr);
		fclose(filePtr);
		uint seedLength = 10;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);
		int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE,
				(int) seedLength);
		int i;
		for (i = 0; i < numRefTuples; ++i)
		{
			if (hashTable[i] != NULL)
				fail("Incorrect behavior when reference size is less "
						"than seed length.\n");
		}
		remove(refFile);
		free(hashTable);
	}

	/* Test behavior when reference size is equal to seed length. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref4.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCA", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);
		RefPosList *node = refPosListCreate(0, 0);
		int nodeHash = getHash("ACGTTGCA", 8);
		if (hashTable == NULL
				|| refPosListIsEqual(hashTable[nodeHash], node) != 0)
			fail("Incorrect behavior when reference size is equal to "
					"seed length.\n");
		remove(refFile);
		if (hashTable != NULL)
		{
			int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE,
					(int) seedLength);
			int i;
			for (i = 0; i < numRefTuples; ++i)
				refPosListDelete(hashTable[i]);
			free(hashTable);
		}
		refPosListDelete(node);
	}

	/* Test behavior when reference size is greater than seed length. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref5.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCATTACGACGACTG\n", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);
		RefPosList *node1 = refPosListCreate(0, 0);
		int node1Hash = getHash("ACGTTGCA", 8);
		RefPosList *node2 = refPosListCreate(0, 8);
		int node2Hash = getHash("TTACGACG", 8);
		if (hashTable == NULL
				|| refPosListIsEqual(hashTable[node1Hash], node1) != 0
				|| refPosListIsEqual(hashTable[node2Hash], node2) != 0)
			fail("Incorrect behavior when reference size is greater than "
					"seed length.\n");
		remove(refFile);
		if (hashTable != NULL)
		{
			int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE,
					(int) seedLength);
			int i;
			for (i = 0; i < numRefTuples; ++i)
				refPosListDelete(hashTable[i]);
			free(hashTable);
		}
		refPosListDelete(node1);
		refPosListDelete(node2);
	}

	/* Test behavior when there are multiple references and reference
	 * sequences are all in 1 line. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref6.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCATTACGACGACTG\n"
				">ref2\nTTACTACCGGAGGTACGTAC", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);
		RefPosList *node1 = refPosListCreate(0, 0);
		int node1Hash = getHash("ACGTTGCA", 8);
		RefPosList *node2 = refPosListCreate(0, 8);
		int node2Hash = getHash("TTACGACG", 8);
		RefPosList *node3 = refPosListCreate(1, 0);
		int node3Hash = getHash("TTACTACC", 8);
		RefPosList *node4 = refPosListCreate(1, 8);
		int node4Hash = getHash("GGAGGTAC", 8);
		if (hashTable == NULL
				|| refPosListIsEqual(hashTable[node1Hash], node1) != 0
				|| refPosListIsEqual(hashTable[node2Hash], node2) != 0
				|| refPosListIsEqual(hashTable[node3Hash], node3) != 0
				|| refPosListIsEqual(hashTable[node4Hash], node4) != 0)
			fail("Incorrect behavior when there are multiple references "
					"and reference sequences are all in 1 line.\n");
		remove(refFile);
		if (hashTable != NULL)
		{
			int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE,
					(int) seedLength);
			int i;
			for (i = 0; i < numRefTuples; ++i)
				refPosListDelete(hashTable[i]);
			free(hashTable);
		}
		refPosListDelete(node1);
		refPosListDelete(node2);
		refPosListDelete(node3);
		refPosListDelete(node4);
	}

	/* Test behavior when there are multiple references and reference
	 * sequences span multiple lines. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref7.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCATT\nACGACGACTG\n"
				">ref2\nTTACTACCGG\nAGGTACGTAC", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);
		RefPosList *node1 = refPosListCreate(0, 0);
		int node1Hash = getHash("ACGTTGCA", 8);
		RefPosList *node2 = refPosListCreate(0, 8);
		int node2Hash = getHash("TTACGACG", 8);
		RefPosList *node3 = refPosListCreate(1, 0);
		int node3Hash = getHash("TTACTACC", 8);
		RefPosList *node4 = refPosListCreate(1, 8);
		int node4Hash = getHash("GGAGGTAC", 8);
		if (hashTable == NULL
				|| refPosListIsEqual(hashTable[node1Hash], node1) != 0
				|| refPosListIsEqual(hashTable[node2Hash], node2) != 0
				|| refPosListIsEqual(hashTable[node3Hash], node3) != 0
				|| refPosListIsEqual(hashTable[node4Hash], node4) != 0)
			fail("Incorrect behavior when there are multiple references "
					"and reference sequences span multiple lines.\n");
		remove(refFile);
		if (hashTable != NULL)
		{
			int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE,
					(int) seedLength);
			int i;
			for (i = 0; i < numRefTuples; ++i)
				refPosListDelete(hashTable[i]);
			free(hashTable);
		}
		refPosListDelete(node1);
		refPosListDelete(node2);
		refPosListDelete(node3);
		refPosListDelete(node4);
	}

	/* Test behavior when there is an empty line between subsequent
	 * references. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref8.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCATT\nACGACGACTG\n\n"
				">ref2\nTTACTACCGG\nAGGTACGTAC", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);
		RefPosList *node1 = refPosListCreate(0, 0);
		int node1Hash = getHash("ACGTTGCA", 8);
		RefPosList *node2 = refPosListCreate(0, 8);
		int node2Hash = getHash("TTACGACG", 8);
		RefPosList *node3 = refPosListCreate(1, 0);
		int node3Hash = getHash("TTACTACC", 8);
		RefPosList *node4 = refPosListCreate(1, 8);
		int node4Hash = getHash("GGAGGTAC", 8);
		if (hashTable == NULL
				|| refPosListIsEqual(hashTable[node1Hash], node1) != 0
				|| refPosListIsEqual(hashTable[node2Hash], node2) != 0
				|| refPosListIsEqual(hashTable[node3Hash], node3) != 0
				|| refPosListIsEqual(hashTable[node4Hash], node4) != 0)
			fail("Incorrect behavior when there is an empty line "
					"between subsequent references.\n");
		remove(refFile);
		if (hashTable != NULL)
		{
			int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE,
					(int) seedLength);
			int i;
			for (i = 0; i < numRefTuples; ++i)
				refPosListDelete(hashTable[i]);
			free(hashTable);
		}
		refPosListDelete(node1);
		refPosListDelete(node2);
		refPosListDelete(node3);
		refPosListDelete(node4);
	}
}
END_TEST


/**
 * Tests @a serializeHashTable2 function.
 */
START_TEST(serializeHashTable2_wrap)
{
	/* Number of reference tuples is 1. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCATT", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength) + 1;
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < (numKeys - 1); ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}

		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		uint tupleIgnoreThres = 10;
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		int hash = getHash("ACGTTGCA", seedLength);

		if (numKeys != 65537 || numVals != 1)
		{
			fail("Incorrect behavior when number of reference tuples is 1 "
					"(case 0).\n");
		}

		for (i = 0; i < numKeys; ++i)
		{
			if (i == hash)
			{
				if (keys[i] != 0 || vals[0] != 0)
					fail("Incorrect behavior  when number of reference tuples "
							"is 1 (case 1).\n");
			}
			else if (i == (numKeys - 1))
			{
				if (keys[i] == UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 1 (case 2).\n");
			}
			else
			{
				if (keys[i] != UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 1 (case 3).\n");
			}
		}

		free(keys);
		free(vals);
		refDeleteHashTable(hashTable, (numKeys - 1));
	}

	/* Number of reference tuples is 1 and the tuple is 'AAAAA'. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nAAAAAATT", filePtr);
		fclose(filePtr);
		uint seedLength = 5;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength) + 1;
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < (numKeys - 1); ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}

		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		uint tupleIgnoreThres = 10;
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		int hash = getHash("AAAAA", seedLength);

		if (numKeys != 1025 || numVals != 1)
		{
			fail("Incorrect behavior when number of reference tuples is 1 "
					"and the tuple is 'AAAAA' (case 0).\n");
		}

		for (i = 0; i < numKeys; ++i)
		{
			if (i == hash)
			{
				if (keys[i] != 0 || vals[0] != 0)
					fail("Incorrect behavior  when number of reference tuples "
							"is 1 and the tuple is 'AAAAA' (case 1).\n");
			}
			else if (i == (numKeys - 1))
			{
				if (keys[i] == UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 1 and the tuple is 'AAAAA' (case 2).\n");
			}
			else
			{
				if (keys[i] != UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 1 and the tuple is 'AAAAA' (case 3).\n");
			}
		}

		free(keys);
		free(vals);
		refDeleteHashTable(hashTable, (numKeys - 1));
	}

	/* Number of reference tuples is 2. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATTCCCA", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength) + 1;
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < (numKeys - 1); ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}

		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		uint tupleIgnoreThres = 10;
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		int hash1 = getHash("ACGTACGT", seedLength);
		int hash2 = getHash("AATTCCCA", seedLength);

		if (numKeys != 65537 || numVals != 2)
		{
			fail("Incorrect behavior when number of reference tuples is 2 "
					"(case 0).\n");
		}

		for (i = 0; i < numKeys; ++i)
		{
			if (i == hash1)
			{
				if (keys[i] != 1 || vals[1] != 0)
					fail("Incorrect behavior  when number of reference tuples "
							"is 2 (case 1).\n");
			}
			else if (i == hash2)
			{
				if (keys[i] != 0 || vals[0] != 1)
					fail("Incorrect behavior  when number of reference tuples "
							"is 2 (case 2).\n");
			}
			else if (i == (numKeys - 1))
			{
				if (keys[i] == UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 2 (case 3).\n");
			}
			else
			{
				if (keys[i] != UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 2 (case 4).\n");
			}
		}

		free(keys);
		free(vals);
		refDeleteHashTable(hashTable, (numKeys - 1));
	}

	/* Number of reference tuples is 4. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCATTACGACGACTG\n"
				">ref2\nTTACTACCGG\nAGGTACGTAC", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength) + 1;
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < (numKeys - 1); ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}

		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		uint tupleIgnoreThres = 10;
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		int hash1 = getHash("ACGTTGCA", seedLength);
		int hash2 = getHash("TTACGACG", seedLength);
		int hash3 = getHash("TTACTACC", seedLength);
		int hash4 = getHash("GGAGGTAC", seedLength);

		if (numKeys != 65537 || numVals != 4)
		{
			fail("Incorrect behavior when number of reference tuples is 4 "
					"(case 0).\n");
		}

		for (i = 0; i < numKeys; ++i)
		{
			if (i == hash1)
			{
				if (keys[i] != 0 || vals[0] != 0)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 (case 1).\n");
			}
			else if (i == hash2)
			{
				if (keys[i] != 3 || vals[3] != 1)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 (case 2).\n");
			}
			else if (i == hash3)
			{
				if (keys[i] != 2 || vals[2] != 134217728)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 (case 3).\n");
			}
			else if (i == hash4)
			{
				if (keys[i] != 1 || vals[1] != 134217729)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 (case 4).\n");
			}
			else if (i == (numKeys - 1))
			{
				if (keys[i] == UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 2 (case 3).\n");
			}
			else
			{
				if (keys[i] != UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 (case 5).\n");
			}
		}

		free(keys);
		free(vals);
		refDeleteHashTable(hashTable, (numKeys - 1));
	}

	/* Number of reference tuples is 4 and seed length is 5. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nAAAAACCCCC\n"
				">ref2\nGGGGGAAAAA\n", filePtr);
		fclose(filePtr);
		uint seedLength = 5;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength) + 1;
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < (numKeys - 1); ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}

		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		uint tupleIgnoreThres = 10;
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		int hash1 = getHash("AAAAA", seedLength);
		int hash2 = getHash("CCCCC", seedLength);
		int hash3 = getHash("GGGGG", seedLength);

		if (numKeys != 1025 || numVals != 4)
		{
			fail("Incorrect behavior when number of reference tuples is 4 "
					"and seed length is 5 (case 0).\n");
		}

		for (i = 0; i < numKeys; ++i)
		{
			if (i == hash1)
			{
				if (keys[i] != 0 || vals[0] != 0)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 and seed length is 5 (case 1).\n");
			}
			else if (i == hash2)
			{
				if (keys[i] != 2 || vals[2] != 1)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 and seed length is 5 (case 2).\n");
			}
			else if (i == hash3)
			{
				if (keys[i] != 3 || vals[3] != 134217728)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 and seed length is 5 (case 3).\n");
			}
			else if (i == (numKeys - 1))
			{
				if (keys[i] == UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 and seed length is 5 (case 4).\n");
			}
			else
			{
				if (keys[i] != UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 and seed length is 5 (case 5).\n");
			}
		}

		free(keys);
		free(vals);
		refDeleteHashTable(hashTable, (numKeys - 1));
	}

	/* Number of reference tuples is 9. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCATTACGACGACTG\n"
				">ref2\nTTACTACCGG\nAGGTACTTACGACGGGATCTAA\n"
				"TTACGACGTTACGACGTTACGACG", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength) + 1;
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < (numKeys - 1); ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}

		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		uint tupleIgnoreThres = 10;
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		int hash1 = getHash("ACGTTGCA", seedLength);
		int hash2 = getHash("TTACGACG", seedLength);
		int hash3 = getHash("TTACTACC", seedLength);
		int hash4 = getHash("GGAGGTAC", seedLength);
		int hash5 = getHash("GGATCTAA", seedLength);

		if (numKeys != 65537 || numVals != 9)
		{
			fail("Incorrect behavior when number of reference tuples is 9 "
					"(case 0).\n");
		}

		for (i = 0; i < numKeys; ++i)
		{
			if (i == hash1)
			{
				if (keys[i] != 1 || vals[1] != 0)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 (case 1).\n");
			}
			else if (i == hash2)
			{
				if (keys[i] != 4 || vals[4] != 1 || vals[5] != 134217730
						|| vals[6] != 134217732 || vals[7] != 134217733
						|| vals[8] != 134217734)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 (case 2).\n");
			}
			else if (i == hash3)
			{
				if (keys[i] != 3 || vals[3] != 134217728)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 (case 3).\n");
			}
			else if (i == hash4)
			{
				if (keys[i] != 2 || vals[2] != 134217729)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 (case 4).\n");
			}
			else if (i == hash5)
			{
				if (keys[i] != 0 || vals[0] != 134217731)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 (case 5).\n");
			}
			else if (i == (numKeys - 1))
			{
				if (keys[i] == UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 (case 6).\n");
			}
			else
			{
				if (keys[i] != UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 (case 7).\n");
			}
		}

		free(keys);
		free(vals);
		refDeleteHashTable(hashTable, (numKeys - 1));
	}

	/* Number of reference tuples is 9; one of the tuples has more number
	 * of repeats than allowed by the threshold. */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCATTACGACGACGTTGCAATT\n"
				">ref2\nTTACTACCGG\nAGGTACTTACGACGGGATCTAA\n"
				"TTACGACGTTACGACGTTACGACG", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength) + 1;
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < (numKeys - 1); ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}

		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		uint tupleIgnoreThres = 4;
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		int hash1 = getHash("ACGTTGCA", seedLength);
		int hash2 = getHash("TTACGACG", seedLength);
		int hash3 = getHash("TTACTACC", seedLength);
		int hash4 = getHash("GGAGGTAC", seedLength);
		int hash5 = getHash("GGATCTAA", seedLength);

		if (numKeys != 65537 || numVals != 10)
		{
			fail("Incorrect behavior when number of reference tuples is 9 "
					"and one of the tuples has more number of repeats "
					"than allowed by the threshold (case 0).\n");
		}

		for (i = 0; i < numKeys; ++i)
		{
			if (i == hash1)
			{
				if (keys[i] != 1 || vals[1] != 0 || vals[2] != 2)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 and one of the tuples has more number of "
							"repeats than allowed by the threshold (case 1)."
							"\n");
			}
			else if (i == hash3)
			{
				if (keys[i] != 4 || vals[4] != 134217728)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 and one of the tuples has more number of "
							"repeats than allowed by the threshold (case 2). "
							"\n");
			}
			else if (i == hash4)
			{
				if (keys[i] != 3 || vals[3] != 134217729)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 and one of the tuples has more number of "
							"repeats  than allowed by the threshold (case 3). "
							"\n");
			}
			else if (i == hash5)
			{
				if (keys[i] != 0 || vals[0] != 134217731)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 and one of the tuples has more number of "
							"repeats than allowed by the threshold (case 4)."
							"\n");
			}
			else if (i == (numKeys - 1))
			{
				if (keys[i] == UINT_MAX)
					fail("Incorrect behavior when number of reference tuples "
							"is 9 and one of the tuples has more number of "
							"repeats than allowed by the threshold (case 5)."
							"\n");
			}
			else
			{
				if (keys[i] != UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 9 and one of the tuples has more number of "
							"repeats than allowed by the threshold (case 6). "
							"\n");
			}
		}

		free(keys);
		free(vals);
		refDeleteHashTable(hashTable, (numKeys - 1));
	}

	/* Number of reference tuples is 4 (2 unique tuples and 2 duplicates). */
	{
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">r1\nACGTACGTAATTCCCA\n>r2\nAATTCCCAACGTACGT", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength) + 1;
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < (numKeys - 1); ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}

		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		uint tupleIgnoreThres = 4;
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		int hash1 = getHash("AATTCCCA", seedLength);
		int hash2 = getHash("ACGTACGT", seedLength);

		if (numKeys != 65537 && numVals != 4)
		{
			fail("Incorrect behavior when number of reference tuples is 4 "
					"(2 unique tuples and 2 duplicates) (case 0).\n");
		}

		for (i = 0; i < numKeys; ++i)
		{
			if (i == hash1)
			{
				if (keys[i] != 0 || vals[0] != 1 || vals[1] != 134217728)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 (2 unique tuples and 2 duplicates) (case 1)."
							"\n");
			}
			else if (i == hash2)
			{
				if (keys[i] != 2 || vals[2] != 0 || vals[3] != 134217729)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 (2 unique tuples and 2 duplicates) (case 2)."
							"\n");
			}
			else if (i == (numKeys - 1))
			{
				if (keys[i] == UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 (2 unique tuples and 2 duplicates) (case 3)."
							"\n");
			}
			else
			{
				if (keys[i] != UINT_MAX)
					fail("Incorrect behavior  when number of reference tuples "
							"is 4 (2 unique tuples and 2 duplicates) (case 4)."
							"\n");
			}
		}

		free(keys);
		free(vals);
		refDeleteHashTable(hashTable, (numKeys - 1));
	}
}
END_TEST


/**
 * Tests @a pow_gpu function.
 */
START_TEST(pow_gpu)
{
	/* Base is 2 and exponent is 5. */
	{
		int base = 2;
		int exponent = 5;
		int *p_d;
		cudaMalloc(&p_d, sizeof(int));
		PRINT_CUDA_ERROR()

		pow_gpu_wrap<<<1, 1>>>(base, exponent, p_d);
		PRINT_CUDA_ERROR()

		int p;
		cudaMemcpy(&p, p_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (p != 32)
			fail("Incorrect behavior when base is 2 and exponent is 5.\n");

		cudaFree(p_d);
	}

	/* Base is 3 and exponent is 3. */
	{
		int base = 3;
		int exponent = 3;
		int *p_d;
		cudaMalloc(&p_d, sizeof(int));
		PRINT_CUDA_ERROR()

		pow_gpu_wrap<<<1, 1>>>(base, exponent, p_d);
		PRINT_CUDA_ERROR()

		int p;
		cudaMemcpy(&p, p_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (p != 27)
			fail("Incorrect behavior when base is 3 and exponent is 3.\n");

		cudaFree(p_d);
	}
}
END_TEST


/**
 * Tests @a refSearchQuery2_gpu function.
 */
START_TEST(refSearchQuery2_gpu)
{
	/* 2 references and 4 queries. */
	{
		/* Preprocess reference file and create a hash table. */
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCAATTACGTAGCTTACGACGACTGTC\nACCATGCACTATGC\n"
				">ref2\nTTACTACCTTACCAGAGG\nCAATTTCAAC", filePtr);
		fclose(filePtr);
		uint seedLength = 8;
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		/* Serialize the hash table. */
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength);
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < numKeys; ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}
		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		uint tupleIgnoreThres = 2;
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		/* Copy serialized hash table to GPU. */
		uint *keys_d;
		cudaMalloc(&keys_d, numKeys * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(keys_d, keys, numKeys * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		uint *vals_d;
		cudaMalloc(&vals_d, numVals * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(vals_d, vals, numVals * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Copy queries to GPU. */
		int numQrs = 4;
		int qryLen = 10;
		int totalQrySeqLen = (numQrs * qryLen * sizeof(char)) + sizeof(char);
		char *qrySeq = (char *) malloc(totalQrySeqLen);
		strcpy(qrySeq, "ATTACGTAGC");
		strcat(qrySeq, "ACCATGCACT");
		strcat(qrySeq, "GCAATTTCAG");
		strcat(qrySeq, "TTTTACCAGA");
		char *qrySeq_d;
		cudaMalloc(&qrySeq_d, totalQrySeqLen);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qrySeq_d, qrySeq, totalQrySeqLen, cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Allocate GPU memory for results. */
		uint maxNumHits = 1;
		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, numQrs * maxNumHits * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		/* Search queries. */
		dim3 dimGrid, dimBlock;
		dimBlock.x = 1;
		dimBlock.y = 1;
		dimBlock.z = 1;
		dimGrid.x = numQrs;
		dimGrid.y = 1;
		int randNum = rand();
		refSearchQuery2_kernel<<<dimGrid, dimBlock>>>(keys_d, vals_d, numKeys,
				numVals, qrySeq_d, refIdxArr_d, shiftArr_d, posArr_d, qryLen,
				seedLength, maxNumHits, randNum);
		PRINT_CUDA_ERROR()

		/* Copy results from GPU to CPU. */
		char *refIdxArr = (char *) malloc(numQrs * maxNumHits * sizeof(char));
		cudaMemcpy(refIdxArr, refIdxArr_d, (numQrs * maxNumHits * sizeof(char)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *shiftArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(shiftArr, shiftArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *posArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(posArr, posArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		/* Check results. */
		if (refIdxArr[0] != 0 || shiftArr[0] != 8 || posArr[0] != 8
				|| refIdxArr[1] != 0 || shiftArr[1] != 32 || posArr[1] != 32
				|| refIdxArr[2] != -1 || shiftArr[2] != -1 || posArr[2] != -1
				|| refIdxArr[3] != 1 || shiftArr[3] != 6 || posArr[3] != 8)
			fail("Incorrect behavior when there are 2 references and "
					"4 queries.\n");

		remove(refFile);
		refDeleteHashTable(hashTable, numKeys);
		free(keys);
		free(vals);
		free(qrySeq);
		free(refIdxArr);
		free(shiftArr);
		free(posArr);
		cudaFree(keys_d);
		cudaFree(vals_d);
		cudaFree(qrySeq_d);
		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* 2 references and 3 queries. */
	{
		uint seedLength = 5;
		int numQrs = 3;
		int qryLen = 15;
		uint maxNumHits = 1;
		uint tupleIgnoreThres = 2;

		/* Preprocess reference file and create a hash table. */
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCAATTACGTAGCTTACGAC\n"
				">ref2\nACGTTGCAATTACGTAGCTTACGAC",
				filePtr);
		fclose(filePtr);
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		/* Serialize the hash table. */
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength);
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < numKeys; ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}
		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		/* Copy serialized hash table to GPU. */
		uint *keys_d;
		cudaMalloc(&keys_d, numKeys * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(keys_d, keys, numKeys * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		uint *vals_d;
		cudaMalloc(&vals_d, numVals * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(vals_d, vals, numVals * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Copy queries to GPU. */
		int totalQrySeqLen = (numQrs * qryLen * sizeof(char)) + sizeof(char);
		char *qrySeq = (char *) malloc(totalQrySeqLen);
		strcpy(qrySeq, "GCAATCCACCAGCTT");
		strcat(qrySeq, "GCAATTACGTGGGAT");
		strcat(qrySeq, "TACGTAGCTTACGAC");
		char *qrySeq_d;
		cudaMalloc(&qrySeq_d, totalQrySeqLen);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qrySeq_d, qrySeq, totalQrySeqLen, cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Allocate GPU memory for results. */
		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, numQrs * maxNumHits * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		/* Search queries. */
		dim3 dimGrid, dimBlock;
		dimBlock.x = 1;
		dimBlock.y = 1;
		dimBlock.z = 1;
		dimGrid.x = numQrs;
		dimGrid.y = 1;
		int randNum = rand();
		refSearchQuery2_kernel<<<dimGrid, dimBlock>>>(keys_d, vals_d, numKeys,
				numVals, qrySeq_d, refIdxArr_d, shiftArr_d, posArr_d, qryLen,
				seedLength, maxNumHits, randNum);
		PRINT_CUDA_ERROR()

		/* Copy results from GPU to CPU. */
		char *refIdxArr = (char *) malloc(numQrs * maxNumHits * sizeof(char));
		cudaMemcpy(refIdxArr, refIdxArr_d, (numQrs * maxNumHits * sizeof(char)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *shiftArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(shiftArr, shiftArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *posArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(posArr, posArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		/* Check results. */
		if (((refIdxArr[0] != 0 || shiftArr[0] != 5 || posArr[0] != 5)
				&& (refIdxArr[0] != 1 || shiftArr[0] != 5 || posArr[0] != 5))
				|| ((refIdxArr[1] != 0 || shiftArr[1] != 5 || posArr[1] != 5)
						&& (refIdxArr[1] != 1 || shiftArr[1] != 5
								|| posArr[1] != 5))
				|| ((refIdxArr[2] != 0 || shiftArr[2] != 10 || posArr[2] != 10)
						&& (refIdxArr[2] != 1 || shiftArr[2] != 10
								|| posArr[2] != 10)))
			fail("Incorrect behavior when there are 2 references and 3 "
					"queries.\n");

		remove(refFile);
		refDeleteHashTable(hashTable, numKeys);
		free(keys);
		free(vals);
		free(qrySeq);
		free(refIdxArr);
		free(shiftArr);
		free(posArr);
		cudaFree(keys_d);
		cudaFree(vals_d);
		cudaFree(qrySeq_d);
		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* 4 references and 3 queries. */
	{
		uint seedLength = 5;
		int numQrs = 3;
		int qryLen = 15;
		uint maxNumHits = 2;
		uint tupleIgnoreThres = 2;

		/* Preprocess reference file and create a hash table. */
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTTGCAATTACGTAGCTTACGAC\n"
				">ref2\nACGTTGCAATTACGTAGCTTACGAC\n"
				">ref3\nACGTTGCAATTAGGTCGCTTACGAC\n"
				">ref4\nACGTTGCAATTAAGTTGCTTACGAC\n",
				filePtr);
		fclose(filePtr);
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		/* Serialize the hash table. */
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength);
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < numKeys; ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}
		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		/* Copy serialized hash table to GPU. */
		uint *keys_d;
		cudaMalloc(&keys_d, numKeys * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(keys_d, keys, numKeys * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		uint *vals_d;
		cudaMalloc(&vals_d, numVals * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(vals_d, vals, numVals * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Copy queries to GPU. */
		int totalQrySeqLen = (numQrs * qryLen * sizeof(char)) + sizeof(char);
		char *qrySeq = (char *) malloc(totalQrySeqLen);
		strcpy(qrySeq, "GCAATCCACCAGCTT");
		strcat(qrySeq, "GCAATTACGTGGGAT");
		strcat(qrySeq, "TACGTACGTTACGAC");
		char *qrySeq_d;
		cudaMalloc(&qrySeq_d, totalQrySeqLen);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qrySeq_d, qrySeq, totalQrySeqLen, cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Allocate GPU memory for results. */
		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, numQrs * maxNumHits * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		/* Search queries. */
		dim3 dimGrid, dimBlock;
		dimBlock.x = 1;
		dimBlock.y = 1;
		dimBlock.z = 1;
		dimGrid.x = numQrs;
		dimGrid.y = 1;
		int randNum = rand();
		refSearchQuery2_kernel<<<dimGrid, dimBlock>>>(keys_d, vals_d, numKeys,
				numVals, qrySeq_d, refIdxArr_d, shiftArr_d, posArr_d, qryLen,
				seedLength, maxNumHits, randNum);
		PRINT_CUDA_ERROR()

		/* Copy results from GPU to CPU. */
		char *refIdxArr = (char *) malloc(numQrs * maxNumHits * sizeof(char));
		cudaMemcpy(refIdxArr, refIdxArr_d, (numQrs * maxNumHits * sizeof(char)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *shiftArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(shiftArr, shiftArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *posArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(posArr, posArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		/* Check results. */
		if (((refIdxArr[0] != 0 || shiftArr[0] != 5 || posArr[0] != 15)
				&& (refIdxArr[0] != 1 || shiftArr[0] != 5 || posArr[0] != 15))
				|| ((refIdxArr[1] != 0 || shiftArr[1] != 5 || posArr[1] != 15)
						&& (refIdxArr[1] != 1 || shiftArr[1] != 5
								|| posArr[1] != 15))
				|| ((refIdxArr[2] != 0 || shiftArr[2] != 5 || posArr[2] != 10)
						&& (refIdxArr[2] != 1 || shiftArr[2] != 5
								|| posArr[2] != 10))
				|| ((refIdxArr[3] != 0 || shiftArr[3] != 5 || posArr[3] != 10)
						&& (refIdxArr[3] != 1 || shiftArr[3] != 5
								|| posArr[3] != 10))
				|| ((refIdxArr[4] != 0 || shiftArr[4] != 10 || posArr[4] != 10)
						&& (refIdxArr[4] != 1 || shiftArr[4] != 10
								|| posArr[4] != 10)
						&& (refIdxArr[4] != 0 || shiftArr[4] != 6
								|| posArr[4] != 10)
						&& (refIdxArr[4] != 1 || shiftArr[4] != 6
								|| posArr[4] != 10))
				|| ((refIdxArr[5] != 0 || shiftArr[5] != 10 || posArr[5] != 10)
						&& (refIdxArr[5] != 1 || shiftArr[5] != 10
								|| posArr[3] != 10)
						&& (refIdxArr[5] != 0 || shiftArr[5] != 6
								|| posArr[3] != 10)
						&& (refIdxArr[5] != 1 || shiftArr[5] != 6
								|| posArr[3] != 10)))
			fail("Incorrect behavior when there are 4 references and 3 "
					"queries.\n");

		remove(refFile);
		refDeleteHashTable(hashTable, numKeys);
		free(keys);
		free(vals);
		free(qrySeq);
		free(refIdxArr);
		free(shiftArr);
		free(posArr);
		cudaFree(keys_d);
		cudaFree(vals_d);
		cudaFree(qrySeq_d);
		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* 1 reference and 1 query. */
	{
		uint seedLength = 8;
		int numQrs = 1;
		int qryLen = 10;
		uint maxNumHits = 1;
		uint tupleIgnoreThres = 2;

		/* Preprocess reference file and create a hash table. */
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATTCCCA", filePtr);
		fclose(filePtr);
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		/* Serialize the hash table. */
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength) + 1;
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < numKeys; ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}
		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		/* Copy serialized hash table to GPU. */
		uint *keys_d;
		cudaMalloc(&keys_d, numKeys * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(keys_d, keys, numKeys * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		uint *vals_d;
		cudaMalloc(&vals_d, numVals * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(vals_d, vals, numVals * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Copy queries to GPU. */
		int totalQrySeqLen = (numQrs * qryLen * sizeof(char)) + sizeof(char);
		char *qrySeq = (char *) malloc(totalQrySeqLen);
		strcpy(qrySeq, "ACGTACGTGG");
		char *qrySeq_d;
		cudaMalloc(&qrySeq_d, totalQrySeqLen);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qrySeq_d, qrySeq, totalQrySeqLen, cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Allocate GPU memory for results. */
		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, numQrs * maxNumHits * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		/* Search queries. */
		dim3 dimGrid, dimBlock;
		dimBlock.x = 1;
		dimBlock.y = 1;
		dimBlock.z = 1;
		dimGrid.x = numQrs;
		dimGrid.y = 1;
		int randNum = rand();
		refSearchQuery2_kernel<<<dimGrid, dimBlock>>>(keys_d, vals_d,
				numKeys, numVals, qrySeq_d, refIdxArr_d, shiftArr_d, posArr_d,
				qryLen, seedLength, maxNumHits, randNum);
		PRINT_CUDA_ERROR()

		/* Copy results from GPU to CPU. */
		char *refIdxArr = (char *) malloc(numQrs * maxNumHits * sizeof(char));
		cudaMemcpy(refIdxArr, refIdxArr_d, (numQrs * maxNumHits * sizeof(char)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *shiftArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(shiftArr, shiftArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *posArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(posArr, posArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		/* Check results. */
		if (refIdxArr[0] != 0 || shiftArr[0] != 0 || posArr[0] != 0)
			fail("Incorrect behavior when there is 1 reference and 1 query.\n");

		remove(refFile);
		refDeleteHashTable(hashTable, numKeys);
		free(keys);
		free(vals);
		free(qrySeq);
		free(refIdxArr);
		free(shiftArr);
		free(posArr);
		cudaFree(keys_d);
		cudaFree(vals_d);
		cudaFree(qrySeq_d);
		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* 1 reference and 2 queries. */
	{
		uint seedLength = 8;
		int numQrs = 2;
		int qryLen = 10;
		uint maxNumHits = 1;
		uint tupleIgnoreThres = 2;

		/* Preprocess reference file and create a hash table. */
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATTCCCA", filePtr);
		fclose(filePtr);
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		/* Serialize the hash table. */
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength);
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < numKeys; ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}
		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		/* Copy serialized hash table to GPU. */
		uint *keys_d;
		cudaMalloc(&keys_d, numKeys * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(keys_d, keys, numKeys * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		uint *vals_d;
		cudaMalloc(&vals_d, numVals * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(vals_d, vals, numVals * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Copy queries to GPU. */
		int totalQrySeqLen = (numQrs * qryLen * sizeof(char)) + sizeof(char);
		char *qrySeq = (char *) malloc(totalQrySeqLen);
		strcpy(qrySeq, "ACGTACGTGG");
		strcat(qrySeq, "ACAATTCCCA");
		char *qrySeq_d;
		cudaMalloc(&qrySeq_d, totalQrySeqLen);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qrySeq_d, qrySeq, totalQrySeqLen, cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Allocate GPU memory for results. */
		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, numQrs * maxNumHits * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		/* Search queries. */
		dim3 dimGrid, dimBlock;
		dimBlock.x = 1;
		dimBlock.y = 1;
		dimBlock.z = 1;
		dimGrid.x = numQrs;
		dimGrid.y = 1;
		int randNum = rand();
		refSearchQuery2_kernel<<<dimGrid, dimBlock>>>(keys_d, vals_d,
				numKeys, numVals, qrySeq_d, refIdxArr_d, shiftArr_d, posArr_d,
				qryLen, seedLength, maxNumHits, randNum);
		PRINT_CUDA_ERROR()

		/* Copy results from GPU to CPU. */
		char *refIdxArr = (char *) malloc(numQrs * maxNumHits * sizeof(char));
		cudaMemcpy(refIdxArr, refIdxArr_d, (numQrs * maxNumHits * sizeof(char)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *shiftArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(shiftArr, shiftArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *posArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(posArr, posArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		/* Check results. */
		if (refIdxArr[0] != 0 || shiftArr[0] != 0 || posArr[0] != 0
				|| refIdxArr[1] != 0 || shiftArr[1] != 6 || posArr[1] != 8)
			fail("Incorrect behavior when there is 1 reference and 2 queries."
					"\n");

		remove(refFile);
		refDeleteHashTable(hashTable, numKeys);
		free(keys);
		free(vals);
		free(qrySeq);
		free(refIdxArr);
		free(shiftArr);
		free(posArr);
		cudaFree(keys_d);
		cudaFree(vals_d);
		cudaFree(qrySeq_d);
		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* 1 reference and 2 queries and max number of hits is 3. */
	{
		uint seedLength = 8;
		int numQrs = 2;
		int qryLen = 10;
		uint maxNumHits = 3;
		uint tupleIgnoreThres = 2;

		/* Preprocess reference file and create a hash table. */
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATTCCCA", filePtr);
		fclose(filePtr);
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		/* Serialize the hash table. */
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength);
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < numKeys; ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}
		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		/* Copy serialized hash table to GPU. */
		uint *keys_d;
		cudaMalloc(&keys_d, numKeys * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(keys_d, keys, numKeys * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		uint *vals_d;
		cudaMalloc(&vals_d, numVals * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(vals_d, vals, numVals * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Copy queries to GPU. */
		int totalQrySeqLen = (numQrs * qryLen * sizeof(char)) + sizeof(char);
		char *qrySeq = (char *) malloc(totalQrySeqLen);
		strcpy(qrySeq, "ACGTACGTGG");
		strcat(qrySeq, "ACAATTCCCA");
		char *qrySeq_d;
		cudaMalloc(&qrySeq_d, totalQrySeqLen);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qrySeq_d, qrySeq, totalQrySeqLen, cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Allocate GPU memory for results. */
		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, numQrs * maxNumHits * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		/* Search queries. */
		dim3 dimGrid, dimBlock;
		dimBlock.x = 1;
		dimBlock.y = 1;
		dimBlock.z = 1;
		dimGrid.x = numQrs;
		dimGrid.y = 1;
		int randNum = rand();
		refSearchQuery2_kernel<<<dimGrid, dimBlock>>>(keys_d, vals_d,
				numKeys, numVals, qrySeq_d, refIdxArr_d, shiftArr_d, posArr_d,
				qryLen, seedLength, maxNumHits, randNum);
		PRINT_CUDA_ERROR()

		/* Copy results from GPU to CPU. */
		char *refIdxArr = (char *) malloc(numQrs * maxNumHits * sizeof(char));
		cudaMemcpy(refIdxArr, refIdxArr_d, (numQrs * maxNumHits * sizeof(char)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *shiftArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(shiftArr, shiftArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *posArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(posArr, posArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		/* Check results. */
		if (((refIdxArr[0] != 0 || shiftArr[0] != 0 || posArr[0] != 0)
				&& (refIdxArr[1] != -1 || shiftArr[1] != -1 || posArr[1] != -1)
				&& (refIdxArr[2] != -1 || shiftArr[2] != -1 || posArr[2] != -1))
				|| ((refIdxArr[3] != 0 || shiftArr[3] != 6 || posArr[3] != 8)
						&& (refIdxArr[4] != -1 || shiftArr[4] != -1
								|| posArr[4] != -1)
						&& (refIdxArr[5] != -1 || shiftArr[5] != -1
								|| posArr[5] != -1)))
			fail("Incorrect behavior when there is 1 reference and 2 queries "
					"and max number of hits is 3.\n");

		remove(refFile);
		refDeleteHashTable(hashTable, numKeys);
		free(keys);
		free(vals);
		free(qrySeq);
		free(refIdxArr);
		free(shiftArr);
		free(posArr);
		cudaFree(keys_d);
		cudaFree(vals_d);
		cudaFree(qrySeq_d);
		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* 1 reference and 4 queries. */
	{
		uint seedLength = 8;
		int numQrs = 4;
		int qryLen = 10;
		uint maxNumHits = 1;
		uint tupleIgnoreThres = 2;

		/* Preprocess reference file and create a hash table. */
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATTCCCA", filePtr);
		fclose(filePtr);
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		/* Serialize the hash table. */
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength);
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < numKeys; ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}
		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		/* Copy serialized hash table to GPU. */
		uint *keys_d;
		cudaMalloc(&keys_d, numKeys * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(keys_d, keys, numKeys * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		uint *vals_d;
		cudaMalloc(&vals_d, numVals * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(vals_d, vals, numVals * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Copy queries to GPU. */
		int totalQrySeqLen = (numQrs * qryLen * sizeof(char)) + sizeof(char);
		char *qrySeq = (char *) malloc(totalQrySeqLen);
		strcpy(qrySeq, "ACGTACGTGG");
		strcat(qrySeq, "CCACGTACGT");
		strcat(qrySeq, "ACAATTCCCA");
		strcat(qrySeq, "TGGGAATTGT");
		char *qrySeq_d;
		cudaMalloc(&qrySeq_d, totalQrySeqLen);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qrySeq_d, qrySeq, totalQrySeqLen, cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Allocate GPU memory for results. */
		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, numQrs * maxNumHits * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		/* Search queries. */
		dim3 dimGrid, dimBlock;
		dimBlock.x = 1;
		dimBlock.y = 1;
		dimBlock.z = 1;
		dimGrid.x = numQrs;
		dimGrid.y = 1;
		int randNum = rand();
		refSearchQuery2_kernel<<<dimGrid, dimBlock>>>(keys_d, vals_d, numKeys,
				numVals, qrySeq_d, refIdxArr_d, shiftArr_d, posArr_d, qryLen,
				seedLength, maxNumHits, randNum);
		PRINT_CUDA_ERROR()

		/* Copy results from GPU to CPU. */
		char *refIdxArr = (char *) malloc(numQrs * maxNumHits * sizeof(char));
		cudaMemcpy(refIdxArr, refIdxArr_d, (numQrs * maxNumHits * sizeof(char)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *shiftArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(shiftArr, shiftArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *posArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(posArr, posArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		/* Check results. */
		if (refIdxArr[0] != 0 || shiftArr[0] != 0 || posArr[0] != 0
				|| refIdxArr[1] != 0 || shiftArr[1] != -2 || posArr[1] != 0
				|| refIdxArr[2] != 0 || shiftArr[2] != 6 || posArr[2] != 8
				|| refIdxArr[3] != -1 || shiftArr[3] != -1 || posArr[3] != -1)
			fail("Incorrect behavior when there is 1 reference and 4 queries."
					"\n");

		remove(refFile);
		refDeleteHashTable(hashTable, numKeys);
		free(keys);
		free(vals);
		free(qrySeq);
		free(refIdxArr);
		free(shiftArr);
		free(posArr);
		cudaFree(keys_d);
		cudaFree(vals_d);
		cudaFree(qrySeq_d);
		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* 1 reference and 4 queries and max number of hits is 3. */
	{
		uint seedLength = 8;
		int numQrs = 4;
		int qryLen = 10;
		uint maxNumHits = 3;
		uint tupleIgnoreThres = 2;

		/* Preprocess reference file and create a hash table. */
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">ref1\nACGTACGTAATTCCCA", filePtr);
		fclose(filePtr);
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		/* Serialize the hash table. */
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength);
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < numKeys; ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}
		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		/* Copy serialized hash table to GPU. */
		uint *keys_d;
		cudaMalloc(&keys_d, numKeys * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(keys_d, keys, numKeys * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		uint *vals_d;
		cudaMalloc(&vals_d, numVals * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(vals_d, vals, numVals * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Copy queries to GPU. */
		int totalQrySeqLen = (numQrs * qryLen * sizeof(char)) + sizeof(char);
		char *qrySeq = (char *) malloc(totalQrySeqLen);
		strcpy(qrySeq, "ACGTACGTGG");
		strcat(qrySeq, "CCACGTACGT");
		strcat(qrySeq, "ACAATTCCCA");
		strcat(qrySeq, "TGGGAATTGT");
		char *qrySeq_d;
		cudaMalloc(&qrySeq_d, totalQrySeqLen);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qrySeq_d, qrySeq, totalQrySeqLen, cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Allocate GPU memory for results. */
		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, numQrs * maxNumHits * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		/* Search queries. */
		dim3 dimGrid, dimBlock;
		dimBlock.x = 1;
		dimBlock.y = 1;
		dimBlock.z = 1;
		dimGrid.x = numQrs;
		dimGrid.y = 1;
		int randNum = rand();
		refSearchQuery2_kernel<<<dimGrid, dimBlock>>>(keys_d, vals_d, numKeys,
				numVals, qrySeq_d, refIdxArr_d, shiftArr_d, posArr_d, qryLen,
				seedLength, maxNumHits, randNum);
		PRINT_CUDA_ERROR()

		/* Copy results from GPU to CPU. */
		char *refIdxArr = (char *) malloc(numQrs * maxNumHits * sizeof(char));
		cudaMemcpy(refIdxArr, refIdxArr_d, (numQrs * maxNumHits * sizeof(char)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *shiftArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(shiftArr, shiftArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *posArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(posArr, posArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		/* Check results. */
		if (((refIdxArr[0] != 0 || shiftArr[0] != 0 || posArr[0] != 0)
				&& (refIdxArr[1] != -1 || shiftArr[1] != -1 || posArr[1] != -1)
				&& (refIdxArr[2] != -1 || shiftArr[2] != -1 || posArr[2] != -1))
				|| ((refIdxArr[3] != 0 || shiftArr[3] != -2 || posArr[3] != 0)
						&& (refIdxArr[4] != -1 || shiftArr[4] != -1
								|| posArr[4] != -1)
						&& (refIdxArr[5] != -1 || shiftArr[5] != -1
								|| posArr[5] != -1))
				|| ((refIdxArr[6] != 0 || shiftArr[6] != 6 || posArr[6] != 8)
						&& (refIdxArr[7] != -1 || shiftArr[7] != -1
								|| posArr[7] != -1)
						&& (refIdxArr[8] != -1 || shiftArr[8] != -1
								|| posArr[8] != -1))
				|| ((refIdxArr[9] != -1 || shiftArr[9] != -1 || posArr[9] != -1)
						&& (refIdxArr[10] != -1 || shiftArr[10] != -1
								|| posArr[10] != -1)
						&& (refIdxArr[11] != -1 || shiftArr[11] != -1
								|| posArr[11] != -1)))
			fail("Incorrect behavior when there is 1 reference and 4 queries "
					"and max number of hits is 3.\n");

		remove(refFile);
		refDeleteHashTable(hashTable, numKeys);
		free(keys);
		free(vals);
		free(qrySeq);
		free(refIdxArr);
		free(shiftArr);
		free(posArr);
		cudaFree(keys_d);
		cudaFree(vals_d);
		cudaFree(qrySeq_d);
		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* 2 references and 4 queries. */
	{
		(void) signal(SIGSEGV, signalHandler);

		uint seedLength = 8;
		int numQrs = 4;
		int qryLen = 10;
		uint maxNumHits = 1;
		uint tupleIgnoreThres = 5;

		/* Preprocess reference file and create a hash table. */
		char refFile[MAX_FILE_NAME_LENGTH];
		sprintf(refFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "test_ref.fa");
		FILE *filePtr = fopen(refFile, "w");
		fputs(">r1\nACGTACGTAATTCCCA\n>r2\nAATTCCCAACGTACGT", filePtr);
		fclose(filePtr);
		RefPosList **hashTable = refPreprocess(refFile, seedLength);

		/* Serialize the hash table. */
		int numKeys = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength);
		int i;
		uint numVals = 0;
		RefPosList *tmpRefPosList;
		for (i = 0; i < numKeys; ++i)
		{
			tmpRefPosList = hashTable[i];
			while (tmpRefPosList != NULL)
			{
				++numVals;
				tmpRefPosList = tmpRefPosList->next;
			}
		}
		uint *keys = (uint *) malloc(numKeys * sizeof(uint));
		uint *vals = (uint *) malloc(numVals * sizeof(uint));
		serializeHashTable2_wrap(hashTable, numKeys, numVals, keys, vals,
				seedLength, tupleIgnoreThres);

		/* Copy serialized hash table to GPU. */
		uint *keys_d;
		cudaMalloc(&keys_d, numKeys * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(keys_d, keys, numKeys * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()
		uint *vals_d;
		cudaMalloc(&vals_d, numVals * sizeof(uint));
		PRINT_CUDA_ERROR()
		cudaMemcpy(vals_d, vals, numVals * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Copy queries to GPU. */
		int totalQrySeqLen = (numQrs * qryLen * sizeof(char)) + sizeof(char);
		char *qrySeq = (char *) malloc(totalQrySeqLen);
		strcpy(qrySeq, "ACGTACGTGG");
		strcat(qrySeq, "CCACGTACGT");
		strcat(qrySeq, "ACAATTCCCA");
		strcat(qrySeq, "TGGGAATTGT");
		char *qrySeq_d;
		cudaMalloc(&qrySeq_d, totalQrySeqLen);
		PRINT_CUDA_ERROR()
		cudaMemcpy(qrySeq_d, qrySeq, totalQrySeqLen, cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		/* Allocate GPU memory for results. */
		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, numQrs * maxNumHits * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, numQrs * maxNumHits * sizeof(int));
		PRINT_CUDA_ERROR()

		/* Search queries. */
		dim3 dimGrid, dimBlock;
		dimBlock.x = 1;
		dimBlock.y = 1;
		dimBlock.z = 1;
		dimGrid.x = numQrs;
		dimGrid.y = 1;
		int randNum = rand();
		refSearchQuery2_kernel<<<dimGrid, dimBlock>>>(keys_d, vals_d, numKeys,
				numVals, qrySeq_d, refIdxArr_d, shiftArr_d, posArr_d, qryLen,
				seedLength, maxNumHits, randNum);
		PRINT_CUDA_ERROR()

		/* Copy results from GPU to CPU. */
		char *refIdxArr = (char *) malloc(numQrs * maxNumHits * sizeof(char));
		cudaMemcpy(refIdxArr, refIdxArr_d, (numQrs * maxNumHits * sizeof(char)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *shiftArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(shiftArr, shiftArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int *posArr = (int *) malloc(numQrs * maxNumHits * sizeof(int));
		cudaMemcpy(posArr, posArr_d, (numQrs * maxNumHits * sizeof(int)),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		/* Check results. */
		if (((refIdxArr[0] != 0 || shiftArr[0] != 0 || posArr[0] != 0)
				&& (refIdxArr[0] != 1 || shiftArr[0] != 8 || posArr[0] != 8))
				|| ((refIdxArr[1] != 0 || shiftArr[1] != -2 || posArr[1] != 0)
						&& (refIdxArr[1] != 1 || shiftArr[1] != 6
								|| posArr[1] != 8))
				|| ((refIdxArr[2] != 0 || shiftArr[2] != 6 || posArr[2] != 8)
						&& (refIdxArr[2] != 1 || shiftArr[2] != -2
								|| posArr[2] != 0))
				|| (refIdxArr[3] != -1 || shiftArr[3] != -1 || posArr[3] != -1))
			fail("Incorrect behavior when there are 2 references and 4 queries."
					"\n");

		remove(refFile);
		refDeleteHashTable(hashTable, numKeys);
		free(keys);
		free(vals);
		free(qrySeq);
		free(refIdxArr);
		free(shiftArr);
		free(posArr);
		cudaFree(keys_d);
		cudaFree(vals_d);
		cudaFree(qrySeq_d);
		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}
}
END_TEST


/**
 * Tests @a bubbleSort2 function.
 */
START_TEST(bubbleSort2)
{
	/* Arrays have 11 elements (case 1). */
	{
		char refIdxArr[] = {8, 22, 11, 19, 22, 4, 1, 22, 10, 8, 8};
		int shiftArr[] = {5, -2, 10, 3, 0, 33, 57, 12, 3, 41, -60};
		int posArr[] = {100, 3, 79, 17, 26, 1005, 89, 52, 91, 563, 79};
		int arrSize = 11;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		bubbleSort2_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr[0] != 1 || shiftArr[0] != 57 || posArr[0] != 89
				|| refIdxArr[1] != 4 || shiftArr[1] != 33 || posArr[1] != 1005
				|| refIdxArr[2] != 8 || shiftArr[2] != -60 || posArr[2] != 79
				|| refIdxArr[3] != 8 || shiftArr[3] != 5 || posArr[3] != 100
				|| refIdxArr[4] != 8 || shiftArr[4] != 41 || posArr[4] != 563
				|| refIdxArr[5] != 10 || shiftArr[5] != 3 || posArr[5] != 91
				|| refIdxArr[6] != 11 || shiftArr[6] != 10 || posArr[6] != 79
				|| refIdxArr[7] != 19 || shiftArr[7] != 3 || posArr[7] != 17
				|| refIdxArr[8] != 22 || shiftArr[8] != -2 || posArr[8] != 3
				|| refIdxArr[9] != 22 || shiftArr[9] != 0 || posArr[9] != 26
				|| refIdxArr[10] != 22 || shiftArr[10] != 12 || posArr[10] != 52)
			fail("Incorrect behavior when arrays have 11 elements (case 1).\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* Arrays have 11 elements (case 2). */
	{
		char refIdxArr[] = {8, 22, 11, 19, 22, 4, 1, 22, 10, 8, 8};
		int shiftArr[] = {-5, -2, 10, 3, 0, 33, 57, 12, 3, -5, 5};
		int posArr[] = {100, 3, 79, 17, 26, 1005, 89, 52, 91, 99, 79};
		int arrSize = 11;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		bubbleSort2_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr[0] != 1 || shiftArr[0] != 57 || posArr[0] != 89
				|| refIdxArr[1] != 4 || shiftArr[1] != 33 || posArr[1] != 1005
				|| refIdxArr[2] != 8 || shiftArr[2] != -5 || posArr[2] != 99
				|| refIdxArr[3] != 8 || shiftArr[3] != -5 || posArr[3] != 100
				|| refIdxArr[4] != 8 || shiftArr[4] != 5 || posArr[4] != 79
				|| refIdxArr[5] != 10 || shiftArr[5] != 3 || posArr[5] != 91
				|| refIdxArr[6] != 11 || shiftArr[6] != 10 || posArr[6] != 79
				|| refIdxArr[7] != 19 || shiftArr[7] != 3 || posArr[7] != 17
				|| refIdxArr[8] != 22 || shiftArr[8] != -2 || posArr[8] != 3
				|| refIdxArr[9] != 22 || shiftArr[9] != 0 || posArr[9] != 26
				|| refIdxArr[10] != 22 || shiftArr[10] != 12 || posArr[10] != 52)
			fail("Incorrect behavior when array has 11 elements (case 2).\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}
}
END_TEST


/**
 * Tests @a quickSort2_gpu function.
 */
START_TEST(quickSort2_gpu)
{
	/* Arrays have 11 elements (case 1). */
	{
		char refIdxArr[] = {8, 22, 11, 19, 22, 4, 1, 22, 10, 8, 8};
		int shiftArr[] = {5, -2, 10, 3, 0, 33, 57, 12, 3, 41, -60};
		int posArr[] = {100, 3, 79, 17, 26, 1005, 89, 52, 91, 563, 79};
		int arrSize = 11;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		quickSort2_gpu_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr[0] != 1 || shiftArr[0] != 57 || posArr[0] != 89
				|| refIdxArr[1] != 4 || shiftArr[1] != 33 || posArr[1] != 1005
				|| refIdxArr[2] != 8 || shiftArr[2] != -60 || posArr[2] != 79
				|| refIdxArr[3] != 8 || shiftArr[3] != 5 || posArr[3] != 100
				|| refIdxArr[4] != 8 || shiftArr[4] != 41 || posArr[4] != 563
				|| refIdxArr[5] != 10 || shiftArr[5] != 3 || posArr[5] != 91
				|| refIdxArr[6] != 11 || shiftArr[6] != 10 || posArr[6] != 79
				|| refIdxArr[7] != 19 || shiftArr[7] != 3 || posArr[7] != 17
				|| refIdxArr[8] != 22 || shiftArr[8] != -2 || posArr[8] != 3
				|| refIdxArr[9] != 22 || shiftArr[9] != 0 || posArr[9] != 26
				|| refIdxArr[10] != 22 || shiftArr[10] != 12 || posArr[10] != 52)
			fail("Incorrect behavior when arrays have 11 elements (case 1).\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* Arrays have 11 elements (case 2). */
	{
		char refIdxArr[] = {8, 22, 11, 19, 22, 4, 1, 22, 10, 8, 8};
		int shiftArr[] = {-5, -2, 10, 3, 0, 33, 57, 12, 3, -5, 5};
		int posArr[] = {100, 3, 79, 17, 26, 1005, 89, 52, 91, 99, 79};
		int arrSize = 11;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		quickSort2_gpu_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr[0] != 1 || shiftArr[0] != 57 || posArr[0] != 89
				|| refIdxArr[1] != 4 || shiftArr[1] != 33 || posArr[1] != 1005
				|| refIdxArr[2] != 8 || shiftArr[2] != -5 || posArr[2] != 99
				|| refIdxArr[3] != 8 || shiftArr[3] != -5 || posArr[3] != 100
				|| refIdxArr[4] != 8 || shiftArr[4] != 5 || posArr[4] != 79
				|| refIdxArr[5] != 10 || shiftArr[5] != 3 || posArr[5] != 91
				|| refIdxArr[6] != 11 || shiftArr[6] != 10 || posArr[6] != 79
				|| refIdxArr[7] != 19 || shiftArr[7] != 3 || posArr[7] != 17
				|| refIdxArr[8] != 22 || shiftArr[8] != -2 || posArr[8] != 3
				|| refIdxArr[9] != 22 || shiftArr[9] != 0 || posArr[9] != 26
				|| refIdxArr[10] != 22 || shiftArr[10] != 12 || posArr[10] != 52)
			fail("Incorrect behavior when array has 11 elements (case 2).\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* Arrays have 1 element. */
	{
		char refIdxArr[] = {1};
		int shiftArr[] = {-5};
		int posArr[] = {100};
		int arrSize = 1;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		bubbleSort2_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr[0] != 1 || shiftArr[0] != -5 || posArr[0] != 100)
			fail("Incorrect behavior when array has 1 element.\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* Arrays have no elements. */
	{
		(void) signal(SIGSEGV, signalHandler);

		char refIdxArr[] = {};
		int shiftArr[] = {};
		int posArr[] = {};
		int arrSize = 0;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		bubbleSort2_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* Arrays have 13 elements, but only 11 should be sorted. */
	{
		char refIdxArr[] = {8, 22, 11, 19, 22, 4, 1, 22, 10, 8, 8, 2, 5};
		int shiftArr[] = {-5, -2, 10, 3, 0, 33, 57, 12, 3, -5, 5, 3, 8};
		int posArr[] = {100, 3, 79, 17, 26, 1005, 89, 52, 91, 99, 79, 10, 12};
		int arrSize = 11;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		quickSort2_gpu_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr[0] != 1 || shiftArr[0] != 57 || posArr[0] != 89
				|| refIdxArr[1] != 4 || shiftArr[1] != 33 || posArr[1] != 1005
				|| refIdxArr[2] != 8 || shiftArr[2] != -5 || posArr[2] != 99
				|| refIdxArr[3] != 8 || shiftArr[3] != -5 || posArr[3] != 100
				|| refIdxArr[4] != 8 || shiftArr[4] != 5 || posArr[4] != 79
				|| refIdxArr[5] != 10 || shiftArr[5] != 3 || posArr[5] != 91
				|| refIdxArr[6] != 11 || shiftArr[6] != 10 || posArr[6] != 79
				|| refIdxArr[7] != 19 || shiftArr[7] != 3 || posArr[7] != 17
				|| refIdxArr[8] != 22 || shiftArr[8] != -2 || posArr[8] != 3
				|| refIdxArr[9] != 22 || shiftArr[9] != 0 || posArr[9] != 26
				|| refIdxArr[10] != 22 || shiftArr[10] != 12 || posArr[10] != 52
				|| refIdxArr[11] != 2 || shiftArr[11] != 3 || posArr[11] != 10
				|| refIdxArr[12] != 5 || shiftArr[12] != 8 || posArr[12] != 12)
			fail("Incorrect behavior when array has 13 elements, but "
					"only 11 should be sorted.\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}
}
END_TEST


void getBestMatches2_wrap(char *refIdxArr, int *shiftArr,
		int *posArr, uint size, uint qryLen, int maxNumMatches,
		char *refIdxArr2, int *shiftArr2, int *posArr2);

/**
 * Tests the @a getBestMatches2 function.
 */
START_TEST(getBestMatches2)
{
	/* There is 1 max cluster and 1 max number of matches to be returned
	 * and array size is 1. */
	{
		char refIdxArr[] = {0};
		int shiftArr[] = {0};
		int posArr[] = {0};
		int arrSize = 1;
		int maxNumMatches = 1;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		char *refIdxArr2_d;
		cudaMalloc(&refIdxArr2_d, maxNumMatches * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr2_d;
		cudaMalloc(&shiftArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr2_d;
		cudaMalloc(&posArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		getBestMatches2_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize, maxNumMatches, refIdxArr2_d, shiftArr2_d, posArr2_d);
		PRINT_CUDA_ERROR()

		char refIdxArr2[maxNumMatches];
		cudaMemcpy(refIdxArr2, refIdxArr2_d, maxNumMatches * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int shiftArr2[maxNumMatches];
		cudaMemcpy(shiftArr2, shiftArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int posArr2[maxNumMatches];
		cudaMemcpy(posArr2, posArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr2[0] != 0 || shiftArr2[0] != 0 || posArr2[0] != 0)
			fail("Incorrect behavior when there is 1 max cluster and "
					"max number of matches to be returned is 1 and "
					"array size is 1.\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
		cudaFree(refIdxArr2_d);
		cudaFree(shiftArr2_d);
		cudaFree(posArr2_d);
	}

	/* There is 1 max cluster and 1 max number of matches to be returned
	 * and array size is 2. */
	{
		char refIdxArr[] = {0, 1};
		int shiftArr[] = {0, 8};
		int posArr[] = {0, 8};
		int arrSize = 2;
		int maxNumMatches = 1;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		char *refIdxArr2_d;
		cudaMalloc(&refIdxArr2_d, maxNumMatches * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr2_d;
		cudaMalloc(&shiftArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr2_d;
		cudaMalloc(&posArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		getBestMatches2_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize, maxNumMatches, refIdxArr2_d, shiftArr2_d, posArr2_d);
		PRINT_CUDA_ERROR()

		char refIdxArr2[maxNumMatches];
		cudaMemcpy(refIdxArr2, refIdxArr2_d, maxNumMatches * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int shiftArr2[maxNumMatches];
		cudaMemcpy(shiftArr2, shiftArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int posArr2[maxNumMatches];
		cudaMemcpy(posArr2, posArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if ((refIdxArr2[0] != 0 || shiftArr2[0] != 0 || posArr2[0] != 0)
				&& (refIdxArr2[0] != 1 || shiftArr2[0] != 8 || posArr2[0] != 8))
			fail("Incorrect behavior when there is 1 max cluster and "
					"max number of matches to be returned is 1 and "
					"array size is 2.\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
		cudaFree(refIdxArr2_d);
		cudaFree(shiftArr2_d);
		cudaFree(posArr2_d);
	}

	/* There is 1 max cluster and 1 max number of matches to be returned. */
	{
		char refIdxArr[] = {1, 4, 8, 8, 8, 10, 10, 19, 22, 22, 22};
		int shiftArr[] = {57, 33, -60, -60, -60, 3, 3, 3, -2, 0, 0};
		int posArr[] = {89, 1005, 79, 100, 563, 79, 91, 17, 3, 26, 52};
		int arrSize = 11;
		int maxNumMatches = 1;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		char *refIdxArr2_d;
		cudaMalloc(&refIdxArr2_d, maxNumMatches * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr2_d;
		cudaMalloc(&shiftArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr2_d;
		cudaMalloc(&posArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		getBestMatches2_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize, maxNumMatches, refIdxArr2_d, shiftArr2_d, posArr2_d);
		PRINT_CUDA_ERROR()

		char refIdxArr2[maxNumMatches];
		cudaMemcpy(refIdxArr2, refIdxArr2_d, maxNumMatches * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int shiftArr2[maxNumMatches];
		cudaMemcpy(shiftArr2, shiftArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int posArr2[maxNumMatches];
		cudaMemcpy(posArr2, posArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr2[0] != 8 || shiftArr2[0] != -60 || posArr2[0] != 79)
			fail("Incorrect behavior when there is 1 max cluster and "
					"max number of matches to be returned is 1.\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
		cudaFree(refIdxArr2_d);
		cudaFree(shiftArr2_d);
		cudaFree(posArr2_d);
	}

	/* There are 2 max clusters and 1 max number of matches to be returned. */
	{
		char refIdxArr[] = {1, 4, 8, 8, 8, 10, 10, 19, 22, 22, 22};
		int shiftArr[] = {57, 33, -60, -60, -60, 3, 3, 3, 0, 0, 0};
		int posArr[] = {89, 1005, 79, 100, 563, 79, 91, 17, 3, 26, 52};
		int arrSize = 11;
		int maxNumMatches = 1;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		char *refIdxArr2_d;
		cudaMalloc(&refIdxArr2_d, maxNumMatches * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr2_d;
		cudaMalloc(&shiftArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr2_d;
		cudaMalloc(&posArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		getBestMatches2_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize, maxNumMatches, refIdxArr2_d, shiftArr2_d, posArr2_d);
		PRINT_CUDA_ERROR()

		char refIdxArr2[maxNumMatches];
		cudaMemcpy(refIdxArr2, refIdxArr2_d, maxNumMatches * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int shiftArr2[maxNumMatches];
		cudaMemcpy(shiftArr2, shiftArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int posArr2[maxNumMatches];
		cudaMemcpy(posArr2, posArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if ((refIdxArr2[0] != 8 || shiftArr2[0] != -60 || posArr2[0] != 79)
				&& (refIdxArr2[0] != 22 || shiftArr2[0] != 0 || posArr2[0] != 3))
			fail("Incorrect behavior when there are 2 max clusters and "
					"max number of matches to be returned is 1.\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
		cudaFree(refIdxArr2_d);
		cudaFree(shiftArr2_d);
		cudaFree(posArr2_d);
	}

	/* There are 2 max clusters and 2 max number of matches to be returned. */
	{
		char refIdxArr[] = {1, 4, 8, 8, 8, 10, 10, 19, 22, 22, 22};
		int shiftArr[] = {57, 33, -60, -60, -60, 3, 3, 3, 0, 0, 0};
		int posArr[] = {89, 1005, 79, 100, 563, 79, 91, 17, 3, 26, 52};
		int arrSize = 11;
		int maxNumMatches = 2;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		char *refIdxArr2_d;
		cudaMalloc(&refIdxArr2_d, maxNumMatches * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr2_d;
		cudaMalloc(&shiftArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr2_d;
		cudaMalloc(&posArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		getBestMatches2_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize, maxNumMatches, refIdxArr2_d, shiftArr2_d, posArr2_d);
		PRINT_CUDA_ERROR()

		char refIdxArr2[maxNumMatches];
		cudaMemcpy(refIdxArr2, refIdxArr2_d, maxNumMatches * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int shiftArr2[maxNumMatches];
		cudaMemcpy(shiftArr2, shiftArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int posArr2[maxNumMatches];
		cudaMemcpy(posArr2, posArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (((refIdxArr2[0] != 8 || shiftArr2[0] != -60 || posArr2[0] != 79)
				&& (refIdxArr2[0] != 22 || shiftArr2[0] != 0 || posArr2[0] != 3))
				|| ((refIdxArr2[1] != 8 || shiftArr2[1] != -60
						|| posArr2[1] != 79) && (refIdxArr2[1] != 22
								|| shiftArr2[1] != 0 || posArr2[1] != 3)))
			fail("Incorrect behavior when there are 2 max clusters and "
					"max number of matches to be returned is 2.\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
		cudaFree(refIdxArr2_d);
		cudaFree(shiftArr2_d);
		cudaFree(posArr2_d);
	}

	/* There is 1 max cluster and 2 max number of matches to be returned. */
	{
		char refIdxArr[] = {1, 4, 8, 8, 8, 10, 10, 19, 22, 22, 22};
		int shiftArr[] = {57, 33, -60, -60, -60, 3, 3, 3, -2, 0, 0};
		int posArr[] = {89, 1005, 79, 100, 563, 79, 91, 17, 3, 26, 52};
		int arrSize = 11;
		int maxNumMatches = 2;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		char *refIdxArr2_d;
		cudaMalloc(&refIdxArr2_d, maxNumMatches * sizeof(char));
		PRINT_CUDA_ERROR()

		int *shiftArr2_d;
		cudaMalloc(&shiftArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		int *posArr2_d;
		cudaMalloc(&posArr2_d, maxNumMatches * sizeof(int));
		PRINT_CUDA_ERROR()

		getBestMatches2_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize, maxNumMatches, refIdxArr2_d, shiftArr2_d, posArr2_d);
		PRINT_CUDA_ERROR()

		char refIdxArr2[maxNumMatches];
		cudaMemcpy(refIdxArr2, refIdxArr2_d, maxNumMatches * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int shiftArr2[maxNumMatches];
		cudaMemcpy(shiftArr2, shiftArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		int posArr2[maxNumMatches];
		cudaMemcpy(posArr2, posArr2_d, maxNumMatches * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if ((refIdxArr2[0] != 8 || shiftArr2[0] != -60 || posArr2[0] != 79)
				&& (refIdxArr2[0] != -1 || shiftArr2[0] != -1
						|| posArr2[0] != -1))
			fail("Incorrect behavior when there is 1 max clusters and "
					"max number of matches to be returned is 2.\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
		cudaFree(refIdxArr2_d);
		cudaFree(shiftArr2_d);
		cudaFree(posArr2_d);
	}
}
END_TEST


/**
 * Tests @a getRandNum function.
 */
START_TEST(getRandNum)
{
	/* Array size is 10. */
	{
		int arrSize = 10;
		uint arr[arrSize];

		uint *arr_d;
		cudaMalloc(&arr_d, arrSize * sizeof(uint));
		PRINT_CUDA_ERROR()

		getRandNum_wrap<<<1, 1>>>(arr_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr, arr_d, arrSize * sizeof(uint), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (arr[0] == arr[1] && arr[1] == arr[2] && arr[2] == arr[3]
				&& arr[3] == arr[4] && arr[4] == arr[5] && arr[5] == arr[6]
				&& arr[6] == arr[7] && arr[7] == arr[8] && arr[8] == arr[9]
				&& arr[9] == arr[10])
			fail("Incorrect behavior when array size is 10.\n");

		cudaFree(arr_d);
	}

	/* Array size is 2. */
	{
		int arrSize = 2;
		uint arr[arrSize];

		uint *arr_d;
		cudaMalloc(&arr_d, arrSize * sizeof(uint));
		PRINT_CUDA_ERROR()

		getRandNum_wrap<<<1, 1>>>(arr_d, arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr, arr_d, arrSize * sizeof(uint), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (arr[0] == arr[1])
			fail("Incorrect behavior when array size is 2.\n");

		cudaFree(arr_d);
	}
}
END_TEST


/**
 * Tests @a arrGetRandomNums_gpu function.
 */
START_TEST(arrGetRandomNums)
{
	/* When both lower limit and upper limit are 0 and number of random
	 * numbers needed is 1. */
	{
		uint n = 1;
		uint arr[n];
		uint lowerLimit = 0;
		uint upperLimit = 0;

		uint *arr_d;
		cudaMalloc(&arr_d, n * sizeof(uint));
		PRINT_CUDA_ERROR()

		arrGetRandomNums_gpu_wrap<<<1, 1>>>(n, lowerLimit, upperLimit, arr_d);
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr, arr_d, n * sizeof(uint), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (arr[0] != 0)
			fail("Incorrect behavior when both lower limit and upper limit "
					"are 0 and number of random numbers needed is 1.\n");

		cudaFree(arr_d);
	}

	/* When both lower limit and upper limit are 1 and number of random
	 * numbers needed is 1. */
	{
		uint n = 1;
		uint arr[n];
		uint lowerLimit = 1;
		uint upperLimit = 1;

		uint *arr_d;
		cudaMalloc(&arr_d, n * sizeof(uint));
		PRINT_CUDA_ERROR()

		arrGetRandomNums_gpu_wrap<<<1, 1>>>(n, lowerLimit, upperLimit, arr_d);
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr, arr_d, n * sizeof(uint), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (arr[0] != 1)
			fail("Incorrect behavior when both lower limit and upper limit "
					"are 0 and number of random numbers needed is 1.\n");

		cudaFree(arr_d);
	}

	/* When lower limit is 0 and upper limit is 1 and number of random
	 * numbers needed is 1. */
	{
		uint n = 1;
		uint arr[n];
		uint lowerLimit = 0;
		uint upperLimit = 1;

		uint *arr_d;
		cudaMalloc(&arr_d, n * sizeof(uint));
		PRINT_CUDA_ERROR()

		arrGetRandomNums_gpu_wrap<<<1, 1>>>(n, lowerLimit, upperLimit, arr_d);
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr, arr_d, n * sizeof(uint), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (arr[0] != 0 && arr[0] != 1)
			fail("Incorrect behavior when lower limit is 0 and upper limit "
					"is 1 and number of random numbers needed is 1.\n");

		cudaFree(arr_d);
	}

	/* When lower limit = 1, upper limit = 20, and number of random
	 * numbers needed is 5. */
	{
		uint n = 5;
		uint arr[n];
		uint lowerLimit = 1;
		uint upperLimit = 20;

		uint *arr_d;
		cudaMalloc(&arr_d, n * sizeof(uint));
		PRINT_CUDA_ERROR()

		arrGetRandomNums_gpu_wrap<<<1, 1>>>(n, lowerLimit, upperLimit, arr_d);
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr, arr_d, n * sizeof(uint), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (arr[0] < 1 || arr[0] > 20 || arr[1] < 1 || arr[1] > 20
				|| arr[2] < 1 || arr[2] > 20 || arr[3] < 1 || arr[3] > 20
				|| arr[4] < 1 || arr[4] > 20
				|| (arr[0] == arr[1] && arr[1] == arr[2] && arr[2] == arr[3]
				   && arr[3] == arr[4])
				|| (arr[0] < arr[1] && arr[1] < arr[2] && arr[2] < arr[3]
				   && arr[3] < arr[4])
				|| (arr[0] > arr[1] && arr[1] > arr[2] && arr[2] > arr[3]
				   && arr[3] > arr[4]))
			fail("Incorrect behavior when lower limit = 1, upper limit = 20, "
					"and number of random numbers needed is 5.\n");

		cudaFree(arr_d);
	}
}
END_TEST


/**
 * Tests @ arrSearch_gpu_wrap function.
 */
START_TEST(arrSearch_gpu_wrap)
{
	/* When number to be searched is somewhere in the middle of the array. */
	{
		int arrSize = 10;

		uint arr[] = {3, 6, 1, 27, 83, 5, 8, 31, 99, 64};

		uint *arr_d;
		cudaMalloc(&arr_d, arrSize * sizeof(uint));
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr_d, arr, arrSize * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *index_d;
		cudaMalloc(&index_d, sizeof(int));
		PRINT_CUDA_ERROR()

		uint num = arr[6];
		arrSearch_gpu_wrap<<<1, 1>>>(arr_d, arrSize, num, index_d);
		PRINT_CUDA_ERROR()

		int index;
		cudaMemcpy(&index, index_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (index != 6)
			fail("Incorrect behavior when the number to be searched "
					"is somewhere in the middle of the array.\n");

		cudaFree(arr_d);
		cudaFree(index_d);
	}

	/* When number to be searched is in the beginning of the array. */
	{
		int arrSize = 10;

		uint arr[] = {3, 6, 1, 27, 83, 5, 8, 31, 99, 64};

		uint *arr_d;
		cudaMalloc(&arr_d, arrSize * sizeof(uint));
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr_d, arr, arrSize * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *index_d;
		cudaMalloc(&index_d, sizeof(int));
		PRINT_CUDA_ERROR()

		uint num = arr[0];
		arrSearch_gpu_wrap<<<1, 1>>>(arr_d, arrSize, num, index_d);
		PRINT_CUDA_ERROR()

		int index;
		cudaMemcpy(&index, index_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (index != 0)
			fail("Incorrect behavior when the number to be searched "
					"is in the beginning of the array.\n");

		cudaFree(arr_d);
		cudaFree(index_d);
	}

	/* When number to be searched is in the end of the array. */
	{
		int arrSize = 10;

		uint arr[] = {3, 6, 1, 27, 83, 5, 8, 31, 99, 64};

		uint *arr_d;
		cudaMalloc(&arr_d, arrSize * sizeof(uint));
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr_d, arr, arrSize * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *index_d;
		cudaMalloc(&index_d, sizeof(int));
		PRINT_CUDA_ERROR()

		uint num = arr[9];
		arrSearch_gpu_wrap<<<1, 1>>>(arr_d, arrSize, num, index_d);
		PRINT_CUDA_ERROR()

		int index;
		cudaMemcpy(&index, index_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (index != 9)
			fail("Incorrect behavior when the number to be searched "
					"is in the end of the array.\n");

		cudaFree(arr_d);
		cudaFree(index_d);
	}

	/* When the number to be searched is present in more than 1 place. */
	{
		int arrSize = 10;

		uint arr[] = {3, 6, 1, 27, 83, 5, 8, 1, 99, 64};

		uint *arr_d;
		cudaMalloc(&arr_d, arrSize * sizeof(uint));
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr_d, arr, arrSize * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *index_d;
		cudaMalloc(&index_d, sizeof(int));
		PRINT_CUDA_ERROR()

		uint num = arr[2];
		arrSearch_gpu_wrap<<<1, 1>>>(arr_d, arrSize, num, index_d);
		PRINT_CUDA_ERROR()

		int index;
		cudaMemcpy(&index, index_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (index != 2)
			fail("Incorrect behavior when the number to be searched "
					"is present in more than 1 place.\n");

		cudaFree(arr_d);
		cudaFree(index_d);
	}

	/* When the array has only 2 numbers. */
	{
		int arrSize = 2;

		uint arr[] = {3, 5};

		uint *arr_d;
		cudaMalloc(&arr_d, arrSize * sizeof(uint));
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr_d, arr, arrSize * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *index_d;
		cudaMalloc(&index_d, sizeof(int));
		PRINT_CUDA_ERROR()

		uint num = arr[1];
		arrSearch_gpu_wrap<<<1, 1>>>(arr_d, arrSize, num, index_d);
		PRINT_CUDA_ERROR()

		int index;
		cudaMemcpy(&index, index_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (index != 1)
			fail("Incorrect behavior when the array has only 2 element.\n");

		cudaFree(arr_d);
		cudaFree(index_d);
	}

	/* When the array has only 1 number. */
	{
		int arrSize = 1;

		uint arr[] = {3};

		uint *arr_d;
		cudaMalloc(&arr_d, arrSize * sizeof(uint));
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr_d, arr, arrSize * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *index_d;
		cudaMalloc(&index_d, sizeof(int));
		PRINT_CUDA_ERROR()

		uint num = arr[0];
		arrSearch_gpu_wrap<<<1, 1>>>(arr_d, arrSize, num, index_d);
		PRINT_CUDA_ERROR()

		int index;
		cudaMemcpy(&index, index_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (index != 0)
			fail("Incorrect behavior when the array has only 1 element.\n");

		cudaFree(arr_d);
		cudaFree(index_d);
	}

	/* When the array to be searched is empty. */
	{
		int arrSize = 0;

		uint arr[] = {};

		uint *arr_d;
		cudaMalloc(&arr_d, arrSize * sizeof(uint));
		PRINT_CUDA_ERROR()

		cudaMemcpy(arr_d, arr, arrSize * sizeof(uint), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *index_d;
		cudaMalloc(&index_d, sizeof(int));
		PRINT_CUDA_ERROR()

		uint num = arr[0];
		arrSearch_gpu_wrap<<<1, 1>>>(arr_d, arrSize, num, index_d);
		PRINT_CUDA_ERROR()

		int index;
		cudaMemcpy(&index, index_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (index != -1)
			fail("Incorrect behavior when the array is empty.\n");

		cudaFree(arr_d);
		cudaFree(index_d);
	}
}
END_TEST


/**
 * Tests @a getHash_gpu2 function.
 */
START_TEST(getHash_gpu)
{
	/* String = "ACGT". */
	{
		char s[] = "ACGT";
		int length = strlen(s);
		int hash;

		char *s_d;
		cudaMalloc(&s_d, length * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(s_d, s, length * sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *hash_d;
		cudaMalloc(&hash_d, sizeof(int));
		PRINT_CUDA_ERROR()

		getHash_gpu2_wrap<<<1, 1>>>(s_d, length, hash_d);
		PRINT_CUDA_ERROR()

		cudaMemcpy(&hash, hash_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (hash != 228)
			fail("Incorrect behavior when string is %s\n", s);

		cudaFree(s_d);
		cudaFree(hash_d);
	}

	/* String = "AAA". */
	{
		char s[] = "AAA";
		int length = strlen(s);

		char *s_d;
		cudaMalloc(&s_d, length * sizeof(char));
		PRINT_CUDA_ERROR()

		cudaMemcpy(s_d, s, length * sizeof(char), cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *hash_d;
		cudaMalloc(&hash_d, sizeof(int));
		PRINT_CUDA_ERROR()

		getHash_gpu2_wrap<<<1, 1>>>(s_d, length, hash_d);
		PRINT_CUDA_ERROR()

		int hash;
		cudaMemcpy(&hash, hash_d, sizeof(int), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (hash != 0)
			fail("Incorrect behavior when string is %s\n", s);

		cudaFree(s_d);
		cudaFree(hash_d);
	}
}
END_TEST


/**
 * Tests @a strncpy_gpu function.
 */
START_TEST(strncpy_gpu)
{
	/* Length to be copied is equal to the length of the source string. */
	{
		char s[] = "hello world";
		int n = strlen(s);
		char t[n + 1];

		char *s_d;
		cudaMalloc(&s_d, (n + 1) * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(s_d, s, (n + 1) * sizeof(char), cudaMemcpyHostToDevice);

		char *t_d;
		cudaMalloc(&t_d, (n + 1) * sizeof(char));
		PRINT_CUDA_ERROR()

		strncpy_gpu_wrap<<<1, 1>>>(t_d, s_d, n);
		PRINT_CUDA_ERROR()

		cudaMemcpy(t, t_d, (n + 1) * sizeof(char), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (strcmp(t, s) != 0)
			fail("Incorrect behavior when length to be copied is equal to "
					"the length of the source string.\n");

		cudaFree(t_d);
		cudaFree(s_d);
	}

	/* Length to be copied is less than the length of the source string and
	 * target string length is equal to source string length. */
	{
		char s[] = "hello world";
		int sLen = strlen(s);
		int n = 7;
		char t[sLen];

		char *s_d;
		cudaMalloc(&s_d, (sLen + 1) * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(s_d, s, (sLen + 1) * sizeof(char), cudaMemcpyHostToDevice);

		char *t_d;
		cudaMalloc(&t_d, sLen * sizeof(char));
		PRINT_CUDA_ERROR()

		strncpy_gpu_wrap<<<1, 1>>>(t_d, s_d, n);
		PRINT_CUDA_ERROR()

		cudaMemcpy(t, t_d, sLen * sizeof(char), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (strcmp(t, "hello w") != 0)
			fail("Incorrect behavior when length to be copied is less than "
					"the length of the source string and target string length "
					"is equal to source string length.\n");

		cudaFree(t_d);
		cudaFree(s_d);
	}

	/* Length to be copied is equal to the length of the source string and
	 * target string length is greater than the source string length. */
	{
		char s[] = "hello world";
		int sLen = strlen(s);
		int n = sLen;
		char t[2 * sLen];

		char *s_d;
		cudaMalloc(&s_d, (sLen + 1) * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(s_d, s, (sLen + 1) * sizeof(char), cudaMemcpyHostToDevice);

		char *t_d;
		cudaMalloc(&t_d, (2 * sLen) * sizeof(char));
		PRINT_CUDA_ERROR()

		strncpy_gpu_wrap<<<1, 1>>>(t_d, s_d, n);
		PRINT_CUDA_ERROR()

		cudaMemcpy(t, t_d, (2 * sLen) * sizeof(char), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (strcmp(t, "hello world") != 0)
			fail("Incorrect behavior when length to be copied is equal to "
					"the length of the source string and target string length "
					"is greater than source string length.\n");

		cudaFree(t_d);
		cudaFree(s_d);
	}

	/* Length to be copied is less than the length of the source string and
	 * target string length is greater than the source string length. */
	{
		char s[] = "hello world";
		int sLen = strlen(s);
		int n = 7;
		char t[2 * sLen];

		char *s_d;
		cudaMalloc(&s_d, (sLen + 1) * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(s_d, s, (sLen + 1) * sizeof(char), cudaMemcpyHostToDevice);

		char *t_d;
		cudaMalloc(&t_d, (sLen + 1) * sizeof(char));
		PRINT_CUDA_ERROR()

		strncpy_gpu_wrap<<<1, 1>>>(t_d, s_d, n);
		PRINT_CUDA_ERROR()

		cudaMemcpy(t, t_d, (2 * sLen) * sizeof(char), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (strcmp(t, "hello w") != 0)
			fail("Incorrect behavior when length to be copied is less than "
					"the length of the source string and target string length "
					"is greater than source string length.\n");

		cudaFree(t_d);
		cudaFree(s_d);
	}

	/* Length to be copied is less than the length of the source string and
	 * target string length is less than the source string length. */
	{
		char s[] = "hello world";
		int sLen = strlen(s);
		int n = 7;
		char t[n + 1];

		char *s_d;
		cudaMalloc(&s_d, (sLen + 1) * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(s_d, s, (sLen + 1) * sizeof(char), cudaMemcpyHostToDevice);

		char *t_d;
		cudaMalloc(&t_d, (n + 1) * sizeof(char));
		PRINT_CUDA_ERROR()

		strncpy_gpu_wrap<<<1, 1>>>(t_d, s_d, n);
		PRINT_CUDA_ERROR()

		cudaMemcpy(t, t_d, (n + 1) * sizeof(char), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (strcmp(t, "hello w") != 0)
			fail("Incorrect behavior when length to be copied is less than "
					"the length of the source string and target string length "
					"is less than source string length.\n");

		cudaFree(t_d);
		cudaFree(s_d);
	}

	/* Length to be copied is greater than the length of the source string and
	 * target string length is greater than the source string length. */
	{
		char s[] = "hello world";
		int sLen = strlen(s);
		int n = sLen + 10;
		char t[n + 1];

		char *s_d;
		cudaMalloc(&s_d, (sLen + 1) * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(s_d, s, (sLen + 1) * sizeof(char), cudaMemcpyHostToDevice);

		char *t_d;
		cudaMalloc(&t_d, (n + 1) * sizeof(char));
		PRINT_CUDA_ERROR()

		strncpy_gpu_wrap<<<1, 1>>>(t_d, s_d, n);
		PRINT_CUDA_ERROR()

		cudaMemcpy(t, t_d, (n + 1) * sizeof(char), cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (strcmp(t, "hello world") != 0)
			fail("Incorrect behavior when length to be copied is greater than "
					"the length of the source string and target string length "
					"is greater than source string length.\n");

		cudaFree(t_d);
		cudaFree(s_d);
	}
}
END_TEST


/**
 * Tests @a insertionSort3_gpu function.
 */
START_TEST(insertionSort3_gpu)
{
	/* Arrays have 2 elements and the last element is not in order. */
	{
		char refIdxArr[] = {22, 8};
		int shiftArr[] = {5, -2};
		int posArr[] = {100, 114};
		int arrSize = 2;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		insertionSort3_gpu_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr[0] != 8 || shiftArr[0] != -2 || posArr[0] != 114
				|| refIdxArr[1] != 22 || shiftArr[1] != 5 || posArr[1] != 100)
			fail("Incorrect behavior when arrays have 2 elements and the last "
					"element is not in order.\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* Arrays have 2 elements and the last element is in order. */
	{
		char refIdxArr[] = {8, 22};
		int shiftArr[] = {5, -2};
		int posArr[] = {100, 114};
		int arrSize = 2;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		insertionSort3_gpu_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr[0] != 8 || shiftArr[0] != 5 || posArr[0] != 100
				|| refIdxArr[1] != 22 || shiftArr[1] != -2 || posArr[1] != 114)
			fail("Incorrect behavior when arrays have 2 elements and the last "
					"element is in order.\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* Arrays have 5 elements and the last element is supposed to be
	 * in the beginning of the array. */
	{
		char refIdxArr[] = {8, 22, 24, 24, 8};
		int shiftArr[] = {5, -2, 10, 12, 3};
		int posArr[] = {100, 114, 15, 42, 677};
		int arrSize = 5;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		insertionSort3_gpu_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr[0] != 8 || shiftArr[0] != 3 || posArr[0] != 677
				|| refIdxArr[1] != 8 || shiftArr[1] != 5 || posArr[1] != 100
				|| refIdxArr[2] != 22 || shiftArr[2] != -2 || posArr[2] != 114
				|| refIdxArr[3] != 24 || shiftArr[3] != 10 || posArr[3] != 15
				|| refIdxArr[4] != 24 || shiftArr[4] != 12 || posArr[4] != 42)
			fail("Incorrect behavior when arrays have 5 elements and the "
					"last element is supposed to be in the beginning of the "
					"array.\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* Arrays have 5 elements and the last element is supposed to be
	 * in the end of the array. */
	{
		char refIdxArr[] = {8, 22, 24, 24, 24};
		int shiftArr[] = {5, -2, 10, 12, 12};
		int posArr[] = {100, 114, 15, 42, 677};
		int arrSize = 5;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		insertionSort3_gpu_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr[0] != 8 || shiftArr[0] != 5 || posArr[0] != 100
				|| refIdxArr[1] != 22 || shiftArr[1] != -2 || posArr[1] != 114
				|| refIdxArr[2] != 24 || shiftArr[2] != 10 || posArr[2] != 15
				|| refIdxArr[3] != 24 || shiftArr[3] != 12 || posArr[3] != 42
				|| refIdxArr[4] != 24 || shiftArr[4] != 12 || posArr[4] != 677)
			fail("Incorrect behavior when arrays have 5 elements and the "
					"last element is supposed to be in the end of the array."
					"\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}

	/* Arrays have 5 elements and the last element is supposed to be
	 * in the middle of the array. */
	{
		char refIdxArr[] = {8, 22, 24, 24, 22};
		int shiftArr[] = {5, -2, 10, 12, 12};
		int posArr[] = {100, 114, 15, 42, 677};
		int arrSize = 5;

		char *refIdxArr_d;
		cudaMalloc(&refIdxArr_d, arrSize * sizeof(char));
		PRINT_CUDA_ERROR()
		cudaMemcpy(refIdxArr_d, refIdxArr, arrSize * sizeof(char),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *shiftArr_d;
		cudaMalloc(&shiftArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(shiftArr_d, shiftArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		int *posArr_d;
		cudaMalloc(&posArr_d, arrSize * sizeof(int));
		PRINT_CUDA_ERROR()
		cudaMemcpy(posArr_d, posArr, arrSize * sizeof(int),
				cudaMemcpyHostToDevice);
		PRINT_CUDA_ERROR()

		insertionSort3_gpu_wrap<<<1, 1>>>(refIdxArr_d, shiftArr_d, posArr_d,
				arrSize);
		PRINT_CUDA_ERROR()

		cudaMemcpy(refIdxArr, refIdxArr_d, arrSize * sizeof(char),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(shiftArr, shiftArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		cudaMemcpy(posArr, posArr_d, arrSize * sizeof(int),
				cudaMemcpyDeviceToHost);
		PRINT_CUDA_ERROR()

		if (refIdxArr[0] != 8 || shiftArr[0] != 5 || posArr[0] != 100
				|| refIdxArr[1] != 22 || shiftArr[1] != -2 || posArr[1] != 114
				|| refIdxArr[2] != 22 || shiftArr[2] != 12 || posArr[2] != 677
				|| refIdxArr[3] != 24 || shiftArr[3] != 10 || posArr[3] != 15
				|| refIdxArr[4] != 24 || shiftArr[4] != 12 || posArr[4] != 42)
			fail("Incorrect behavior when arrays have 5 elements and the "
					"last element is supposed to be in the middle of the array."
					"\n");

		cudaFree(refIdxArr_d);
		cudaFree(shiftArr_d);
		cudaFree(posArr_d);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *referenceSuite(void)
{
	Suite *s = suite_create("reference");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, searchQuery);
	tcase_add_test(testCaseCore, searchQuery_paired);
	tcase_add_test(testCaseCore, preprocess);
	tcase_add_test(testCaseCore, serializeHashTable2_wrap);
	tcase_add_test(testCaseCore, pow_gpu);
	tcase_add_test(testCaseCore, getBestMatches2);
	tcase_add_test(testCaseCore, refSearchQuery2_gpu);
	tcase_add_test(testCaseCore, bubbleSort2);
	tcase_add_test(testCaseCore, quickSort2_gpu);
	tcase_add_test(testCaseCore, getRandNum);
	tcase_add_test(testCaseCore, arrGetRandomNums);
	tcase_add_test(testCaseCore, arrSearch_gpu_wrap);
	tcase_add_test(testCaseCore, getHash_gpu);
	tcase_add_test(testCaseCore, strncpy_gpu);
	tcase_add_test(testCaseCore, insertionSort3_gpu);
	suite_add_tcase(s, testCaseCore);

	return s;
}
