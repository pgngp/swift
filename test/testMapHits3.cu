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
#include "../src/mapHits3.h"
#include "testMapHits3.h"


/**
 * Tests @a isEqual function.
 */
START_TEST(isEqual)
{
	/* First hit is less than second hit. */
	{
		char refIdx1 = 1;
		int shift1 = 10;
		int refPos1 = 20;
		BST *node1 = createBSTNode_wrap(refIdx1, shift1, refPos1);

		char refIdx2 = 2;
		int shift2 = 5;
		int refPos2 = 10;
		BST *node2 = createBSTNode_wrap(refIdx2, shift2, refPos2);

		int result = isEqual_wrap(node1, node2);
		if (result >= 0)
			fail("Incorrect behavior when first hit is less than second hit.\n");

		free(node1);
		free(node2);
	}

	/* First hit is equal to second hit. */
	{
		char refIdx1 = 1;
		int shift1 = 10;
		int refPos1 = 20;
		BST *node1 = createBSTNode_wrap(refIdx1, shift1, refPos1);

		char refIdx2 = 1;
		int shift2 = 10;
		int refPos2 = 30;
		BST *node2 = createBSTNode_wrap(refIdx2, shift2, refPos2);

		int result = isEqual_wrap(node1, node2);
		if (result != 0)
			fail("Incorrect behavior when first hit is equal to second hit.\n");

		free(node1);
		free(node2);
	}

	/* First hit is greater than second hit. */
	{
		char refIdx1 = 3;
		int shift1 = 10;
		int refPos1 = 20;
		BST *node1 = createBSTNode_wrap(refIdx1, shift1, refPos1);

		char refIdx2 = 1;
		int shift2 = 20;
		int refPos2 = 30;
		BST *node2 = createBSTNode_wrap(refIdx2, shift2, refPos2);

		int result = isEqual_wrap(node1, node2);
		if (result <= 0)
			fail("Incorrect behavior when first hit is greater than second hit.\n");

		free(node1);
		free(node2);
	}
}
END_TEST


/**
 * Tests @a findBiggestClusters function.
 */
START_TEST(findBiggestClusters)
{
	/* There is 1 biggest cluster. */
	{
		int seedLen = 8;
		mapHits3Create(seedLen);

		char refIdx = 1;
		int shift = 10;
		int refPos = 20;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 1;
		shift = 5;
		refPos = 30;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 1;
		shift = 10;
		refPos = 30;
		mapHits3AddHit(refIdx, shift, refPos);

		BST *bsTree = NULL, **biggestClusters = NULL;
		int numBiggestClusters = 0;
		findBiggestClusters_wrap(&bsTree, &biggestClusters, &numBiggestClusters);

		if (numBiggestClusters != 1 || biggestClusters[0]->refIdx != 1
				|| biggestClusters[0]->shift != 10 || biggestClusters[1] != NULL)
			fail("Incorrect behavior when there is 1 biggest cluster.\n");

		mapHits3Delete();
	}

	/* There are 2 biggest clusters. */
	{
		int seedLen = 8;
		mapHits3Create(seedLen);

		char refIdx = 1;
		int shift = 10;
		int refPos = 20;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 1;
		shift = 5;
		refPos = 30;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 1;
		shift = 10;
		refPos = 30;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 1;
		shift = 5;
		refPos = 10;
		mapHits3AddHit(refIdx, shift, refPos);

		BST *bsTree = NULL, **biggestClusters = NULL;
		int numBiggestClusters = 0;
		findBiggestClusters_wrap(&bsTree, &biggestClusters, &numBiggestClusters);

		if (biggestClusters[0]->refIdx != 1
				|| biggestClusters[0]->shift != 10
				|| biggestClusters[1]->refIdx != 1
				|| biggestClusters[1]->shift != 5
				|| biggestClusters[2] != NULL)
			fail("Incorrect behavior when there are 2 biggest clusters.\n");

		mapHits3Delete();
	}
}
END_TEST


/**
 * Tests @a search function.
 */
START_TEST(search)
{
	/* Node to be searched for is already present in the tree. */
	{
		char refIdx = 1;
		int shift = 10;
		int refPos = 20;
		BST *node1 = createBSTNode_wrap(refIdx, shift, refPos);

		refIdx = 2;
		shift = 10;
		refPos = 20;
		BST *node2 = createBSTNode_wrap(refIdx, shift, refPos);
		node1->right = node2;

		refIdx = 1;
		shift = 10;
		refPos = 15;
		BST *node3 = createBSTNode_wrap(refIdx, shift, refPos);

		BST *parent = NULL;
		BST *node4 = search_wrap(node1, node3, &parent);
		if (node4 == NULL || node4->refIdx != node3->refIdx
				|| node4->shift != node3->shift)
			fail("Incorrect behavior when the node to be searched for is "
					"already present in the tree.\n");

		free(node3);
		deleteTree_wrap(node1);
	}

	/* Node to be searched for is already present in the tree, but not
	 * at the root. */
	{
		char refIdx = 1;
		int shift = 15;
		int refPos = 20;
		BST *node1 = createBSTNode_wrap(refIdx, shift, refPos);

		refIdx = 1;
		shift = 10;
		refPos = 20;
		BST *node2 = createBSTNode_wrap(refIdx, shift, refPos);
		node1->left = node2;

		refIdx = 1;
		shift = 10;
		refPos = 15;
		BST *node3 = createBSTNode_wrap(refIdx, shift, refPos);

		BST *parent = NULL;
		BST *node4 = search_wrap(node1, node3, &parent);
		if (node4 == NULL)
			fprintf(stderr, "node4 is null\n");
		if (node4 == NULL || node4->refIdx != node3->refIdx
				|| node4->shift != node3->shift)
			fail("Incorrect behavior when the node to be searched for is "
					"already present in the tree, but not at the root.\n");

		free(node3);
		deleteTree_wrap(node1);
	}

	/* Node to be searched for is not already present in the tree. */
	{
		char refIdx = 1;
		int shift = 10;
		int refPos = 20;
		BST *node1 = createBSTNode_wrap(refIdx, shift, refPos);

		refIdx = 2;
		shift = 10;
		refPos = 20;
		BST *node2 = createBSTNode_wrap(refIdx, shift, refPos);
		node1->right = node2;

		refIdx = 1;
		shift = 20;
		refPos = 15;
		BST *node3 = createBSTNode_wrap(refIdx, shift, refPos);

		BST *parent = NULL;
		BST *node4 = search_wrap(node1, node3, &parent);
		if (node4 != NULL || parent->refIdx != node2->refIdx
				|| parent->shift != node2->shift)
			fail("Incorrect behavior when the node to be searched for is "
					"not already present in the tree.\n");

		free(node3);
		deleteTree_wrap(node1);
	}
}
END_TEST


/**
 * Tests @a mapHits3AddHit function.
 */
START_TEST(mapHits3AddHit)
{
	/* Test 1. */
	{
		int seedLen = 8;
		mapHits3Create(seedLen);

		char refIdx = 1;
		int shift = 10;
		int refPos = 20;
		int biggestClusterSize = 0;
		BST *tree = mapHits3AddHit_wrap(refIdx, shift, refPos,
				&biggestClusterSize);
		BST *node1 = createBSTNode_wrap(refIdx, shift, refPos);
		if (biggestClusterSize != 1 || isEqual_wrap(node1, tree) != 0)
			fail("Incorrect behavior test 1 (case 1).\n");

		refIdx = 2;
		shift = 10;
		refPos = 20;
		tree = mapHits3AddHit_wrap(refIdx, shift, refPos, &biggestClusterSize);
		BST *node2 = createBSTNode_wrap(refIdx, shift, refPos);
		if (biggestClusterSize != 1 || isEqual_wrap(tree->right, node2) != 0)
			fail("Incorrect behavior in test 1 (case 2).\n");

		refIdx = 1;
		shift = 10;
		refPos = 30;
		tree = mapHits3AddHit_wrap(refIdx, shift, refPos, &biggestClusterSize);
		if (biggestClusterSize != 2 || tree->numHits != 2)
			fail("Incorrect behavior in test 1 (case 3).\n");

		mapHits3Delete();
		free(node1);
		free(node2);
	}

	/* Test 2. */
	{
		int seedLen = 8;
		mapHits3Create(seedLen);

		char refIdx = 2;
		int shift = 10;
		int refPos = 10;
		int biggestClusterSize = 0;
		BST *tree = mapHits3AddHit_wrap(refIdx, shift, refPos,
				&biggestClusterSize);
		BST *node1 = createBSTNode_wrap(refIdx, shift, refPos);
		if (biggestClusterSize != 1 || isEqual_wrap(node1, tree) != 0)
			fail("Incorrect behavior in test 2 (case 1).\n");

		refIdx = 1;
		shift = 10;
		refPos = 20;
		tree = mapHits3AddHit_wrap(refIdx, shift, refPos, &biggestClusterSize);
		BST *node2 = createBSTNode_wrap(refIdx, shift, refPos);
		if (biggestClusterSize != 1 || isEqual_wrap(tree->left, node2) != 0)
			fail("Incorrect behavior in test 2 (case 2).\n");

		refIdx = 1;
		shift = 10;
		refPos = 30;
		tree = mapHits3AddHit_wrap(refIdx, shift, refPos, &biggestClusterSize);
		if (biggestClusterSize != 2 || tree->left->numHits != 2)
			fail("Incorrect behavior in test 2 (case 3).\n");

		mapHits3Delete();
		free(node1);
		free(node2);
	}

}
END_TEST


/**
 * Tests @a mapHits3GetBestHits function.
 */
START_TEST(mapHits3GetBestHits)
{
	/* There is 1 best hit and max number of hits is 3. */
	{
		int seedLen = 8;
		mapHits3Create(seedLen);

		char refIdx = 2;
		int shift = 10;
		int refPos = 10;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 1;
		shift = 10;
		refPos = 20;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 1;
		shift = 10;
		refPos = 30;
		mapHits3AddHit(refIdx, shift, refPos);

		int numBestHits = 3;
		char *refIdx_bestHit = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestHit = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestHit = (int *) calloc(numBestHits, sizeof(int));
		mapHits3GetBestHits(numBestHits, refIdx_bestHit, shift_bestHit,
				refPos_bestHit);

		if (refIdx_bestHit[0] != 1 || shift_bestHit[0] != 10
				|| refPos_bestHit[0] != 20)
			fail("Incorrect behavior when there is 1 best hit and max number "
					"of hits is 3.\n");

		mapHits3Delete();
		free(refIdx_bestHit);
		free(shift_bestHit);
		free(refPos_bestHit);
	}

	/* There are 3 best hits and max number of hits is 2. */
	{
		int seedLen = 8;
		mapHits3Create(seedLen);

		char refIdx = 2;
		int shift = 10;
		int refPos = 10;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 1;
		shift = 10;
		refPos = 20;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 1;
		shift = 10;
		refPos = 30;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 3;
		shift = 20;
		refPos = 40;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 3;
		shift = 20;
		refPos = 60;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 4;
		shift = 5;
		refPos = 10;
		mapHits3AddHit(refIdx, shift, refPos);

		refIdx = 4;
		shift = 5;
		refPos = 30;
		mapHits3AddHit(refIdx, shift, refPos);

		int numBestHits = 2;
		char *refIdx_bestHit = (char *) calloc(numBestHits, sizeof(char));
		int *shift_bestHit = (int *) calloc(numBestHits, sizeof(int));
		int *refPos_bestHit = (int *) calloc(numBestHits, sizeof(int));
		mapHits3GetBestHits(numBestHits, refIdx_bestHit, shift_bestHit,
				refPos_bestHit);

		if (((refIdx_bestHit[0] == 1 && refIdx_bestHit[1] == 3)
				|| (refIdx_bestHit[0] == 1 && refIdx_bestHit[1] == 4)
				|| (refIdx_bestHit[0] == 3 && refIdx_bestHit[1] == 4))
				&& ((shift_bestHit[0] == 10 && shift_bestHit[1] == 20)
						|| (shift_bestHit[0] == 10 && shift_bestHit[1] == 5)
						|| (shift_bestHit[0] == 20 && shift_bestHit[1] == 5))
				&& ((refPos_bestHit[0] == 20 && refPos_bestHit[1] == 40)
						|| (refPos_bestHit[0] == 20 && refPos_bestHit[1] == 10)
						|| (refPos_bestHit[0] == 40 && refPos_bestHit[1] == 10)))
		{ }
		else
			fail("Incorrect behavior when there are 3 best hits and max number "
					"of hits is 2.\n");

		mapHits3Delete();
		free(refIdx_bestHit);
		free(shift_bestHit);
		free(refPos_bestHit);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *mapHits3Suite(void)
{
	Suite *s = suite_create("mapHits3");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, isEqual);
	tcase_add_test(testCaseCore, findBiggestClusters);
	tcase_add_test(testCaseCore, search);
	tcase_add_test(testCaseCore, mapHits3AddHit);
	tcase_add_test(testCaseCore, mapHits3GetBestHits);
	suite_add_tcase (s, testCaseCore);

	return s;
}
