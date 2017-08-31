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
#include "../src/hitList.h"
#include "testHitList.h"


/**
 * Tests creation of HitList node.
 */
START_TEST(createNode)
{
	HitList *node = hitListCreateNode(0, 0, 0);
	if (node == NULL || node->next != NULL || node->prev != NULL)
		fail("Incorrect behavior when a HitList node is created.\n");
}
END_TEST


/**
 * Tests hitListSwapAdjNodes function.
 */
START_TEST(swapAdjNodes)
{
	/* Test behavior when both node1 and node2 are NULL. */
	{
		HitList *node1 = NULL;
		HitList *node2 = NULL;
		HitList *head = hitListSwapAdjNodes(node1, node2);
		if (head != NULL || node1 != NULL || node2 != NULL)
			fail("Incorrect behavior when both node1 and node2 are NULL.\n");
	}

	/* Test behavior when node1 is NULL and node2 is not NULL. */
	{
		HitList *node1 = NULL;
		HitList *node2 = hitListCreateNode(0, 0, 0);
		HitList *head = hitListSwapAdjNodes(node1, node2);
		if (head != NULL || node1 != NULL || node2->index != 0
				|| node2->shift != 0 || node2->offset != 0
				|| node2->prev != NULL || node2->next != NULL)
			fail("Incorrect behavior when node1 is NULL and node2 is not "
					"NULL.\n");
		hitListDelete(node2);
	}


	/* Test behavior when node1 is not NULL and node2 is NULL. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = NULL;
		HitList *head = hitListSwapAdjNodes(node1, node2);
		if (head != NULL || node1->index != 0 || node1->shift != 0
				|| node1->offset != 0 || node1->prev != NULL
				|| node1->next != NULL || node2 != NULL)
			fail("Incorrect behavior when node1 is not NULL and node2 is "
					"NULL.\n");
		hitListDelete(node1);
	}


	/* Test behavior when node1 and node2 are ordinary adjacent nodes. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(1, 1, 1);
		node1->next = node2;
		node2->prev = node1;

		HitList *head = hitListSwapAdjNodes(node1, node2);

		if (head != node2 || head->next != node1
				|| node2->next != node1 || node2->prev != NULL
				|| node2->index != 1 || node2->shift != 1
				|| node2->offset != 1 || node1->index != 0 || node1->shift != 0
				|| node1->offset != 0 || node1->next != NULL
				|| node1->prev != node2)
			fail("Incorrect behavior when node1 and node2 are ordinary "
					"adjacent nodes.\n");
		hitListDelete(node1);
		hitListDelete(node2);
	}
}
END_TEST


/**
 * Tests @a hitListSort function.
 */
START_TEST(sort)
{
	/* Test behavior when list is empty. */
	{
		HitList *head = NULL;
		head = hitListSort(head, 0);
		if (head != NULL)
			fail("Incorrect behavior when list is empty.\n");
	}

	/* Test behavior when list size is 1. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListSort(node1, 1);
		if (node2 != node1)
			fail("Incorrect behavior when list size is 1.\n");
		hitListDelete(node1);
	}

	/* Test behavior when list size is 2. */
	{
		HitList *node1 =  hitListCreateNode(1, 1, 1);
		HitList *node2 =  hitListCreateNode(0, 0, 0);
		node1->next = node2;
		node2->prev = node1;
		HitList *head = node1;
		head = hitListSort(head, 2);
		if (head != node2 || head->next != node1)
			fail("Incorrect behavior when list size is 2.\n");
		hitListDelete(head);
	}

	/* Test behavior when some nodes have the same index. */
	{
		HitList *node1 =  hitListCreateNode(1, 1, 1);
		HitList *node2 =  hitListCreateNode(0, 0, 0);
		HitList *node3 =  hitListCreateNode(0, 2, 2);
		node1->next = node2;
		node2->next = node3;
		node3->prev = node2;
		node2->prev = node1;
		HitList *head = node1;
		head = hitListSort(head, 3);
		if (head != node2 || head->next != node3 || head->next->next != node1)
			fail("Incorrect behavior when list size is 3.\n");
		hitListDelete(head);
	}

	/* Test behavior when some nodes have the same index and same shift. */
	{
		HitList *node1 =  hitListCreateNode(1, 1, 1);
		HitList *node2 =  hitListCreateNode(0, 0, 0);
		HitList *node3 =  hitListCreateNode(0, 0, 2);
		HitList *node4 =  hitListCreateNode(1, 1, 0);
		node1->next = node2;
		node2->next = node3;
		node3->next = node4;
		node4->prev = node3;
		node3->prev = node2;
		node2->prev = node1;
		HitList *head = node1;
		head = hitListSort(head, 4);
		HitList *tmp1 = head;
		HitList *tmp2 = tmp1->next;
		HitList *tmp3 = tmp2->next;
		HitList *tmp4 = tmp3->next;
		if (tmp1 != node2 || tmp2 != node3 || tmp3 != node4 || tmp4 != node1)
			fail("Incorrect behavior when some nodes have the same index "
					"and the same shift.\n");
		hitListDelete(head);
	}

	/* Test behavior when some nodes have the same index and same shift and
	 * same offset. */
	{
		HitList *node1 =  hitListCreateNode(1, 1, 1);
		HitList *node2 =  hitListCreateNode(0, 0, 0);
		HitList *node3 =  hitListCreateNode(0, 0, 0);
		HitList *node4 =  hitListCreateNode(1, 1, 1);
		HitList *node5 =  hitListCreateNode(0, 1, 1);
		node1->next = node2;
		node2->next = node3;
		node3->next = node4;
		node4->next = node5;
		node5->prev = node4;
		node4->prev = node3;
		node3->prev = node2;
		node2->prev = node1;
		HitList *head = node1;
		head = hitListSort(head, 5);
		HitList *tmpA = head;
		HitList *tmpB = tmpA->next;
		HitList *tmpC = tmpB->next;
		HitList *tmpD = tmpC->next;
		HitList *tmpE = tmpD->next;
		if ((tmpA != node2 && tmpA != node3)
				|| (tmpB != node2 && tmpB != node3)
				|| tmpC != node5 || (tmpD != node4 && tmpD != node1)
				|| (tmpE != node4 && tmpE != node1))
			fail("Incorrect behavior when some nodes have the same index "
					"and same shift and same offset.\n");
		hitListDelete(head);
	}
}
END_TEST


/**
 * Tests the hitListCmpNodes function.
 */
START_TEST(cmpNodes)
{
	/* Test behavior when both node1 and node2 are NULL. */
	{
		int result = hitListCmpNodes(NULL, NULL);
		if (result == 0 || result == -1 || result == 1)
			fail("Incorrect behavior when both node1 and node2 are NULL.\n");
	}

	/* Test behavior when node1 is NULL. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		int result = hitListCmpNodes(node1, NULL);
		if (result == 0 || result == -1 || result == 1)
			fail("Incorrect behavior when node1 is NULL.\n");
		hitListDelete(node1);
	}

	/* Test behavior when node2 is NULL. */
	{
		HitList *node2 = hitListCreateNode(0, 0, 0);
		int result = hitListCmpNodes(NULL, node2);
		if (result == 0 || result == -1 || result == 1)
			fail("Incorrect behavior when node2 is NULL.\n");
		hitListDelete(node2);
	}

	/* Test behavior when node1 and node2 are equal. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 0, 0);
		int result1 = hitListCmpNodes(node1, node2);
		int result2 = hitListCmpNodes(node2, node1);
		if (result1 != 0 || result2 != 0)
			fail("Incorrect behavior when both nodes are equal.\n");
		hitListDelete(node1);
		hitListDelete(node2);
	}

	/* Test behavior when node1 is smaller than node2. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 1, 0);
		HitList *node3 = hitListCreateNode(0, 0, 1);
		HitList *node4 = hitListCreateNode(1, 0, 0);
		HitList *node5 = hitListCreateNode(1, 0, 0);
		int result = hitListCmpNodes(node1, node2);
		if (result != -1)
			fail("Incorrect behavior when node1 is smaller than node2.\n");
		result = hitListCmpNodes(node1, node3);
		if (result != -1)
			fail("Incorrect behavior when node1 is smaller than node2.\n");
		result = hitListCmpNodes(node1, node4);
		if (result != -1)
			fail("Incorrect behavior when node1 is smaller than node2.\n");
		result = hitListCmpNodes(node1, node5);
		if (result != -1)
			fail("Incorrect behavior when node1 is smaller than node2.\n");
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		hitListDelete(node4);
		hitListDelete(node5);
	}

	/* Test behavior when node1 is greater than node2. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 1, 0);
		HitList *node3 = hitListCreateNode(0, 0, 1);
		HitList *node4 = hitListCreateNode(1, 0, 0);
		HitList *node5 = hitListCreateNode(1, 0, 0);
		int result = hitListCmpNodes(node2, node1);
		if (result != 1)
			fail("Incorrect behavior when node1 is greater than node2.\n");
		result = hitListCmpNodes(node3, node1);
		if (result != 1)
			fail("Incorrect behavior when node1 is greater than node2.\n");
		result = hitListCmpNodes(node4, node1);
		if (result != 1)
			fail("Incorrect behavior when node1 is greater than node2.\n");
		result = hitListCmpNodes(node5, node1);
		if (result != 1)
			fail("Incorrect behavior when node1 is greater than node2.\n");
		hitListDelete(node1);
		hitListDelete(node2);
		hitListDelete(node3);
		hitListDelete(node4);
		hitListDelete(node5);
	}
}
END_TEST


/**
 * Tests hitListJoin function.
 */
START_TEST(join)
{
	/* Test behavior when both lists are NULL. */
	{
		HitList *list1 = NULL;
		HitList *list2 = NULL;
		HitList *head = hitListJoin(list1, list2);
		if (head != NULL)
			fail("Incorrect behavior when both list1 and list2 are NULL.\n");
	}

	/* Test behavior when list1 is NULL. */
	{
		HitList *list1 = NULL;
		HitList *list2 = hitListCreateNode(0, 0, 0);
		HitList *head = hitListJoin(list1, list2);
		if (head != list2 || head->prev != NULL)
			fail("Incorrect behavior when list1 is NULL.\n");
		hitListDelete(list2);
	}

	/* Test behavior when list2 is NULL. */
	{
		HitList *list1 = hitListCreateNode(0, 0, 0);
		HitList *list2 = NULL;
		HitList *head = hitListJoin(list1, list2);
		if (head != list1 || head->prev != NULL)
			fail("Incorrect behavior when list2 is NULL.\n");
		hitListDelete(list1);
	}

	/* Test behavior when both list1 and list2 have 1 node. */
	{
		HitList *list1 = hitListCreateNode(0, 0, 0);
		HitList *list2 = hitListCreateNode(0, 0, 0);
		HitList *head = hitListJoin(list1, list2);
		if (head != list1 || list1->prev != NULL || list1->next != list2
				|| list2->prev != list1 || list2->next != NULL)
			fail("Incorrect behavior when both list1 and list2 have 1 node "
					"each.\n");
		hitListDelete(head);
	}

	/* Test behavior when list1 has 1 node and list2 has more than 1 node. */
	{
		HitList *list1 = hitListCreateNode(0, 0, 0);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 0, 0);
		HitList *node3 = hitListCreateNode(0, 0, 0);
		HitList *list2 = node1;
		node1->next = node2;
		node2->next = node3;
		node3->prev = node2;
		node2->prev = node1;
		HitList *head = hitListJoin(list1, list2);
		if (head != list1 || head->next != list2 || list2->prev != list1)
			fail("Incorrect behavior when list1 has 1 node and list2 has "
					"more than 1 node.\n");
		hitListDelete(head);
	}

	/* Test behavior when list1 has more than 1 node and list2 has 1 node. */
	{
		HitList *list1 = hitListCreateNode(0, 0, 0);
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 0, 0);
		HitList *node3 = hitListCreateNode(0, 0, 0);
		HitList *list2 = node1;
		node1->next = node2;
		node2->next = node3;
		node3->prev = node2;
		node2->prev = node1;
		HitList *head = hitListJoin(list2, list1);
		if (head != list2 || node3->next != list1 || list1->prev != node3)
			fail("Incorrect behavior when list1 has more than 1 node and "
					"list2 has 1 node.\n");
		hitListDelete(head);
	}

	/* Test behavior when both list1 and list2 have more than 1 node. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 0, 0);
		HitList *node3 = hitListCreateNode(0, 0, 0);
		HitList *node4 = hitListCreateNode(0, 0, 0);
		HitList *node5 = hitListCreateNode(0, 0, 0);
		HitList *node6 = hitListCreateNode(0, 0, 0);
		node1->next = node2;
		node2->next = node3;
		node3->prev = node2;
		node2->prev = node1;
		HitList *list1 = node1;

		node4->next = node5;
		node5->next = node6;
		node6->prev = node5;
		node5->prev = node4;
		HitList *list2 = node4;

		HitList *head = hitListJoin(list1, list2);
		if (head != list1 || node3->next != list2 || list2->prev != node3)
			fail("Incorrect behavior when both list1 and list2 have more  "
					"than 1 node.\n");
		hitListDelete(head);
	}
}
END_TEST


/**
 * Tests hitListGetSize function.
 */
START_TEST(getSize)
{
	/* Test behavior when list is NULL. */
	{
		HitList *list = NULL;
		int size = hitListGetSize(list);
		if (size != 0)
			fail("Incorrect behavior when list is NULL.\n");
	}

	/* Test behavior when list has 1 node. */
	{
		HitList *list = hitListCreateNode(0, 0, 0);
		int size = hitListGetSize(list);
		if (size != 1)
			fail("Incorrect behavior when list has 1 node.\n");
		hitListDelete(list);
	}

	/* Test behavior when list has more than 1 node. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 0, 0);
		HitList *node3 = hitListCreateNode(0, 0, 0);
		node1->next = node2;
		node2->next = node3;
		node3->prev = node2;
		node2->prev = node1;
		HitList *list = node1;
		int size = hitListGetSize(list);
		if (size != 3)
			fail("Incorrect behavior when list has 3 nodes.\n");
		hitListDelete(list);
	}
}
END_TEST


/**
 * Tests hitListCopyNode function.
 */
START_TEST(copyNode)
{
	/* Test whether source node is copied to destination node correctly. */
	{
		HitList *source = hitListCreateNode(2, 2, 2);
		HitList *destination = hitListCreateNode(0, 0, 0);
		destination = hitListCopyNode(destination, source);
		if (destination->index != 2 || destination->shift != 2
				|| destination->offset != 2
				|| destination->next != NULL
				|| destination->prev != NULL)
			fail("Source node was not copied correctly to destination node.\n");
		hitListDelete(source);
		hitListDelete(destination);
	}
}
END_TEST


/**
 * Tests the hitListDuplicateNode function.
 */
START_TEST(duplicateNode)
{
	/* The given node is NULL. */
	{
		HitList *orig = NULL;
		HitList *duplicate = hitListDuplicateNode(orig);
		if (duplicate != NULL)
			fail("Incorrect behavior when original node is NULL.\n");
	}

	/* The given node is not NULL. */
	{
		HitList *orig = hitListCreateNode(2, 3, 4);
		HitList *duplicate = hitListDuplicateNode(orig);
		if (duplicate == NULL || duplicate->index != orig->index
				|| duplicate->shift != orig->shift
				|| duplicate->offset != orig->offset
				|| duplicate->next != NULL || duplicate->prev != NULL)
			fail("Incorrect behaviour when original node is not NULL.\n");
	}

}
END_TEST


/**
 * Tests the hitListPopFront function.
 */
START_TEST(popFront)
{
//	/* Test behavior when list is empty. */
//	{
//		HitList *list = NULL;
//		HitList *result = hitListPopFront(&list);
//		if (result != NULL)
//			fail("Incorrect behavior when list is empty.\n");
//	}

	/* Test behavior when list has only one node. */
	{
		HitList *list = hitListCreateNode(2, 2, 2);
		HitList *result = hitListPopFront(&list);
		if (result == NULL || result->index != 2 || result->shift != 2
					|| result->offset != 2 || result->next
					!= NULL || result->prev != NULL || list != NULL)
			fail("Incorrect behavior when list has only one node.\n");
		hitListDelete(result);
	}

	/* Test behavior when list has more than one node. */
	{
		HitList *node1 = hitListCreateNode(3, 3, 3);
		HitList *node2 = hitListCreateNode(4, 4, 4);
		HitList *node3 = hitListCreateNode(5, 5, 5);
		node1->next = node2;
		node2->next = node3;
		node3->prev = node2;
		node2->prev = node1;
		HitList *list = node1;
		HitList *result = hitListPopFront(&list);
		if (result == NULL || result->index != 3 || result->shift != 3
					|| result->offset != 3 || result->next != NULL
					|| result->prev != NULL || list->prev != NULL
					|| hitListGetSize(list) != 2)
			fail("Incorrect behavior when list has 3 nodes.\n");
		hitListDelete(result);
		hitListDelete(list);
	}
}
END_TEST


/**
 * Tests hitListPushBack function.
 */
START_TEST(pushBack)
{
	/* Test behavior when the list is NULL. */
	{
		HitList *list = NULL;
		HitList *tail = NULL;
		HitList *newNode = hitListCreateNode(0, 0, 0);
		HitList *result = hitListPushBack(list, tail, newNode);
		if (result != newNode || result->prev != NULL || result->next != NULL)
			fail("Incorrect behavior when list is NULL.\n");
		hitListDelete(result);
	}

	/* Test behavior when the new node is NULL. */
	{
		HitList *list = hitListCreateNode(0, 0, 0);
		HitList *tail = NULL;
		HitList *newNode = NULL;
		HitList *result = hitListPushBack(list, tail, newNode);
		if (result != list || result->prev != NULL || result->next != NULL)
			fail("Incorrect behavior when node that is to be appended is "
					"NULL.\n");
		hitListDelete(result);
	}

	/* Test behavior when tail does not point to end of the list. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 0, 0);
		HitList *node3 = hitListCreateNode(0, 0, 0);
		node1->next = node2;
		node2->next = node3;
		node3->prev = node2;
		node2->prev = node1;
		HitList *list = node1;
		HitList *tail = node2;
		HitList *newNode = hitListCreateNode(0, 0, 0);
		HitList *result = hitListPushBack(list, tail, newNode);
		if (result != node1 || node3->next != newNode
				|| newNode->prev != node3)
			fail("Incorrect behavior when tail does not point to the "
					"end of the list.\n");
		hitListDelete(result);
	}

	/* Test behavior when tail points to end of the list. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 0, 0);
		HitList *node3 = hitListCreateNode(0, 0, 0);
		node1->next = node2;
		node2->next = node3;
		node3->prev = node2;
		node2->prev = node1;
		HitList *list = node1;
		HitList *tail = node3;
		HitList *newNode = hitListCreateNode(0, 0, 0);
		HitList *result = hitListPushBack(list, tail, newNode);
		if (result != node1 || node3->next != newNode || newNode->prev != node3)
			fail("Incorrect behavior when tail points to the "
					"end of the list.\n");
		hitListDelete(result);
	}

	/* Test behavior when tail is NULL. */
	{
		HitList *node1 = hitListCreateNode(0, 0, 0);
		HitList *node2 = hitListCreateNode(0, 0, 0);
		HitList *node3 = hitListCreateNode(0, 0, 0);
		node1->next = node2;
		node2->next = node3;
		node3->prev = node2;
		node2->prev = node1;
		HitList *list = node1;
		HitList *tail = NULL;
		HitList *newNode = hitListCreateNode(0, 0, 0);
		HitList *result = hitListPushBack(list, tail, newNode);
		if (result != node1 || node3->next != newNode || newNode->prev != node3)
			fail("Incorrect behavior when tail is NULL.\n");
		hitListDelete(result);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *hitListSuite(void)
{
	Suite *s = suite_create("hitList");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, createNode);
	tcase_add_test(testCaseCore, swapAdjNodes);
	tcase_add_test(testCaseCore, cmpNodes);
	tcase_add_test(testCaseCore, sort);
	tcase_add_test(testCaseCore, join);
	tcase_add_test(testCaseCore, getSize);
	tcase_add_test(testCaseCore, copyNode);
	tcase_add_test(testCaseCore, duplicateNode);
	tcase_add_test(testCaseCore, popFront);
	tcase_add_test(testCaseCore, pushBack);
	suite_add_tcase (s, testCaseCore);

	return s;
}
