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
#include "testRefPosList.h"
#include "../src/refPosList.h"


/**
 * Tests list creation.
 */
START_TEST(create)
{
	/* Check whether a list is created successfully. */
	RefPosList *list = refPosListCreate(0, 0);
	if (list == NULL)
		fail("RefPosList list not created successfully.");

	/* Check whether the 'next' property of the list
	 * that was created is set to NULL or not. */
	if (list->next != NULL)
		fail("The 'next' member of RefPosList is not pointing to NULL.");
	refPosListDelete(list);
}
END_TEST


/**
 * Tests list deletion.
 */
START_TEST(deleteList)
{
	/* Check whether the program crashes if we try to delete
	 * a NULL list. */
	RefPosList *list = NULL;
	refPosListDelete(list);
}
END_TEST


/**
 * Tests refPosListPushBack function.
 */
START_TEST(pushBack)
{
	/* Test behavior when both head and new node point to NULL. */
	{
		RefPosList *head = NULL;
		RefPosList *tail = NULL;
		RefPosList *newNode = NULL;
		RefPosList *result = refPosListPushBack(head, tail, newNode);
		if (result != NULL)
			fail("Incorrect behavior when both head and new node are NULL.\n");
	}

	/* Test behavior when head points to NULL. */
	{
		RefPosList *head = NULL;
		RefPosList *tail = NULL;
		RefPosList *newNode = refPosListCreate(0, 0);
		RefPosList *result = refPosListPushBack(head, tail, newNode);
		if (result == NULL || result != newNode)
			fail("Incorrect behavior when head is NULL.\n");
		refPosListDelete(result);
	}

	/* Test behavior when new node points to NULL. */
	{
		RefPosList *head = refPosListCreate(0, 0);
		RefPosList *tail = head;
		RefPosList *newNode = NULL;
		RefPosList *result = refPosListPushBack(head, tail, newNode);
		if (result == NULL || result != head)
			fail("Incorrect behavior when new node is NULL.\n");
		refPosListDelete(result);
	}

	/* Test behavior when tail is NULL. */
	{
		RefPosList *node1 = refPosListCreate(0, 0);
		RefPosList *node2 = refPosListCreate(0, 10);
		RefPosList *node3 = refPosListCreate(0, 15);
		node1->next = node2;
		node2->next = node3;
		RefPosList *head = node1;
		RefPosList *tail = NULL;
		RefPosList *newNode = refPosListCreate(1, 10);
		RefPosList *result = refPosListPushBack(head, tail, newNode);
		if (result == NULL || result != node1 || result->next != node2
				|| result->next->next != node3
				|| result->next->next->next != newNode)
			fail("Incorrect behavior when tail is NULL.\n");
		refPosListDelete(result);
	}

	/* Test behavior when tail does not point to the last node of the list. */
	{
		RefPosList *node1 = refPosListCreate(0, 0);
		RefPosList *node2 = refPosListCreate(0, 10);
		RefPosList *node3 = refPosListCreate(0, 15);
		node1->next = node2;
		node2->next = node3;
		RefPosList *head = node1;
		RefPosList *tail = node2;
		RefPosList *newNode = refPosListCreate(1, 10);
		RefPosList *result = refPosListPushBack(head, tail, newNode);
		if (result == NULL || result != node1 || result->next != node2
				|| result->next->next != node3
				|| result->next->next->next != newNode)
			fail("Incorrect behavior when tail does not point to the last "
					"node of the list.\n");
		refPosListDelete(result);
	}

	/* Test behavior when list has only 1 node. */
	{
		RefPosList *head = refPosListCreate(0, 0);
		RefPosList *tail = head;
		RefPosList *newNode = refPosListCreate(1, 10);
		RefPosList *result = refPosListPushBack(head, tail, newNode);
		if (result == NULL || result != head || result->next != newNode)
			fail("Incorrect behavior when list has only 1 node.\n");
		refPosListDelete(result);
	}

	/* Test behavior when list has more than 1 nodes. */
	{
		RefPosList *node1 = refPosListCreate(0, 0);
		RefPosList *node2 = refPosListCreate(0, 10);
		RefPosList *node3 = refPosListCreate(0, 15);
		node1->next = node2;
		node2->next = node3;
		RefPosList *head = node1;
		RefPosList *tail = node3;
		RefPosList *newNode = refPosListCreate(1, 10);
		RefPosList *result = refPosListPushBack(head, tail, newNode);
		if (result == NULL || result != node1 || result->next != node2
				|| result->next->next != node3
				|| result->next->next->next != newNode)
			fail("Incorrect behavior when list has more than 1 node.\n");
		refPosListDelete(result);
	}

	/* Test behavior when new node has a child. */
	{
		RefPosList *node1 = refPosListCreate(0, 0);
		RefPosList *node2 = refPosListCreate(0, 10);
		RefPosList *node3 = refPosListCreate(0, 15);
		node1->next = node2;
		node2->next = node3;
		RefPosList *head = node1;
		RefPosList *tail = node3;
		RefPosList *newNode = refPosListCreate(1, 10);
		RefPosList *newNode2 = refPosListCreate(2, 10);
		newNode->next = newNode2;
		RefPosList *result = refPosListPushBack(head, tail, newNode);
		if (result == NULL || result != node1 || result->next != node2
				|| result->next->next != node3
				|| result->next->next->next != newNode
				|| result->next->next->next->next != newNode2)
			fail("Incorrect behavior when list has more than 1 node.\n");
		refPosListDelete(result);
	}
}
END_TEST


/**
 * Tests the refPosListIsEqual function.
 */
START_TEST(isEqual)
{
	/* Both nodes are NULL. */
	{
		RefPosList *node1 = NULL;
		RefPosList *node2 = NULL;
		int result = refPosListIsEqual(node1, node2);
		if (result != -1)
			fail("Incorrect behavior when both nodes are NULL.\n");
	}

	/* Only node1 is NULL. */
	{
		RefPosList *node1 = NULL;
		RefPosList *node2 = refPosListCreate(0, 0);
		int result = refPosListIsEqual(node1, node2);
		if (result != -1)
			fail("Incorrect behavior when only node1 is NULL.\n");
		refPosListDelete(node2);
	}

	/* Only node2 is NULL. */
	{
		RefPosList *node1 = refPosListCreate(0, 0);
		RefPosList *node2 = NULL;
		int result = refPosListIsEqual(node1, node2);
		if (result != -1)
			fail("Incorrect behavior when only node2 is NULL.\n");
		refPosListDelete(node1);
	}

	/* Node1's sequence position is equal to node2's sequence position. */
	{
		RefPosList *node1 = refPosListCreate(0, 0);
		RefPosList *node2 = refPosListCreate(0, 0);
		int result = refPosListIsEqual(node1, node2);
		if (result != 0)
			fail("Incorrect behavior when node1's sequence position "
					"is equal to node2's sequence position.\n");
		refPosListDelete(node1);
		refPosListDelete(node2);
	}

	/* Node1's sequence position is less than node2's sequence position. */
	{
		RefPosList *node1 = refPosListCreate(0, 0);
		RefPosList *node2 = refPosListCreate(0, 5);
		int result = refPosListIsEqual(node1, node2);
		if (result != -1)
			fail("Incorrect behavior when node1's sequence position "
					"is less than node2's sequence position.\n");
		refPosListDelete(node1);
		refPosListDelete(node2);
	}

	/* Node1's sequence position is greater than node2's sequence position. */
	{
		RefPosList *node1 = refPosListCreate(0, 8);
		RefPosList *node2 = refPosListCreate(0, 0);
		int result = refPosListIsEqual(node1, node2);
		if (result != -1)
			fail("Incorrect behavior when node1's sequence position "
					"is greater than node2's sequence position.\n");
		refPosListDelete(node1);
		refPosListDelete(node2);
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *refPosListSuite(void)
{
	Suite *s = suite_create("refPosList");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, create);
	tcase_add_test(testCaseCore, deleteList);
	tcase_add_test(testCaseCore, pushBack);
	tcase_add_test(testCaseCore, isEqual);
	suite_add_tcase(s, testCaseCore);

	return s;
}
