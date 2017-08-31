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

#include "refPosList.h"
#include <stdio.h>


/**
 * Creates a @a refPosList list.
 *
 * @param index Index of the reference sequence
 * @param position Position of the sequence fragment with respect to the
 * beginning of the reference sequence.
 * @return A new @a RefPosList list.
 */
RefPosList *refPosListCreate(char index, int position)
{
	RefPosList *list = (RefPosList *) malloc(sizeof(RefPosList));
	list->index = index;
	list->position = position;
	list->next = NULL;

	return list;
}


/**
 * Deletes the given list.
 *
 * @param list List to be deleted.
 */
void refPosListDelete(RefPosList *list)
{
	if (list == NULL)
		return;

	RefPosList *tmp = NULL;
	while (list->next != NULL)
	{
		tmp = list->next;
		free(list);
		list = tmp;
	}
	free(list);
}


/**
 * Prints the nodes in the given list
 *
 * @param list List of nodes to be printed
 */
void refPosListPrint(RefPosList *list)
{
	if (list == NULL)
		return;
	printf("(%d, %d);", list->index, list->position);

	while (list->next != NULL)
	{
		list = list->next;
		printf("(%d, %d);", list->index, list->position);
	}
	printf("\n");
}


/**
 * Joins the given 2 lists and returns the pointer to the head of the
 * combined list
 *
 * @param list1 Pointer to the head of the first list
 * @param list2 Pointer to the head of the second list
 * @return Pointer to the head of the combined list
 */
RefPosList *refPosListJoin(RefPosList *list1, RefPosList *list2)
{
	if (list1 == NULL && list2 == NULL)
		return NULL;
	else if (list1 == NULL)
		return list2;
	else if (list2 == NULL)
		return list1;

	RefPosList *tmp = list1;
	while (tmp->next != NULL)
		tmp = tmp->next;
	tmp->next = list2;

	return list1;
}


/**
 * Returns the size of the list
 *
 * @param list Pointer to the head of the list
 * @return Size of the list
 */
int refPosListGetSize(RefPosList *list)
{
	if (list == NULL)
		return 0;

	int count = 1;
	while (list->next != NULL)
	{
		++count;
		list = list->next;
	}

	return count;
}


/**
 * Appends @a newNode to the end of the list pointed to by @a head and
 * returns the head of the combined list. It is assumed that @a head
 * points to the beginning of the list.
 *
 * @param head Head of the list to which @a newNode is to be appended.
 * If this is NULL, @a newNode will be returned.
 * @param tail Pointer to the last node of the list. If this is NULL,
 * then the whole list pointed to by @a head will be parsed until the
 * last node is found, which will slow down the program.
 * @param newNode Node to be appended to the list pointed to by @a head.
 * @return A pointer to the head of the combined list.
 */
RefPosList *refPosListPushBack(RefPosList *head, RefPosList *tail,
		RefPosList *newNode)
{
	if (head == NULL && newNode == NULL)
		return NULL;
	else if (head == NULL)
		return newNode;
	else if (newNode == NULL)
		return head;

	if (tail != NULL && tail->next == NULL)
		tail->next = newNode;
	else if (tail != NULL && tail->next != NULL)
	{
		while (tail->next != NULL)
			tail = tail->next;
		tail->next = newNode;
	}
	else
	{
		tail = head;
		while (tail->next != NULL)
			tail = tail->next;
		tail->next = newNode;
	}

	return head;
}


/**
 * Compares the sequence positions of @a node1 and @a node2 and
 * returns 0 if the sequence position of @a node1 is equal to the
 * sequence position of @a node2, otherwise returns -1.
 *
 * @param node1 One of the two nodes.
 * @param node2 Other of the two nodes.
 * @return 0 if the sequence position of @a node1 is equal to the
 * sequence position of @a node2, otherwise returns -1.
 */
int refPosListIsEqual(const RefPosList *node1, const RefPosList *node2)
{
	if (node1 == NULL || node2 == NULL)
		return -1;

	if (node1->index == node2->index && node1->position == node2->position)
		return 0;
	else
		return -1;
}


