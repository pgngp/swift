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

#include "hitList.h"
#include <stdio.h>
#include <assert.h>


/**
 * Creates a @a HitList node using the given parameters.
 *
 * @param index Reference index.
 * @param shift The position on the reference sequence where the start.
 * position of the query is supposed to align.
 * @param offset The position on the reference sequence where a tuple of the
 * query sequence aligns.
 * @return A @a HitList node.
 */
HitList *hitListCreateNode(char index, int shift, int offset)
{
	HitList *node = (HitList *) malloc(sizeof(HitList));
	node->index = index;
	node->shift = shift;
	node->offset = offset;
	node->next = NULL;
	node->prev = NULL;

	return node;
}


/**
 * Deletes all nodes in the given list.
 *
 * @param list List to be deleted.
 */
void hitListDelete(HitList *list)
{
	if (list == NULL)
		return;

	HitList *tmp = NULL;
	while (list->next != NULL)
	{
		tmp = list->next;
		free(list);
		list = tmp;
	}
	free(list);
}


/**
 * Joins the given 2 lists and returns the pointer to the head of the
 * combined list.
 *
 * @a list2 is appended to @a list1 and the head of @a list1 is returned.
 * If both lists are NULL, returns NULL. If @a list1 is NULL, returns head
 * of @a list2. If @a list2 is NULL, returns head of @a list1. Otherwise,
 * if both lists are not NULL, appends @a list2 to @a list1 and returns
 * head of @a list1.
 *
 * @param list1 Pointer to the head of the first list.
 * @param list2 Pointer to the head of the second list.
 * @return Pointer to the head of the combined list. If both lists are NULL,
 * returns NULL. If @a list1 is NULL, returns head of @a list2.
 * If @a list2 is NULL, returns head of @a list1. Otherwise, if both lists
 * are not NULL, appends @a list2 to @a list1 and returns head of @a list1.
 */
HitList *hitListJoin(HitList *list1, HitList *list2)
{
	if (list1 == NULL && list2 == NULL)
		return NULL;
	else if (list1 == NULL)
		return list2;
	else if (list2 == NULL)
		return list1;

	HitList *tmp = list1;
	while (tmp->next != NULL)
		tmp = tmp->next;
	tmp->next = list2;
	list2->prev = tmp;
	return list1;
}


/**
 * Returns the number of nodes in the given list.
 *
 * @param list Pointer to the head of the list.
 * @return Number of nodes in the given list.
 */
int hitListGetSize(HitList *list)
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
 * Sorts the given list and returns the pointer to the head of the sorted
 * list.
 *
 * Each node in the list is of type (index, shift, offset). This function
 * first sorts by index and then by shift.
 *
 * @param list Pointer to the head of the list that is to be sorted.
 * @param listSize Number of nodes in @a list.
 * @return Pointer to the head of the sorted list.
 */
HitList *hitListSort(HitList *list, uint listSize)
{
	list = mergeSort(list, listSize);
	return list;
}


/**
 * Prints the nodes in the given list.
 *
 * @param list List of nodes to be printed.
 */
void hitListPrintNodes(HitList *list)
{
	if (list == NULL)
		return;
	printf("(%d, %d, %d);\n", list->index, list->shift, list->offset);

	while (list->next != NULL)
	{
		list = list->next;
		printf("(%d, %d, %d);\n", list->index, list->shift, list->offset);
	}
}


/**
 * Sorts the given list using merge sort and returns a pointer to
 * the head of the sorted list. This function assumes that the @a listSize
 * is equal to the number of nodes in @a list. Otherwise, the result
 * may be incorrect.
 *
 * @param list List to be sorted.
 * @param listSize Number of nodes in the given list. This must be equal
 * to the number of nodes in @a list, otherwise the result may be incorrect.
 * @return Pointer to the head of the sorted list.
 */
static HitList *mergeSort(HitList *list, uint listSize)
{
	if (list == NULL || list->next == NULL)
		return list;

	int middle, i;
	HitList *left, *right;

	middle = listSize / 2;
	i = 0;
	left = list;
	right = list;
	while (i < middle)
	{
		right = right->next;
		++i;
	}
	right->prev->next = NULL;
	right->prev = NULL;

	left = mergeSort(left, middle);
	right = mergeSort(right, listSize - middle);

	return merge(left, middle, right, listSize - middle);
}


/**
 * Merges the given lists and returns a pointer to the head of the merged
 * list.
 *
 * @param list1 One of the lists to be merged.
 * @param list1Size Number of nodes in @a list1.
 * @param list2 Other of the lists to be merged.
 * @param list2Size Number of nodes in @a list2.
 * @return Pointer to the head of the merged list.
 */
static HitList *merge(HitList *list1, uint list1Size, HitList *list2,
		uint list2Size)
{
	static HitList *result, *tmpList, *tail;
	result = NULL;
	tmpList = NULL;
	tail = NULL;

	while (list1Size > 0 || list2Size > 0)
	{
		if (list1Size > 0 && list2Size > 0)
		{
			if (hitListCmpNodes(list1, list2) <= 0)
			{
				tmpList = hitListPopFront(&list1);
				--list1Size;
			}
			else
			{
				tmpList = hitListPopFront(&list2);
				--list2Size;
			}
		}
		else if (list1Size > 0)
		{
			tmpList = hitListPopFront(&list1);
			--list1Size;
		}
		else if (list2Size > 0)
		{
			tmpList = hitListPopFront(&list2);
			--list2Size;
		}
		result = hitListPushBack(result, tail, tmpList);
		tail = tmpList;
	}

	return result;
}


/**
 * Removes and returns the first node of the given list. The size of the given
 * list is reduced by 1.
 *
 * If the pop operation is successful, then @a list will point to the
 * next node in the list if it exists. If the pop operation
 * is succesful but there no other node in @a list, then this pointer
 * will point to NULL.
 *
 * @param[in,out] list List from which the front node is to be removed and
 * returned. If the pop operation is successful, then this pointer will
 * point to the next node in the list if it exists. If the pop operation
 * is succesful but there no other node in the list, then this pointer
 * will point to NULL.
 * the next
 * @return The removed node from the given list.
 */
HitList *hitListPopFront(HitList **list)
{
	if (*list == NULL)
		return NULL;

	static HitList *tmp;
	tmp = *list;
	if ((*list)->next != NULL)
	{
		(*list) = (*list)->next;
		(*list)->prev = NULL;
	}
	else
		(*list) = NULL;

	tmp->next = NULL;
	tmp->prev = NULL;

	return tmp;
}


/**
 * Copies data from the node pointed to by @a source to the node pointed
 * to by @a destination and returns a pointer to the @a destination.
 *
 * It is assumed that both @a source and @a destination will be valid
 * pointers. If @a source contains more than 1 node, only the node pointed
 * by @a source will be copied.
 *
 * @param[out] destination The destination node where the contents of source
 * node would be copied to.
 * @param source The source node from where the contents would be copied to
 * the destination node.
 * @return Pointer to the destination node.
 */
HitList *hitListCopyNode(HitList *destination, const HitList *source)
{
	assert(destination != NULL);
	assert(source != NULL);

	destination->index = source->index;
	destination->shift = source->shift;
	destination->offset = source->offset;
	destination->next = NULL;
	destination->prev = NULL;

	return destination;
}


/**
 * Creates a duplicate of the given node and returns a pointer to the
 * duplicate node.
 *
 * @note The @a prev and @a next data members will not be copied from the
 * original node to the duplicate node.
 *
 * @param orig The node that is to be duplicated.
 * @return The duplicate node.
 */
HitList *hitListDuplicateNode(const HitList *orig)
{
	if (orig == NULL)
		return NULL;

	HitList *node = (HitList *) malloc(sizeof(HitList));
	node->index = orig->index;
	node->shift = orig->shift;
	node->offset = orig->offset;
	node->next = NULL;
	node->prev = NULL;

	return node;
}


/**
 * Swaps adjacent nodes and returns the new left node.
 *
 * @note If the nodes are not adjacent or if @a leftNode and @a rightNode
 * are not in that order in the list, the result might be incorrect.
 *
 * @param leftNode The left node of the two adjacent nodes.
 * @param rightNode The right node of the two adjacent nodes.
 * @return The new left node.
 */
HitList *hitListSwapAdjNodes(HitList *leftNode, HitList *rightNode)
{
	if (leftNode == NULL || rightNode == NULL)
		return NULL;

	static HitList *parent, *child;

	parent = leftNode->prev;
	child = rightNode->next;
	if (parent != NULL)
	{
		parent->next = rightNode;
		rightNode->prev = parent;
	}
	else
		rightNode->prev = NULL;

	if (child != NULL)
	{
		leftNode->next = child;
		child->prev = leftNode;
	}
	else
		leftNode->next = NULL;
	rightNode->next = leftNode;
	leftNode->prev = rightNode;

	return rightNode;
}


/**
 * Compares the given two nodes and returns 0 if they are "equal",
 * -1 if @a node1 is "smaller" than @a node2, and 1 if @a node1 is
 * "greater" than @a node2.
 *
 * In order to compare the two nodes, the index values are looked at first,
 * then the shift values, and finally the offset values of each node.
 *
 * @param node1 One of the two nodes to be compared.
 * @param node2 Other of the two nodes to be compared.
 * @return 0 if the two nodes are "equal", -1 if @a node1 is smaller than
 * @a node2, and 1 if @a node1 is greater than @a node2.
 */
int hitListCmpNodes(HitList *node1, HitList *node2)
{
	if (node1 == NULL || node2 == NULL)
		return -2;

	/* Node1 index greater than node2 index. */
	if (node1->index > node2->index)
		return 1;
	/* Node1 index less than node2 index. */
	else if (node1->index < node2->index)
		return -1;
	/* Node1 index equal to node2 index. */
	else
	{
		/* Node1 shift greater than node2 shift. */
		if (node1->shift > node2->shift)
			return 1;
		/* Node1 shift less than node2 shift. */
		else if (node1->shift < node2->shift)
			return -1;
		/* Node1 shift equal to node2 shift. */
		else
		{
			/* Node1 offset greater than node2 offset. */
			if (node1->offset > node2->offset)
				return 1;
			/* Node1 offset less than node2 offset. */
			else if (node1->offset < node2->offset)
				return -1;
			/* Node1 offset equal to node2 offset. */
			else
				return 0;
		}
	}
}


/**
 * Appends @a node to the end of @a list and returns the head of the combined
 * list.
 *
 * If @a tail is NULL or if it does not point to the end of @ list,
 * then the @a list will be parsed until the end of the @a list is found. It
 * is assumed that @a head points to the head of the list.
 *
 * @param head Head of the list at the end of which the @a node
 * is to be appended. It is assumed that @a head points to the beginning
 * of the list. @a head can be NULL in which case @a node will be returned.
 * @param tail Pointer to the end of @a list. If @a tail is NULL or if it
 * does not point to the end of @a list, then @a list will be parsed
 * until the end of @a list is found.
 * @param node The node to be appended to the end of @a list.
 * @return Head of the combined list. If @a head is NULL, @a node will
 * be returned.
 */
HitList *hitListPushBack(HitList *head, HitList *tail, HitList *node)
{
	if (head == NULL && node == NULL)
		return NULL;
	else if (head == NULL)
		return node;
	else if (node == NULL)
		return head;

	/* If tail is not NULL. */
	if (tail != NULL && tail->next == NULL)
	{
		tail->next = node;
		node->prev = tail;
	}
	/* If tail is NULL. */
	else
	{
		tail = head;
		while (tail->next != NULL)
			tail = tail->next;
		tail->next = node;
		node->prev = tail;
	}

	return head;
}

