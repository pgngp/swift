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

#ifndef HISTLIST_H_
#define HISTLIST_H_

#include "refPosList.h"
#include "common.h"

typedef struct hitNode
{
	char index; /**< Reference sequence ID. */
	int shift; /**< Reference tuple position - query tuple position. */
	int offset; /**< Reference tuple position. */
	struct hitNode *next; /**< Pointer to the next node. */
	struct hitNode *prev; /**< Pointer to the previous node. */
} HitList;

HitList *hitListCreateNode(char index, int shift, int offset);
void hitListDelete(HitList *list);
HitList *hitListJoin(HitList *list1, HitList *list2);
int hitListGetSize(HitList *list);
HitList *hitListSort(HitList *list, uint listSize);
void hitListPrintNodes(HitList *list);
static HitList *mergeSort(HitList *list, uint listSize);
static HitList *merge(HitList *list1, uint list1Size,
		HitList *list2, uint list2Size);
HitList *hitListCopyNode(HitList *destination, const HitList *source);
HitList *hitListDuplicateNode(const HitList *orig);
HitList *hitListSwapAdjNodes(HitList *leftNode, HitList *rightNode);
int hitListCmpNodes(HitList *node1, HitList *node2);
HitList *hitListPopFront(HitList **list);
HitList *hitListPushBack(HitList *head, HitList *tail, HitList *node);

#endif /* HISTLIST_H_ */
