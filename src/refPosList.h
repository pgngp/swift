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

#ifndef REFPOSLIST_H_
#define REFPOSLIST_H_


/**
 * @struct refPosNode
 *
 * Stores the location of the reference sequence.
 */
typedef struct refPosNode
{
	char index; /**< Reference sequence index. */
	int position; /**< Position in the reference sequence. */
	struct refPosNode *next; /**< Points to the next @em refPosNode node. */
} RefPosList;


RefPosList *refPosListCreate(char index, int position);
void refPosListDelete(RefPosList *list);
void refPosListPrint(RefPosList *list);
RefPosList *refPosListJoin(RefPosList *list1, RefPosList *list2);
int refPosListGetSize(RefPosList *list);
RefPosList *refPosListPushBack(RefPosList *head, RefPosList *tail,
		RefPosList *newNode);
int refPosListIsEqual(const RefPosList *node1, const RefPosList *node2);

#endif /* REFPOSLIST_H_ */
