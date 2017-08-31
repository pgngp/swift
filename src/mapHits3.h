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

#ifndef MAPHITS3_H_
#define MAPHITS3_H_


/**
 * Binary search tree.
 */
typedef struct bstNode
{
	char refIdx;
	int shift;
	int refPos;
	int numHits;
	struct bstNode *left;
	struct bstNode *right;
} BST;

BST *createBSTNode_wrap(char refIdx, int shift, int refPos);

static BST *createBSTNode(char refIdx, int shift, int refPos);

int isEqual_wrap(BST *node1, BST *node2);

static int isEqual(BST *node1, BST *node2);

BST *search_wrap(BST *searchTree, BST *node, BST **parent);

static BST *search(BST *tree, BST *node, BST **parent);

void mapHits3Delete();

void deleteTree_wrap(BST *tree);

static void deleteTree(BST *tree);

void mapHits3Create(int seedLen);

void mapHits3Reset();

void findBiggestClusters_wrap(BST **bsTree, BST ***biggestClusters,
		int *numBiggestClusters);

static void findBiggestClusters(BST *tree);

int mapHits3GetBestHits(int numBestHits, char *refIdx, int *shift, int *refPos);

static int compare(const void *a, const void *b);

BST *mapHits3AddHit_wrap(char refIdx, int shift, int refPos,
		int *biggestClusterSize);

void mapHits3AddHit(char refIdx, int shift, int refPos);

#endif /* MAPHITS3_H_ */
