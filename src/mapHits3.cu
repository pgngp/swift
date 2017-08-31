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

#include "mapHits3.h"
#include "common.h"
#include "array.h"
#include <stdio.h>


static BST *tree = NULL; /* Points to the head of the tree. */
static int _biggestClusterSize = 1; /* Biggest cluster size. */
static BST **_biggestClusters = NULL; /* Holds pointers to biggest clusters. */
static int _numBiggestClusters = 0; /* Number of biggest clusters. */
static int *_randNumArr = NULL; /* Random number array. */


/**
 * Initializes the data structures in this file.
 */
void mapHits3Create(int seedLen)
{
	tree = NULL;
	_biggestClusterSize = 1;
	int numQryTuples = MAX_QRY_SEQ_LENGTH - seedLen;
	int numElements = numQryTuples * TUPLE_IGNORE_THRES;
	_biggestClusters = (BST **) calloc(numElements, sizeof(BST *));
	_numBiggestClusters = 0;
	_randNumArr = (int *) calloc(numElements, sizeof(int));
}


/**
 * Resets key variables in this file.
 */
void mapHits3Reset()
{
	deleteTree(tree);
	tree = NULL;
	_biggestClusterSize = 1;
	_numBiggestClusters = 0;
}


/**
 * Deletes resources acquired in this file.
 */
void mapHits3Delete()
{
	deleteTree(tree);
	tree = NULL;
	_biggestClusterSize = 1;
	free(_biggestClusters);
	_biggestClusters = NULL;
	_numBiggestClusters = 0;
	free(_randNumArr);
	_randNumArr = NULL;
}


/**
 * This is a wrapper function that wraps @a deleteTree function. It has been
 * added so that functions in this file can be unit-tested.
 *
 * @param	tree	The tree to be deleted.
 */
void deleteTree_wrap(BST *tree)
{
	deleteTree(tree);
}


/**
 * Deletes the binary search tree pointed to by @a tree. This function uses
 * post-order tree walk.
 *
 * @param	tree	The tree to be deleted.
 */
static void deleteTree(BST *tree)
{
	if (tree != NULL)
	{
		deleteTree(tree->left);
		deleteTree(tree->right);
		free(tree);
	}
}


/**
 * This is a wrapper function that wraps @a mapHits3AddHit_wrap function. It
 * has been added so that @a mapHits3AddHit_wrap can be unit-tested.
 *
 * @param	refIdx	Reference index.
 * @param	shift	Shift (reference position - query position).
 * @param	refPos	Reference position.
 * @param[out]	biggestClusterSize	Biggest cluster size.
 * @return	The root of the tree.
 */
BST *mapHits3AddHit_wrap(char refIdx, int shift, int refPos,
		int *biggestClusterSize)
{
	mapHits3AddHit(refIdx, shift, refPos);
	*biggestClusterSize = _biggestClusterSize;
	return tree;
}


/**
 * Add the given hit to the tree.
 *
 * @param	refIdx	Reference index.
 * @param	shift	Shift (reference position - query position).
 * @param	refPos	Reference position.
 */
void mapHits3AddHit(char refIdx, int shift, int refPos)
{
	static BST *newNode;
	newNode = createBSTNode(refIdx, shift, refPos);
	if (tree == NULL)
	{
		tree = newNode;
		return;
	}

	static BST *parent, *node;
	node = search(tree, newNode, &parent);
	/* The "hit" is already present in the tree. */
	if (node != NULL)
	{
		free(newNode);
		if (node->refPos > refPos)
			node->refPos = refPos;
		++node->numHits;
		if (_biggestClusterSize < node->numHits)
			_biggestClusterSize = node->numHits;
	}
	/* The "hit" is not already present in the tree. */
	else
	{
		if (isEqual(newNode, parent) < 0)
			parent->left = newNode;
		else
			parent->right = newNode;
	}
}


/**
 * This is a wrapper function that wraps the @a search function. It has been
 * added so that @a search function can be unit-tested.
 *
 * @param 		tree	The tree in which @a node will be searched.
 * @param 		node	The node to be searched.
 * @param[out]	parent	Parent node under which a new node could be added.
 * @return		Pointer to the node that was found in the tree, otherwise a
 * pointer to the parent under which the given node can be added.
 */
BST *search_wrap(BST *searchTree, BST *node, BST **parent)
{
	return search(searchTree, node, parent);
}


/**
 * Searches for the given @a node in the given @a searchTree and returns the
 * pointer to the node if it's found in the tree, otherwise returns a pointer to
 * the parent node under which this node could be added.
 *
 * @param 		tree	The tree in which @a node will be searched.
 * @param 		node	The node to be searched.
 * @param[out]	parent	Parent node under which a new node could be added.
 * @return		Pointer to the node that was found in the tree, otherwise a
 * pointer to the parent under which the given node can be added.
 */
static BST *search(BST *searchTree, BST *node, BST **parent)
{
	static int cmpValue;
	while (searchTree != NULL)
	{
		cmpValue = isEqual(node, searchTree);
		*parent = searchTree;
		if (cmpValue < -1)
			break;
		else if (cmpValue < 0)
			searchTree = searchTree->left;
		else if (cmpValue > 0)
			searchTree = searchTree->right;
		else
			return searchTree;
	}
	return NULL;
}


/**
 * This is a wrapper function that wraps @a isEqual function. It has been
 * added so that @a isEqual function can be unit-tested.
 *
 * @param	node1	One of the nodes to be compared.
 * @param	node2	Another of the nodes to be compared.
 * @return	0 if both nodes are equal, -1 if @a node1 is less than @a node2,
 * and 1 if @a node1 is greater than @a node2.
 */
int isEqual_wrap(BST *node1, BST *node2)
{
	return isEqual(node1, node2);
}


/**
 * Compares @a node1 with @a node2 and returns 0 if they are equal, -1 if
 * @a node1 is less than @a node2, and 1 if @a node1 is greater than @a node2.
 *
 * @note This function does not use "reference position" for comparison.
 *
 * @param	node1	One of the nodes to be compared.
 * @param	node2	Another of the nodes to be compared.
 * @return	0 if both nodes are equal, -1 if @a node1 is less than @a node2,
 * and 1 if @a node1 is greater than @a node2.
 */
static int isEqual(BST *node1, BST *node2)
{
	if (node1 == NULL || node2 == NULL)
		return -2;
	else if (node1->refIdx < node2->refIdx || (node1->refIdx == node2->refIdx
			&& node1->shift < node2->shift))
		return -1;
	else if (node1->refIdx > node2->refIdx || (node1->refIdx == node2->refIdx
			&& node1->shift > node2->shift))
		return 1;
	else
		return 0;
}


/**
 * Returns the best-mapped hits.
 *
 * @param 		numBestHits	Max number of best hits to be returned.
 * @param[out]	refIdx		Reference index.
 * @param[out]	shift		Shift (reference position - query position).
 * @param[out]	refPos		Reference position.
 * @return		Number of best hits.
 */
int mapHits3GetBestHits(int numBestHits, char *refIdx, int *shift, int *refPos)
{
	/* Walk the tree and find the largest clusters. */
	findBiggestClusters(tree);

	/* Randomly select the desired number of best matches. */
	if (numBestHits > _numBiggestClusters)
		numBestHits = _numBiggestClusters;
	arrGetRandomNums((uint) numBestHits, 0, (uint) (_numBiggestClusters - 1),
			(uint *) _randNumArr);
	qsort(_randNumArr, numBestHits, sizeof(int), compare);
	static int i;
	static BST *tmpNode;
	for (i = 0; i < numBestHits; ++i)
	{
		tmpNode = _biggestClusters[_randNumArr[i]];
		refIdx[i] = tmpNode->refIdx;
		shift[i] = tmpNode->shift;
		refPos[i] = tmpNode->refPos;
	}
	return numBestHits;
}


/**
 * This is a wrapper function that wraps @a findBiggestClusters function. It has
 * been added so that @a findBiggestClusters can be unit-tested.
 *
 * @param[out]	bsTree			The tree in which biggest clusters are to be
 * looked for.
 * @param[out]	biggestClusters	Array containing pointers to the biggest
 * clusters.
 * @param[out]	numBiggestClusters	Number of biggest clusters.
 */
void findBiggestClusters_wrap(BST **bstTree, BST ***biggestClusters,
		int *numBiggestClusters)
{
	findBiggestClusters(tree);
	*bstTree = tree;
	*biggestClusters = _biggestClusters;
	*numBiggestClusters = _numBiggestClusters;
}


/**
 * Finds the biggest clusters. This function uses in-order traversal to parse
 * the tree.
 *
 * @param bsTree	The tree in which biggest clusters are to be looked for.
 */
static void findBiggestClusters(BST *bsTree)
{
	if (bsTree != NULL)
	{
		if (_biggestClusterSize == bsTree->numHits)
		{
			_biggestClusters[_numBiggestClusters] = bsTree;
			++_numBiggestClusters;
		}
		findBiggestClusters(bsTree->left);
		findBiggestClusters(bsTree->right);
	}
}


/**
 * This is a wrapper function that wraps @a createBSTNode. It has been
 * added so that functions in this file can be unit-tested.
 *
 * @param 	refIdx	Reference index.
 * @param 	shift	Shift (reference position - query position).
 * @param	refPos	Reference position.
 * @return	A new BST node.
 */
BST *createBSTNode_wrap(char refIdx, int shift, int refPos)
{
	return createBSTNode(refIdx, shift, refPos);
}


/**
 * Creates a BST node using the given parameters.
 *
 * @param 	refIdx	Reference index.
 * @param 	shift		Shift (reference position - query position).
 * @param	refPos	Reference position.
 * @return	A new BST node.
 */
static BST *createBSTNode(char refIdx, int shift, int refPos)
{
	BST *node = (BST *) malloc(sizeof(BST));
	node->refIdx = refIdx;
	node->shift = shift;
	node->refPos = refPos;
	node->numHits = 1;
	node->left = NULL;
	node->right = NULL;
	return node;
}


/**
 * Returns the integer difference between the given two parameters.
 *
 * @param a First parameter.
 * @param b Second parameter.
 * @return The difference between the given two parameters.
 */
static int compare(const void *a, const void *b)
{
	return (*(int *)a - *(int *)b);
}


