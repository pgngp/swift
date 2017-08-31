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

#include "common.h"
#include "mapHits5.h"
#include "array.h"
#include <stdio.h>

//#define	CLUSTER_DISTANCE	50
//#define	CLUSTER_DISTANCE	75
#define CLUSTER_DISTANCE        25


static char *_refIdx = NULL; /* Reference index */
static int *_refPos = NULL; /* Reference position */
static int _numElements = 0; /* Number of elements in the above arrays. */
static int _size = 0; /* Size of the current mapped list. */
static char _tmp = '\0'; /* Temporary variable. */
static int _tmp2 = 0; /* Temporary variable. */
static int *_cluster = NULL; /* Stores the cluster indexes. */
static int *_clusterSize = NULL; /* Stores the cluster sizes. */
static int _numClusters = 0; /* Number of clusters. */
static int _biggestClusterSize = 0; /* Biggest cluster size. */
static int *_biggestCluster = NULL; /* Stores index of the biggest clusters. */
static int _numBiggestClusters = 0; /* Number of biggest clusters. */
static int *_randNumArr = NULL; /* Random number array. */


/**
 * Creates a list of mapped hits.
 *
 * @param seedLen	Seed length.
 */
void mapHits5Create(int seedLen)
{
	int numQryTuples = MAX_QRY_SEQ_LENGTH - seedLen;
	_numElements = numQryTuples * TUPLE_IGNORE_THRES;
	_refIdx = (char *) calloc(_numElements, sizeof(char));
	_refPos = (int *) calloc(_numElements, sizeof(int));
	_cluster = (int *) calloc(_numElements, sizeof(int));
	_clusterSize = (int *) calloc(_numElements, sizeof(int));
	_biggestCluster = (int *) calloc(_numElements, sizeof(int));
	_randNumArr = (int *) calloc(_numElements, sizeof(int));
}


/**
 * Reset key variables in this file.
 */
void mapHits5Reset()
{
	_size = 0;
	_numClusters = 0;
	_biggestClusterSize = 0;
	_numBiggestClusters = 0;
}


/**
 * Deletes mapped hits.
 */
void mapHits5Delete()
{
	free(_refIdx);
	_refIdx = NULL;
	free(_refPos);
	_refPos = NULL;
	free(_cluster);
	_cluster = NULL;
	free(_clusterSize);
	_clusterSize = NULL;
	free(_biggestCluster);
	_biggestCluster = NULL;
	free(_randNumArr);
	_randNumArr = NULL;
	_numElements = 0;
	_size = 0;
	_tmp = '\0';
	_tmp2 = 0;
	_numClusters = 0;
	_biggestClusterSize = 0;
	_numBiggestClusters = 0;
}


/**
 * This is a wrapper function that wraps @a mapHits5AddHit function. It has
 * been added so that @a mapHits5AddHit can be unit-tested.
 *
 * @param[out]	refIdxArr	Reference index array.
 * @param[out]	refPosArr	Reference position array.
 * @param[out]	size		Number of elements in the @a refIdxArr.
 * @param 		refIdx		Reference index to be added.
 * @param 		refPos		Reference position to be added.
 */
void mapHits5AddHit_wrap(char **refIdxArr, int **refPosArr, int *size,
		char refIdx, int refPos)
{
	mapHits5AddHit(refIdx, refPos);
	*refIdxArr = _refIdx;
	*refPosArr = _refPos;
	*size = _size;
}


/**
 * Adds the given hit to the list.
 *
 * @param refIdx	Reference index.
 * @param refPos	Reference position.
 */
void mapHits5AddHit(char refIdx, int refPos)
{
	_refIdx[_size] = refIdx;
	_refPos[_size] = refPos;
	++_size;
}


/**
 * Returns the best mapped hits.
 *
 * @param 		numBestHits			Number of best hits to be returned.
 * @param[out]	refIdx				Best-matching reference index. The number of
 * elements in this array must be equal to @a numBestHits.
 * @param[out]	refPos				Best-matching reference position. The
 * number of elements in this array must be equal to @a numBestHits.
 * @return		Number of best-hits.
 */
int mapHits5GetBestHits(int numBestHits, char *refIdx, int *refPos)
{
	/* Return 0 if there are no hits. */
	if (_size == 0)
		return 0;

	/* Sort the list of mapped hits. */
	quickSort(0, _size - 1);

	/* Create clusters. */
	createClusters();

	/* Find the biggest cluster size. */
	findBiggestClusterSize();

	/* Find the first hit of the biggest clusters. */
	return findBiggestClusters(numBestHits, refIdx, refPos);
}


/**
 * This is a wrapper function that wraps @a findBiggestClusters function. It
 * has been added so that @a findBiggestClusters can be unit-tested.
 *
 * @param 		numBestHits			Number of biggest clusters that are to
 * searched for.
 * @param[out]	refIdx				Best-matching reference index. The number
 * of elements in this array must be equal to @a numBestHits.
 * @param[out]	refPos				Best-matching reference position. The number
 * of elements in this array must be equal to @a numBestHits.
 * @return 		Number of best-hits in @a refIdx_bestMatch.
 */
int findBiggestClusters_wrap5(int numBestHits, char *refIdx, int *refPos)
{
	createClusters();
	findBiggestClusterSize();
	return findBiggestClusters(numBestHits, refIdx, refPos);
}


/**
 * Finds the biggest clusters.
 *
 * @param 		numBestHits	Number of biggest clusters that are to
 * searched for.
 * @param[out]	refIdx		Best-matching reference index. The number
 * of elements in this array must be equal to @a numBestHits.
 * @param[out]	refPos		Best-matching reference position. The number
 * of elements in this array must be equal to @a numBestHits.
 * @return 		Number of best-hits in @a refIdx_bestMatch.
 */
static int findBiggestClusters(int numBestHits, char *refIdx, int *refPos)
{
	static int i, j;
	j = 0;
	for (i = 0; i < _numClusters; ++i)
	{
		if (_clusterSize[i] == _biggestClusterSize)
		{
			_biggestCluster[j] = _cluster[i];
			++j;
		}
	}
	_numBiggestClusters = j;

	if (numBestHits > _numBiggestClusters)
		numBestHits = _numBiggestClusters;

	arrGetRandomNums((uint) numBestHits, 0, (uint) (_numBiggestClusters - 1),
			(uint *) _randNumArr);
	qsort(_randNumArr, numBestHits, sizeof(int), compare);
	static int tmp;
	for (i = 0; i < numBestHits; ++i)
	{
		tmp = _biggestCluster[_randNumArr[i]];
		refIdx[i] = _refIdx[tmp];
		refPos[i] = _refPos[tmp];
	}
	return numBestHits;
}


/**
 * This is a wrapper function that wraps @a findBiggestClusterSize function.
 * It has been added so that @a findBiggestClusterSize can be unit-tested.
 */
int findBiggestClusterSize_wrap5()
{
	createClusters();
	findBiggestClusterSize();
	return _biggestClusterSize;
}


/**
 * Find the biggest cluster size.
 */
static void findBiggestClusterSize()
{
	if (_numClusters == 1)
	{
		_biggestClusterSize = _size;
		_clusterSize[0] = _size;
	}
	else
	{
		int i;
		for (i = 1; i < _numClusters; ++i)
		{
			_clusterSize[i - 1] = _cluster[i] - _cluster[i - 1];
			if (_biggestClusterSize < _clusterSize[i - 1])
				_biggestClusterSize = _clusterSize[i - 1];
		}
		_clusterSize[i - 1] = _size - _cluster[i - 1];
		if (_biggestClusterSize < _clusterSize[i - 1])
			_biggestClusterSize = _clusterSize[i - 1];
	}
}


/**
 * This is a wrapper function that wraps @a createClusters function. It has
 * been added so that @a createClusters can be unit-tested.
 */
void createClusters_wrap5(int **cluster)
{
	createClusters();
	*cluster = _cluster;
}


/**
 * Creates clusters.
 */
static void createClusters()
{
	_cluster[0] = 0;
	++_numClusters;
	static int i, j;
	j = 1;
	for (i = 1; i < _size; ++i)
	{
//		if (_refIdx[i] == _refIdx[i - 1] && _shift[i] == _shift[i - 1])
		if (_refIdx[i] == _refIdx[i - 1]
		    && (_refPos[i] - _refPos[i - 1]) <= CLUSTER_DISTANCE)
		{ }
		else
		{
			_cluster[j] = i;
			++j;
		}
	}
	_numClusters = j;
}


/**
 * This is a wrapper function that wraps @a quickSort function. It has
 * been added so that @a quickSort function can be unit-tested.
 *
 * @param[out]	refIdx	Reference index.
 * @param[out]	refPos	Reference position.
 */
void quickSort_wrap5(char **refIdx, int **refPos)
{
	quickSort(0, _size - 1);
	*refIdx = _refIdx;
	*refPos = _refPos;
}


/**
 * Sorts the list of mapped hits using the QuickSort algorithm.
 *
 * @param left		Beginning index
 * @param right		Ending index
 */
static void quickSort(int left, int right)
{
	int i, last;

	if (left >= right)
		return;
	swap(left, (left + right) / 2);
	last = left;
	for (i = left + 1; i <= right; ++i)
	{
		if ((_refIdx[i] < _refIdx[left])
				|| (_refIdx[i] == _refIdx[left] && _refPos[i] < _refPos[left]))
			swap(++last, i);
	}
	swap(left, last);
	quickSort(left, last - 1);
	quickSort(last + 1, right);
}


/**
 * This is a wrapper function that wraps @a swap function. It has been added
 * so that @a swap function can be unit-tested.
 *
 * @param 		i		Index of one of the elements to be swapped.
 * @param 		j		Index of the other element to be swapped.
 * @param[out]	refIdx	Reference index.
 * @param[out]	refPos	Reference position.
 */
void swap_wrap5(int i, int j, char **refIdx, int **refPos)
{
	swap(i, j);
	*refIdx = _refIdx;
	*refPos = _refPos;
}


/**
 * Swaps element with index @a i with element with index @a j.
 *
 * @param i	Index of the first element to be swapped.
 * @param j	Index of the second element to be swapped.
 */
static void swap(int i, int j)
{
	_tmp = _refIdx[i];
	_refIdx[i] = _refIdx[j];
	_refIdx[j] = _tmp;

	_tmp2 = _refPos[i];
	_refPos[i] = _refPos[j];
	_refPos[j] = _tmp2;
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


/**
 * Prints the mapped hits.
 */
void mapHits5Print()
{
	fprintf(stderr, "************************************\n");
	fprintf(stderr, "refIdx\tshift\trefPos\n");
	int i;
	for (i = 0; i < _size; ++i)
		fprintf(stderr, "%d\t%d\n", _refIdx[i], _refPos[i]);
	fprintf(stderr, "************************************\n");
}
