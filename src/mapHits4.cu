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

#include "mapHits4.h"
#include "common.h"
#include "array.h"
#include <stdio.h>


static ulong *_hits = NULL; /* Composite reference position containing
reference index and position. */
static int _maxNumElements = 0; /* Max number of elements in @a _hits. */
static int _size = 0; /* Number of valid elements in @a _hits. */
static int *_clusters = NULL; /* Stores the cluster indexes. */
static int _numClusters = 0; /* Number of clusters. */
static int *_clusterSize = NULL; /* Stores the cluster sizes. */
static int _biggestClusterSize = 0; /* Biggest cluster size. */
static int *_biggestCluster = NULL; /* Stores index of the biggest clusters. */
static int _numBiggestClusters = 0; /* Number of biggest clusters. */
static int *_randNumArr = NULL; /* Random number array. */


/**
 * Initializes data structures in this file.
 *
 * @param seedLen	Seed length.
 */
void mapHits4Create(int seedLen)
{
	int numQryTuples = MAX_QRY_SEQ_LENGTH - seedLen;
	_maxNumElements = numQryTuples * TUPLE_IGNORE_THRES;
	_hits = (ulong *) calloc(_maxNumElements, sizeof(ulong));
	_clusters = (int *) calloc(_maxNumElements, sizeof(int));
	_clusterSize = (int *) calloc(_maxNumElements, sizeof(int));
	_biggestCluster = (int *) calloc(_maxNumElements, sizeof(int));
	_randNumArr = (int *) calloc(_maxNumElements, sizeof(int));
}


/**
 * Deletes data structures in this file.
 */
void mapHits4Delete()
{
	_maxNumElements = 0;
	_size = 0;
	_numClusters = 0;
	_biggestClusterSize = 0;
	_numBiggestClusters = 0;

	free(_hits);
	_hits = NULL;
	free(_clusters);
	_clusters = NULL;
	free(_clusterSize);
	_clusterSize = NULL;
	free(_biggestCluster);
	_biggestCluster = NULL;
	free(_randNumArr);
	_randNumArr = NULL;
}


/**
 * Resets the variables data structures in this file to their default values.
 */
void mapHits4Reset()
{
	_size = 0;
	_numClusters = 0;
	_biggestClusterSize = 0;
	_numBiggestClusters = 0;
}


/**
 * This is a wrapper function that wraps @a mapHits4AddHit function. It has
 * been added so that @a mapHits4AddHit can be unit-tested.
 *
 * @param[out]	hits		Array of hits.
 * @param[out]	size		Number of elements in the @a refIdxArr.
 * @param 		refIdx		Reference index to be added.
 * @param 		shift		Shift to be added.
 * @param 		refPos		Reference position to be added.
 */
void mapHits4AddHit_wrap(ulong **hits, int *size, char refIdx, int shift,
		int refPos)
{
	mapHits4AddHit(refIdx, shift, refPos);
	*hits = _hits;
	*size = _size;

	/* Subtract MAX_QRY_SEQ_LENGTH from the "shift" part of every hit.*/
	int i;
	for (i = 0; i < _size; ++i)
		_hits[i] -= 53687091200; /* 200 << 28 = 53687091200. */
}


/**
 * Adds the given hit to the list of hits.
 *
 * @param refIdx	Reference index.
 * @param shift		Shift (reference position - query position).
 * @param refPos	Reference position.
 */
void mapHits4AddHit(char refIdx, int shift, int refPos)
{
	static ulong hit, tmp;

	hit = (ulong) refIdx;
	hit = (hit << REF_IDX_SHIFT_BITS);
	tmp = (ulong) (shift + MAX_QRY_SEQ_LENGTH);
	tmp = (tmp << SHIFT_SHIFT_BITS);
	hit = (hit | tmp);
	hit += refPos;

	_hits[_size] = hit;
	++_size;
}


/**
 * Returns the best mapped hits.
 *
 * @param 		numBestHits	Number of best hits to be returned.
 * @param[out]	refIdx		Best-matching reference index. The number of
 * elements in this array must be equal to @a numBestHits.
 * @param[out]	shift		Best-matching shift. The number of elements
 * in this array must be equal to @a numBestHits.
 * @param[out]	refPos		Best-matching reference position. The
 * number of elements in this array must be equal to @a numBestHits.
 * @return		Number of best-hits.
 */
int mapHits4GetBestHits(int numBestHits, char *refIdx, int *shift, int *refPos)
{
	/* Return 0 if there are no hits. */
	if (_size == 0)
		return 0;

	/* Sort hits. */
	qsort(_hits, _size, sizeof(ulong), compare_ulong);

	/* Create clusters. */
	createClusters();

	/* Find the biggest cluster size. */
	findBiggestClusterSize();

	/* Find the first hit of the biggest clusters. */
	return findBiggestClusters(numBestHits, refIdx, shift, refPos);
}


/**
 * This is a wrapper function that wraps @a createClusters function. It has
 * been added so that @a createClusters can be unit-tested.
 */
void createClusters_wrap4(int **cluster)
{
	createClusters();
	*cluster = _clusters;
}


/**
 * Creates clusters.
 */
static void createClusters()
{
	_clusters[0] = 0;
	++_numClusters;
	static int i, j;
	static ulong hit1, hit2;
	j = 1;
	hit1 = _hits[0] >> SHIFT_SHIFT_BITS;
	for (i = 1; i < _size; ++i)
	{
//		hit1 = _hits[i - 1] >> SHIFT_SHIFT_BITS;
//		hit2 = _hits[i] >> SHIFT_SHIFT_BITS;
		hit2 = _hits[i] >> SHIFT_SHIFT_BITS;
		/* Ref index and/or shift of both hits are not same. */
		if ((hit1 ^ hit2) != 0)
		{
			_clusters[j] = i;
			++j;
		}
		hit1 = hit2;
	}
	_numClusters = j;
}


/**
 * This is a wrapper function that wraps @a findBiggestClusterSize function.
 * It has been added so that @a findBiggestClusterSize can be unit-tested.
 */
int findBiggestClusterSize_wrap4()
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
		static int i;
		for (i = 1; i < _numClusters; ++i)
		{
			_clusterSize[i - 1] = _clusters[i] - _clusters[i - 1];
			if (_biggestClusterSize < _clusterSize[i - 1])
				_biggestClusterSize = _clusterSize[i - 1];
		}
		_clusterSize[i - 1] = _size - _clusters[i - 1];
		if (_biggestClusterSize < _clusterSize[i - 1])
			_biggestClusterSize = _clusterSize[i - 1];
	}
}


/**
 * This is a wrapper function that wraps @a findBiggestClusters function. It
 * has been added so that @a findBiggestClusters can be unit-tested.
 *
 * @param 		numBestHits			Number of biggest clusters that are to
 * searched for.
 * @param[out]	refIdx_bestMatch	Best-matching reference index. The number
 * of elements in this array must be equal to @a numBestHits.
 * @param[out]	shift_bestMatch		Best-matching shift. The number of elements
 * in this array must be equal to @a numBestHits.
 * @param[out]	refPos_bestMatch	Best-matching reference position. The number
 * of elements in this array must be equal to @a numBestHits.
 * @return 		Number of best-hits in @a refIdx_bestMatch.
 */
int findBiggestClusters_wrap4(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch)
{
	createClusters();
	findBiggestClusterSize();
	return findBiggestClusters(numBestHits, refIdx_bestMatch, shift_bestMatch,
			refPos_bestMatch);
}


/**
 * Finds the biggest clusters.
 *
 * @param 		numBestHits			Number of biggest clusters that are to
 * searched for.
 * @param[out]	refIdx_bestMatch	Best-matching reference index. The number
 * of elements in this array must be equal to @a numBestHits.
 * @param[out]	shift_bestMatch		Best-matching shift. The number of elements
 * in this array must be equal to @a numBestHits.
 * @param[out]	refPos_bestMatch	Best-matching reference position. The number
 * of elements in this array must be equal to @a numBestHits.
 * @return 		Number of best-hits in @a refIdx_bestMatch.
 */
static int findBiggestClusters(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch)
{
	static int i, j;
	j = 0;
	for (i = 0; i < _numClusters; ++i)
	{
		if (_clusterSize[i] == _biggestClusterSize)
		{
			_biggestCluster[j] = _clusters[i];
			++j;
		}
	}
	_numBiggestClusters = j;

	if (numBestHits > _numBiggestClusters)
		numBestHits = _numBiggestClusters;

	arrGetRandomNums((uint) numBestHits, 0, (uint) (_numBiggestClusters - 1),
			(uint *) _randNumArr);
	qsort(_randNumArr, numBestHits, sizeof(int), compare_int);
	static ulong hit;
	for (i = 0; i < numBestHits; ++i)
	{
		hit = _hits[_biggestCluster[_randNumArr[i]]];
		refIdx_bestMatch[i] = (char) (hit >> REF_IDX_SHIFT_BITS);
		shift_bestMatch[i] = (int) (((hit & SHIFT_MASK) >> SHIFT_SHIFT_BITS)
				- MAX_QRY_SEQ_LENGTH);
		refPos_bestMatch[i] = (int) (hit & REF_POS_MASK);
	}
	return numBestHits;
}


/**
 * Returns -1 if @a a is less than @a b, 1 if @a a is greater than @a b, and
 * 0 if they are equal.
 *
 * @param a First parameter.
 * @param b Second parameter.
 * @return -1 if @a a is less than @a b, 1 if @a a is greater than @a b, and
 * 0 if they are equal.
 */
static int compare_ulong(const void *a, const void *b)
{
//	return (*(ulong *)a - *(ulong *)b);
	static ulong x, y;
	x = *(ulong *)a;
	y = *(ulong *)b;
	if (x < y)
		return -1;
	else if (x > y)
		return 1;
	else
		return 0;
}


/**
 * Returns the integer difference between the given two parameters.
 *
 * @param a First parameter.
 * @param b Second parameter.
 * @return The difference between the given two parameters.
 */
static int compare_int(const void *a, const void *b)
{
	return (*(int *)a - *(int *)b);
}


