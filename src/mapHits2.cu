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

#include "mapHits2.h"
#include "array.h"
#include <stdio.h>
#include <limits.h>

#define MAX_NUM_CHR	32	/* Max number of chromosomes. */


static int _seedLen = 0; /* Seed length. */
static char *_refIdx = NULL; /* Reference index */
static int *_shift = NULL; /* Shift */
static int *_refPos = NULL; /* Reference position */
static int *_order = NULL; /* Maintains the sorted order of @a _refIdx,
@a _shift, and @a _refPos. */
static int _size = 0; /* Size of the current mapped list. */
static int *_numClusterHits = NULL; /* Number of hits in each cluster. */
static int *_randNumArr = NULL; /* Random number array. */
static int *_biggestCluster = NULL; /* Stores index of the biggest clusters. */
static int _biggestClusterSize = 0; /* Biggest cluster size. */
static int _numBiggestClusters = 0; /* Number of biggest clusters. */
static uchar *_hits = NULL; /* Keeps track of the number of
hits with the same reference index and shift. */
static char *_ref = NULL; /* Stores the reference index. */
static int *_position = NULL; /* Stores the smallest reference position
for any given reference index and shift. */
static int *_lastPos = NULL; /* Stores the last position of each reference. */
static int _maxNumShifts = 0; /* Max number of shifts possible. */



/**
 * Creates a list of mapped hits.
 *
 * @param seedLen	Seed length.
 */
void mapHits2Create(int seedLen)
{
	_seedLen = seedLen;
	int numQryTuples = MAX_QRY_SEQ_LENGTH - seedLen;
	int numElements = numQryTuples * TUPLE_IGNORE_THRES;
	_refIdx = (char *) calloc(numElements, sizeof(char));
	_shift = (int *) calloc(numElements, sizeof(int));
	_refPos = (int *) calloc(numElements, sizeof(int));
	_order = (int *) calloc(numElements, sizeof(int));
	_numClusterHits = (int *) calloc(numElements, sizeof(int));
	_randNumArr = (int *) calloc(numElements, sizeof(int));
	_biggestCluster = (int *) calloc(numElements, sizeof(int));
}


/**
 * Creates a list of mapped hits.
 *
 * @param refFile	Reference file.
 * @param seedLen	Seed length.
 */
void mapHits2Create2(const char *refFile, int seedLen)
{
	_seedLen = seedLen;
	int numQryTuples = MAX_QRY_SEQ_LENGTH - seedLen;
	int numElements = numQryTuples * TUPLE_IGNORE_THRES;
	_randNumArr = (int *) calloc(numElements, sizeof(int));
	_biggestCluster = (int *) calloc(numElements, sizeof(int));
	/* I'm adding by MAX_QRY_SEQ_LENGTH to account for that many negative
	 * shifts. */
	_maxNumShifts = (MAX_GENOMIC_SIZE / seedLen) + MAX_QRY_SEQ_LENGTH;
	_hits = (uchar *) calloc(_maxNumShifts, sizeof(uchar));
	_ref = (char *) calloc(_maxNumShifts, sizeof(char *));
	_position = (int *) calloc(_maxNumShifts, sizeof(int));
	_lastPos = (int *) calloc(MAX_NUM_CHR, sizeof(int));
	createRefMap(refFile);

	/* Initialize _position. */
	int i;
	for (i = 0; i < _maxNumShifts; ++i)
		_position[i] = INT_MAX;
}


/**
 * This is a wrapper function that wraps @a createRefMap function. It has
 * been added so that @a createRefMap can be unit-tested.
 *
 * @param refFile	Reference file.
 * @param lastPos	Last position of each reference.
 */
void createRefMap_wrap(int **lastPos)
{
	*lastPos = _lastPos;
}


/**
 * Creates a map of the reference that will be used by other functions
 * in this file.
 *
 * @param refFile	Reference file.
 */
static void createRefMap(const char *refFile)
{
	char line[MAX_LINE_LENGTH];
	int lineLength, refIdx = 0, pos = 0, numBases = 0;
	FILE *filePtr = fopen(refFile, "r");
	while (fgets(line, MAX_LINE_LENGTH, filePtr) != NULL)
	{
		lineLength = strlen(line);
		if (line[lineLength - 1] == '\n')
			--lineLength;

		/* This is an empty line. */
		if (lineLength == 0)
			continue;
		/* This line contains sequence ID. */
		else if (line[0] == '>')
		{
			pos += (numBases / _seedLen);
			_lastPos[refIdx] = pos;
			++refIdx;
			numBases = 0;
		}
		else
			numBases += lineLength;
	}
	pos += (numBases / _seedLen);
	_lastPos[refIdx] = pos;
	fclose(filePtr);
}


/**
 * Reset key variables in this file.
 */
void mapHits2Reset()
{
	_size = 0;
	_biggestClusterSize = 0;
	_numBiggestClusters = 0;

	static int i;
	for (i = 0; i < _maxNumShifts; ++i)
	{
		_ref[i] = 0;
		_hits[i] = 0;
		_position[i] = INT_MAX;
	}
}


/**
 * Deletes mapped hits.
 */
void mapHits2Delete()
{
	free(_refIdx);
	_refIdx = NULL;
	free(_shift);
	_shift = NULL;
	free(_refPos);
	_refPos = NULL;
	free(_numClusterHits);
	_numClusterHits = NULL;
	free(_randNumArr);
	_randNumArr = NULL;
	free(_biggestCluster);
	_biggestCluster = NULL;
	_size = 0;
	_biggestClusterSize = 0;
	_numBiggestClusters = 0;
}


/**
 * Deletes mapped hits.
 */
void mapHits2Delete2()
{
	free(_refIdx);
	_refIdx = NULL;
	free(_shift);
	_shift = NULL;
	free(_refPos);
	_refPos = NULL;
	free(_numClusterHits);
	_numClusterHits = NULL;
	free(_randNumArr);
	_randNumArr = NULL;
	free(_biggestCluster);
	_biggestCluster = NULL;
	free(_hits);
	_hits = NULL;
	free(_position);
	_position = NULL;
	free(_lastPos);
	_lastPos = NULL;
	free(_ref);
	_ref = NULL;
	_size = 0;
	_biggestClusterSize = 0;
	_numBiggestClusters = 0;
}


/**
 * This is a wrapper function that wraps @a mapHits2AddHit function. It has
 * been added so that @a mapHits2AddHit can be unit-tested.
 *
 * @param[out]	refIdxArr			Reference index array.
 * @param[out] 	shiftArr			Shift array.
 * @param[out]	refPosArr			Reference position array.
 * @param[out]	numClusterHits		Number of hits in each cluster.
 * @param[out]	size				Number of elements in the @a refIdxArr.
 * @param 		refIdx				Reference index to be added.
 * @param 		shift				Shift to be added.
 * @param 		refPos				Reference position to be added.
 */
void mapHits2AddHit_wrap(char **refIdxArr, int **shiftArr, int **refPosArr,
		int **numClusterHits, int *size, char refIdx,
		int shift, int refPos)
{
	mapHits2AddHit(refIdx, shift, refPos);
	*refIdxArr = _refIdx;
	*shiftArr = _shift;
	*refPosArr = _refPos;
	*numClusterHits = _numClusterHits;
	*size = _size;
}


/**
 * Adds the given hit to the list.
 *
 * @param refIdx	Reference index.
 * @param shift		Shift (reference position - query position).
 * @param refPos	Reference position.
 */
void mapHits2AddHit(char refIdx, int shift, int refPos)
{
	static int nextHighestIdx, index, i;
	nextHighestIdx = 0;
	index = search(refIdx, shift, &nextHighestIdx);
	if (index == -1)
	{
		_refIdx[_size] = refIdx;
		_shift[_size] = shift;
		_refPos[_size] = refPos;
		_numClusterHits[_size] = 1;
		for (i = _size; i > nextHighestIdx; --i)
			_order[i] = _order[i - 1];
		_order[nextHighestIdx] = _size;
		++_size;
	}
	else
	{
		++_numClusterHits[index];
		if (_refPos[index] > refPos)
			_refPos[index] = refPos;
	}
}


/**
 * This is a wrapper function that wraps @a mapHits2AddHit function. It has
 * been added so that @a mapHits2AddHit can be unit-tested.
 *
 * @param[out]	ref					Reference indexes.
 * @param[out] 	hits				Number of hits in each cluster.
 * @param[out]	pos					Reference positions.
 * @param[out]	biggestClusterSize	Biggest cluster size.
 * @param 		refIdx				Reference index to be added.
 * @param 		qryPos				Query position.
 * @param 		refPos				Reference position to be added.
 */
void mapHits2AddHit2_wrap(char **ref, uchar **hits, int **pos,
		int *biggestClusterSize, char refIdx, int qryPos, int refPos)
{
	mapHits2AddHit2(refIdx, qryPos, refPos);
	*ref = _ref;
	*hits = _hits;
	*pos = _position;
	*biggestClusterSize = _biggestClusterSize;
}


/**
 * Adds the given hit to the list.
 *
 * @param refIdx	Reference index.
 * @param qryPos	Query position.
 * @param refPos	Reference position.
 */
void mapHits2AddHit2(char refIdx, int qryPos, int refPos)
{
	static int pos;
	pos = (refPos + 1) / _seedLen;
	pos += _lastPos[refIdx];
	pos = pos - ((qryPos + 1) / _seedLen) + MAX_QRY_SEQ_LENGTH; /* Adding by
	MAX_QRY_SEQ_LENGTH to move negative shifts into positive. */
	_ref[pos] = refIdx;
	++_hits[pos];
	if (_biggestClusterSize < _hits[pos])
		_biggestClusterSize = _hits[pos];
	if (_position[pos] > refPos)
		_position[pos] = refPos;
}


/**
 * This is a wrapper function that wraps @a search function. It has been
 * added so that @a wrap can be unit-tested.
 *
 * @param		refIdx			Reference index to be searched for.
 * @param		shift			'Shift' to be searched for.
 * @param[out]	nextHighestIdx	Index of the next highest pair of reference
 * index and shift.
 * @return		Returns the index of the element corresponding to the reference
 * index and shift, if found; otherwise, returns -1.
 */
int search_wrap(char refIdx, int shift, int *nextHighestIdx)
{
	return search(refIdx, shift, nextHighestIdx);
}


/**
 * Searches for the given reference index and shift in the list. If found,
 * returns the index of the element associated with the reference index and
 * shift; if not found, returns -1. This function uses the binary search
 * algorithm.
 *
 * @param 		refIdx			Reference index to be searched for.
 * @param 		shift			'Shift' to be searched for.
 * @param[out]	nextHighestIdx	Index of the next highest pair of reference
 * index and shift.
 * @return		Returns the index of the element corresponding to the reference
 * index and shift, if found; otherwise, returns -1.
 */
static int search(char refIdx, int shift, int *nextHighestIdx)
{
	static int low, high, mid, resolvedIdx;
	low = 0;
	high = _size - 1;
	resolvedIdx = -1;
	while (low <= high)
	{
		mid = (low + high) / 2;
		resolvedIdx = _order[mid];
		if (refIdx < _refIdx[resolvedIdx]
				|| ((refIdx == _refIdx[resolvedIdx])
						&& (shift < _shift[resolvedIdx])))
			high = mid - 1;
		else if (refIdx > _refIdx[resolvedIdx]
				|| ((refIdx == _refIdx[resolvedIdx])
						&& (shift > _shift[resolvedIdx])))
			low = mid + 1;
		else
			return resolvedIdx;
	}

	if (resolvedIdx > -1 && (refIdx > _refIdx[resolvedIdx]
			|| (refIdx == _refIdx[resolvedIdx] && shift > _shift[resolvedIdx])))
		*nextHighestIdx = mid + 1;

	return -1;
}


/**
 * Returns the best mapped hits.
 *
 * @param 		numBestHits			Number of best hits to be returned.
 * @param[out]	refIdx_bestMatch	Best-matching reference index. The number of
 * elements in this array must be equal to @a numBestHits.
 * @param[out]	shift_bestMatch		Best-matching shift. The number of elements
 * in this array must be equal to @a numBestHits.
 * @param[out]	refPos_bestMatch	Best-matching reference position. The
 * number of elements in this array must be equal to @a numBestHits.
 * @return		Number of best-hits.
 */
int mapHits2GetBestHits(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch)
{
	/* Return 0 if there are no hits. */
	if (_size == 0)
		return 0;

	/* Find the biggest cluster size. */
	findBiggestClusterSize();

	/* Find the first hit of the biggest clusters. */
	return findBiggestClusters(numBestHits, refIdx_bestMatch, shift_bestMatch,
			refPos_bestMatch);
}


/**
 * Returns the best mapped hits.
 *
 * @param 		numBestHits			Number of best hits to be returned.
 * @param[out]	refIdx_bestMatch	Best-matching reference index. The number of
 * elements in this array must be equal to @a numBestHits.
 * @param[out]	refPos_bestMatch	Best-matching reference position. The
 * number of elements in this array must be equal to @a numBestHits.
 * @return		Number of best-hits.
 */
int mapHits2GetBestHits2(int numBestHits, char *refIdx_bestMatch,
		int *refPos_bestMatch)
{
	/* Return 0 if there are no hits. */
	if (_size == 0)
		return 0;

	/* Find the first hit of the biggest clusters. */
	return findBiggestClusters2(numBestHits, refIdx_bestMatch, refPos_bestMatch);
}


/**
 * This is a wrapper function that wraps @a findBiggestClusterSize function.
 * It has been added so that @a findBiggestClusterSize can be unit-tested.
 */
int findBiggestClusterSize_wrap2()
{
	findBiggestClusterSize();
	return _biggestClusterSize;
}


/**
 * Find the biggest cluster size.
 */
static void findBiggestClusterSize()
{
	if (_size == 1)
		_biggestClusterSize = _numClusterHits[0];
	else
	{
		int i, max = 0;
		for (i = 0; i < _size; ++i)
		{
			if (max < _numClusterHits[i])
				max = _numClusterHits[i];
		}
		_biggestClusterSize = max;
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
int findBiggestClusters_wrap2(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch)
{
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
	int i, j = 0;
	for (i = 0; i < _size; ++i)
	{
		if (_numClusterHits[i] == _biggestClusterSize)
		{
			_biggestCluster[j] = i;
			++j;
		}
	}
	_numBiggestClusters = j;

	if (numBestHits > _numBiggestClusters)
		numBestHits = _numBiggestClusters;

	arrGetRandomNums((uint) numBestHits, 0, (uint) (_numBiggestClusters - 1),
			(uint *) _randNumArr);
	qsort(_randNumArr, numBestHits, sizeof(int), compare);
	int tmp;
	for (i = 0; i < numBestHits; ++i)
	{
		tmp = _biggestCluster[_randNumArr[i]];
		refIdx_bestMatch[i] = _refIdx[tmp];
		shift_bestMatch[i] = _shift[tmp];
		refPos_bestMatch[i] = _position[tmp];
	}
	return numBestHits;
}


/**
 * This is a wrapper function that wraps @a findBiggestClusters function. It
 * has been added so that @a findBiggestClusters can be unit-tested.
 *
 * @param 		numBestHits			Number of biggest clusters that are to
 * searched for.
 * @param[out]	refIdx_bestMatch	Best-matching reference index. The number
 * of elements in this array must be equal to @a numBestHits.
 * @param[out]	refPos_bestMatch	Best-matching reference position. The number
 * of elements in this array must be equal to @a numBestHits.
 * @return 		Number of best-hits in @a refIdx_bestMatch.
 */
int findBiggestClusters2_wrap2(int numBestHits, char *refIdx_bestMatch,
		int *refPos_bestMatch)
{
	return findBiggestClusters2(numBestHits, refIdx_bestMatch, refPos_bestMatch);
}


/**
 * Finds the biggest clusters.
 *
 * @param 		numBestHits			Number of biggest clusters that are to
 * searched for.
 * @param[out]	refIdx_bestMatch	Best-matching reference index. The number
 * of elements in this array must be equal to @a numBestHits.
 * @param[out]	refPos_bestMatch	Best-matching reference position. The number
 * of elements in this array must be equal to @a numBestHits.
 * @return 		Number of best-hits in @a refIdx_bestMatch.
 */
static int findBiggestClusters2(int numBestHits, char *refIdx_bestMatch,
		int *refPos_bestMatch)
{
	static int i, j, tmp;
	j = 0;
	for (i = 0; i < _maxNumShifts; ++i)
	{
		if (_hits[i] == _biggestClusterSize)
		{
			_biggestCluster[j] = i;
			++j;
		}
	}
	_numBiggestClusters = j;

	if (numBestHits > _numBiggestClusters)
		numBestHits = _numBiggestClusters;

	arrGetRandomNums((uint) numBestHits, 0, (uint) (_numBiggestClusters - 1),
			(uint *) _randNumArr);
	qsort(_randNumArr, numBestHits, sizeof(int), compare);
	for (i = 0; i < numBestHits; ++i)
	{
		tmp = _biggestCluster[_randNumArr[i]];
		refIdx_bestMatch[i] = _ref[tmp];
		refPos_bestMatch[i] = _position[tmp];
	}
	return numBestHits;
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

