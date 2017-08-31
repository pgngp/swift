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

#include "refPosMap.h"
#include "common.h"
#include <stdio.h>


/* Array to store reference index and tuple index. Left-most REF_IDX_BITS
 * bits are used to store the reference index and remaining bits are used
 * to store the tuple index. */
static uint *_localPos = NULL;

/* Array to store the file offset corresponding to the reference index
 * and tuple index that are stored in @a _localPos. */
static uint *_globalPos = NULL;

/* Stores the number of elements in each of  @a _localPos and @a _globalPos. */
static int _size = 0;


/**
 * This is a wrapper function that wraps @a refMapCreate function. This
 * function has been added so that @a refMapCreate function can be tested.
 *
 * @param		hashTable		Hashtable using which the map will be created.
 * @param		numKeys			Number of keys in @a hashtable. This number
 * must be equal to the number of keys in @a hashtable, otherwise it will lead
 * to invalid behavior.
 * @param		numVals			Number of values in @a hashtable. This number
 * must be equal to the total number of values in @a hashtable, otherwise it
 * will lead to invalid results.
 * @param		seedLen			Seed length of each reference and query tuple.
 * @param[out]	localPosArr		Array in which reference index and tuple index
 * will be stored.
 * @param[out]	globalPosArr	Array in which the file offset corresponding
 * to the reference index and tuple index in @a localPosArr will be stored.
 * @param[out]	size			Number of elements in each of @a localPosArr and
 * @a globalPosArr.
 */
void refMapCreate_wrap(RefPosList **hashTable, int numKeys, int numVals,
		int seedLen, uint **localPosArr, uint **globalPosArr, int *arrSize)
{
	refMapCreate(hashTable, numKeys, numVals, seedLen);
	*localPosArr = _localPos;
	*globalPosArr = _globalPos;
	*arrSize = _size;
}


/**
 * Creates a table that maps a reference index and tuple index to their file
 * offset.
 *
 * @param	hashtable	Hashtable using which the map will be created.
 * @param	numKeys		Number of keys in @a hashtable. This number must be
 * equal to the number of keys in @a hashtable, otherwise it will lead to
 * invalid behavior.
 * @param	numVals		Number of values in @a hashtable. This number must
 * be equal to the total number of values in @a hashtable, otherwise it will
 * lead to invalid results.
 * @param	seedLen		Seed length of each reference and query tuple.
 */
void refMapCreate(RefPosList **hashTable, int numKeys, int numVals, int seedLen)
{
	if (hashTable == NULL || numKeys == 0 || numVals == 0 || seedLen == 0)
		return;

	_localPos = (uint *) malloc(numVals * sizeof(uint));
	_globalPos = (uint *) malloc(numVals * sizeof(uint));

	/* Copy local and global positions from the hash table to
	 * the parallel arrays. */
	int i, key, index = 0;
	int refIdxShiftBits = (sizeof(int) * NUM_BITS_IN_A_BYTE) - REF_IDX_BITS;
	RefPosList *tmpRefPosList = NULL;
	for (i = 0; i < numKeys; ++i)
	{
		tmpRefPosList = hashTable[i];
		while (tmpRefPosList != NULL)
		{
			key = tmpRefPosList->index;
			key = key << refIdxShiftBits;
			_localPos[index] = key | (tmpRefPosList->position / seedLen);
//			_globalPos[index] = tmpRefPosList->fragOffset;
			++index;
			tmpRefPosList = tmpRefPosList->next;
		}
	}

	/* Sort the parallel arrays. */
	qsort(_localPos, numVals, sizeof(uint), compare);
	qsort(_globalPos, numVals, sizeof(uint), compare);
	_size = numVals;
}


/**
 * Frees resources occupied by @refMapCreate function.
 */
void refMapDelete()
{
	free(_localPos);
	free(_globalPos);
	_localPos = NULL;
	_globalPos = NULL;
	_size = 0;
}


/**
 * Performs integer comparison between the two parameters.
 *
 * @param	a	One of the parameters two be compared.
 * @param	b	Second of the two parameters to be compared.
 * @return		Returns -1 if @a is less than @a b, 0 if they are equal,
 * and 1 if @a is greater than @a b.
 */
static int compare(const void *a, const void *b)
{
//  return (*(uint*)a - *(uint*)b);
	if (*(uint*)a < *(uint*)b)
		return -1;
	else if (*(uint*)a == *(uint*)b)
		return 0;
	else
		return 1;
}


/**
 * Returns a global position (file offset) corresponding to the given
 * local position (reference index and tuple index). Returns UINT_MAX if the
 * global position for the given local position is not available.
 *
 * @param	_localPos	Local position (reference index and tuple index)
 * for which the corresponding global position is to be returned.
 * @return				Global position (file offset) corresponding to
 * the given local position. If the global position for the given local
 * position is not found, then this value will be UINT_MAX.
 */
uint refMapGetGlobalPos(uint local)
{
//	int i;
//	for (i = 0; i < _size && _localPos[i] <= local; ++i)
//	{
//		if (_localPos[i] == local)
//			return _globalPos[i];
//	}
//
//	return UINT_MAX;

	int index = binarySearch(local, _localPos, _size);
	if (index == -1)
		return UINT_MAX;
	else
		return _globalPos[index];
}


/**
 * Prints reference position map.
 */
void refMapPrint()
{
	int i;
	fprintf(stderr, "Local\t->\tGlobal\n");
	for (i = 0; i < _size; ++i)
		fprintf(stderr, "%u\t->\t%u\n", _localPos[i], _globalPos[i]);
}


/**
 * This is a wrapper function that wraps @a binarySearch function. This
 * function has been added so that @a binarySearch function can be unit-tested.
 *
 * @param	x		Integer to be searched in @a arr.
 * @param	arr		Array in which @a x is to be searched.
 * @param	arrSize	Number of elements in @a arr.
 * @return			Array index where @a x appears or -1 if @a x is not found.
 */
int binarySearch_wrap(uint x, uint *arr, uint arrSize)
{
	return binarySearch(x, arr, arrSize);
}


/**
 * Performs a binary search for @a x in @a arr and returns the array index
 * if @a x is found in the array, else returns -1.
 *
 * @param	x		Integer to be searched in @a arr.
 * @param	arr		Array in which @a x is to be searched.
 * @param	arrSize	Number of elements in @a arr.
 * @return			Array index where @a x appears or -1 if @a x is not found.
 *
 * @note This code in this function has been taken from "The C Programming
 * Language" book written by Kernighan and Ritchie.
 */
static int binarySearch(uint x, uint *arr, uint arrSize)
{
	static int low, high, mid;

	low = 0;
	high = arrSize - 1;
	while (low <= high)
	{
		mid = (low + high) / 2;
		if (x < arr[mid])
			high = mid - 1;
		else if (x > arr[mid])
			low = mid + 1;
		else
			return mid;
	}
	return -1;
}
