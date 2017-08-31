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

#include <stdio.h>
#include "refMap.h"
#include "common.h"


static uint *_nameOffset = NULL; /* File-offset of the reference name. */
static uint *_seqOffset = NULL; /* File-offset of the reference sequence. */
static int *_numBasesPerLine = NULL; /* Number of bases per line except
the last line. */
static int _size = 0; /* Number of elements in @a _nameOffset. */


/**
 * This is a wrapper function that wraps @a refMapCreate function. It
 * has been added so that @a refMapCreate can be unit-tested.
 *
 * @param		file				Reference file.
 * @param[out]	nameOffset			Reference name file-offset.
 * @param[out]	seqOffset			Reference sequence file-offset.
 * @param[out]	numBasesPerLine		Number of bases per line except the last
 * line.
 * @param[out]	size				Number of elements in @a nameOffset.
 */
void refMapCreate_wrap(const char *file, uint **nameOffset, uint **seqOffset,
		int **numBasesPerLine, int *size)
{
	refMapCreate(file);
	*nameOffset = _nameOffset;
	*seqOffset = _seqOffset;
	*numBasesPerLine = _numBasesPerLine;
	*size = _size;
}


/**
 * Creates a reference sequence map. Using this map, the user will be
 * able to get the file-offset using the sequence name file-offset and
 * sequence position.
 *
 * @param	file	Reference file.
 */
void refMapCreate(const char *file)
{
	char line[MAX_LINE_LENGTH];
	int lineLength = 0, oldLineLength, numLines;
	uint pos, currentFileOffset;
	FILE *filePtr = fopen(file, "r");
	while (fgets(line, MAX_LINE_LENGTH, filePtr) != NULL)
	{
		if (lineLength > 0)
			oldLineLength = lineLength;
		lineLength = strlen(line);
		if (line[lineLength - 1] == '\n')
			--lineLength;

		/* This is an empty line. */
		if (lineLength == 0)
			continue;
		/* This line contains sequence ID. */
		else if (line[0] == '>')
		{
			++_size;
			_nameOffset = (uint *)
					realloc((uint *) _nameOffset, _size * sizeof(uint));
			_seqOffset = (uint *)
					realloc((uint *) _seqOffset, _size * sizeof(uint));
			_numBasesPerLine = (int *)
					realloc((int *) _numBasesPerLine, _size * sizeof(int));
			currentFileOffset = ftell(filePtr);
			_nameOffset[_size - 1] = currentFileOffset - lineLength - 1;
			_seqOffset[_size - 1] = currentFileOffset;
			if (_size > 1)
			{
				if (numLines > 1)
					_numBasesPerLine[_size - 2] = (pos - oldLineLength)
							/ (numLines - 1);
				else
					_numBasesPerLine[_size - 2] = 0;
			}
			pos = 0;
			numLines = 0;
		}
		/* This is a line containing sequence. */
		else
		{
			++numLines;
			pos += lineLength;
		}
	}
	if (lineLength > 0)
		oldLineLength = lineLength;
	if (_size >= 1)
	{
		if (numLines > 1)
			_numBasesPerLine[_size - 1] = (pos - oldLineLength) / (numLines - 1);
		else
			_numBasesPerLine[_size - 1] = 0;
	}
	fclose(filePtr);
}


/**
 * Returns a sequence file-offset using the name file-offset and sequence
 * position. @a nameOffset and @a pos should be valid, otherwise this function
 * will return an invalid value.
 *
 * @param	nameOffset	Reference name file-offset.
 * @param	pos			Reference sequence position.
 */
uint refMapGetFileOffset(uint nameOffset, uint pos)
{
	int index = binarySearch(nameOffset, _nameOffset, _size);
	int lineNumber = 0;
	if (index >= 0 && _numBasesPerLine[index] > 0)
		lineNumber = pos / _numBasesPerLine[index];
	return (_seqOffset[index] + pos + lineNumber);
}


/**
 * Fetches reference name and sequence file-offset using the reference index
 * and sequence position. @a index and @a pos should be valid, otherwise this
 * function will return an invalid value.
 *
 * @param       index       Reference index.
 * @param       pos         Reference sequence position.
 * @param[out]  nameOffset  File-offset of the reference name.
 * @param[out]  seqOffset   File-offset of the reference sequence fragment.
 */
void refMapGetOffsets(int index, uint pos, uint *nameOffset, uint *seqOffset)
{
    *nameOffset = _nameOffset[index];
    int lineNumber = 0;
    if (_numBasesPerLine[index] > 0)
        lineNumber = pos / _numBasesPerLine[index];
    *seqOffset = (_seqOffset[index] + pos + lineNumber);
}


/**
 * Fetches reference index and position of the sequence fragment using the given
 * file-offset.
 *
 * @param		seqOffset	File-offset of a sequence fragment.
 * @param[out] 	index		Reference index of the sequence fragment.
 * @param[out]	pos			Position of the sequence fragment.
 */
void refMapGetIndexAndPos(uint seqOffset, int *index, int *pos)
{
	static int i, numNewLines;
	static uint tmpSeqOffset;
	*index = -1;
	*pos = -1;

	if (_seqOffset[0] == seqOffset)
	{
		*index = 0;
		*pos = 0;
		return;
	}
	else if (_size == 1)
	{
		*index = 0;
		if (_numBasesPerLine[*index] == 0)
			*pos = seqOffset - _seqOffset[*index];
		else
		{
			numNewLines = 0;
			tmpSeqOffset = _seqOffset[*index];
			while ((tmpSeqOffset + _numBasesPerLine[*index] + 1) <= seqOffset)
			{
				tmpSeqOffset += _numBasesPerLine[*index] + 1;
				++numNewLines;
			}
			*pos = seqOffset - _seqOffset[*index] - numNewLines;
		}
		return;
	}

	for (i = 1; i < _size; ++i)
	{
		if (_seqOffset[i - 1] == seqOffset)
		{
			*index = i - 1;
			*pos = 0;
			return;
		}
		else if (_seqOffset[i] > seqOffset)
		{
			*index = i - 1;
			if (_numBasesPerLine[*index] == 0)
				*pos = seqOffset - _seqOffset[*index];
			else
			{
				numNewLines = 0;
				tmpSeqOffset = _seqOffset[*index];
				while ((tmpSeqOffset + _numBasesPerLine[*index] + 1) <= seqOffset)
				{
					tmpSeqOffset += _numBasesPerLine[*index] + 1;
					++numNewLines;
				}
				*pos = seqOffset - _seqOffset[i - 1] - numNewLines;
			}
			return;
		}
	}
	if (*index == -1)
	{
		*index = _size - 1;
		if (_numBasesPerLine[*index] == 0)
			*pos = seqOffset - _seqOffset[*index];
		else
		{
			numNewLines = 0;
			tmpSeqOffset = _seqOffset[*index];
			while ((tmpSeqOffset + _numBasesPerLine[*index] + 1) <= seqOffset)
			{
				tmpSeqOffset += _numBasesPerLine[*index] + 1;
				++numNewLines;
			}
			*pos = seqOffset - _seqOffset[*index] - numNewLines;
		}
	}
}


/**
 * Releases resources used by reference map.
 */
void refMapFree()
{
	free(_nameOffset);
	free(_seqOffset);
	free(_numBasesPerLine);

	_nameOffset = NULL;
	_seqOffset = NULL;
	_numBasesPerLine = NULL;
	_size = 0;
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
int binarySearch_wrapper(uint x, uint *arr, uint arrSize)
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


/**
 * Prints the reference map.
 */
void refMapPrint2()
{
	fprintf(stderr, "\nReference map:\n");
	int i;
	for (i = 0; i < _size; ++i)
	{
		fprintf(stderr, "nameOffset[%d] = %u, seqOffset[%d] = %u\n",
				i, _nameOffset[i], i, _seqOffset[i]);
	}
}

