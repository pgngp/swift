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

#include "refNameMap.h"
#include <stdio.h>
#include "common.h"

/* Array where reference name file offset will be stored. */
static int *_namePos = NULL;

/* Array where reference sequence-start file offset will be stored. */
static int *_seqStartPos = NULL;

/* Array where reference sequence-end file offset will be stored. */
static int *_seqEndPos = NULL;

/* Variable where number of elements in each of @a _namePos, @a _seqStartPos,
 * and @a _seqEndPos will be stored. */
static int _size = 0;


/**
 * This is a wrapper function that wraps @a refNameMapCreate function. This
 * function has been added so that @a refNameMapCreate can be unit-tested.
 *
 * @param		refFile		Reference file path.
 * @param[out]	namePos		Array where reference name file offset will be
 * stored.
 * @param[out]	seqStartPos	Array where reference sequence-start file offset
 * will be stored.
 * @param[out]	seqEndPos	Array where reference sequence-end file offset
 * will be stored.
 * @param[out]	size		Memory location where number of elements in each of
 * @a namePos, @a seqStartPos, and @a seqEndPos will be stored.
 */
void refNameMapCreate_wrap(const char *refFile, int **namePos,
		int **seqStartPos, int **seqEndPos, int *size)
{
	refNameMapCreate(refFile);
	*namePos = _namePos;
	*seqStartPos = _seqStartPos;
	*seqEndPos = _seqEndPos;
	*size = _size;
}


/**
 * Creates a map that maps each reference name file offset with that
 * reference's sequence start file offset and sequence end file offset.
 *
 * @note In order to de-allocate resources used in this file, you will need
 * to call @a refNameMapDelete function after you are done using the
 * functions in this file.
 *
 * @param	refFile		Reference file path.
 */
void refNameMapCreate(const char *refFile)
{
	if (refFile == NULL)
		return;

	FILE *filePtr = fopen(refFile, "r");
	char line[MAX_LINE_LENGTH];
	int lineLength, currentFileOffset;
	while (fgets(line, MAX_LINE_LENGTH, filePtr) != NULL)
	{
		if (line[0] == '>')
		{
			_namePos = (int *) realloc(_namePos, (_size + 1) * sizeof(int));
			_seqStartPos = (int *) realloc(_seqStartPos,
					(_size + 1) * sizeof(int));
			_seqEndPos = (int *) realloc(_seqEndPos, (_size + 1) * sizeof(int));

			lineLength = strlen(line);
			currentFileOffset = ftell(filePtr);
			_seqStartPos[_size] = currentFileOffset;
			_namePos[_size] = currentFileOffset - lineLength;
			if (_size > 0)
				_seqEndPos[_size - 1] = currentFileOffset - lineLength - 1;
			++_size;
		}
	}
	_seqEndPos[_size - 1] = ftell(filePtr) - 1;
	fclose(filePtr);
}


/**
 * Returns file offset of the name of the reference in which the given
 * sequence file offset belongs. Returns -1 if the file offset of the
 * corresponding reference name is not found.
 *
 * @note Before this function is called, you need to first initialize the
 * data structures in this file by calling @a refNameMapCreate function.
 *
 * @param	seqOffset	File offset of a sequence fragment.
 * @return				File offset of the name of the reference in which the
 * given sequence file offset belongs. If no corresponding name file offset
 * is found, -1 will be returned.
 */
int refNameMapGetNameOffset(uint seqOffset)
{
	int i;
	for (i = 0; i < _size && seqOffset > _namePos[i]; ++i)
	{
		if (seqOffset >= _seqStartPos[i] && seqOffset <= _seqEndPos[i])
			return _namePos[i];
	}

	return -1;
}


/**
 * Frees resrouces allocated by @a refNameMapCreate function.
 */
void refNameMapDelete()
{
	free(_namePos);
	free(_seqStartPos);
	free(_seqEndPos);
	_namePos = NULL;
	_seqStartPos = NULL;
	_seqEndPos = NULL;
	_size = 0;
}


/**
 * Prints reference name map.
 */
void refNameMapPrint()
{
	int i;
	fprintf(stderr, "NamePos\tSeqStart\tSeqEnd\n");
	for (i = 0; i < _size; ++i)
		fprintf(stderr, "%d\t%d\t%d\n", _namePos[i], _seqStartPos[i],
				_seqEndPos[i]);
}
