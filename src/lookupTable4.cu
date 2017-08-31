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
#include "lookupTable4.h"
#include "preprocess.h"
#include "mapHits.h"
#include "mapHits2.h"
#include "mapHits3.h"
#include "mapHits4.h"
#include "mapHits5.h"
#include <stdio.h>
#include <limits.h>
#include <time.h>


static int *_lookupTable = NULL;
static char *_refIdx = NULL;
static int *_refPos = NULL;
static int _seedLen = 0;
static int _maxHitsPerQry = 0;
static int _numDistinctTuples = 0;
static int *_numRepeatsPerTuple = NULL;
static char *_seed = NULL;
static int _tupleIgnoreThreshold = 0;


/**
 * This is a wrapper function that wraps @a lookupTable3Create function.
 * This function has been added so that @a lookupTable3Create can be
 * unit-tested.
 *
 * @param		refFile					Reference file.
 * @param		seedLen					Seed length.
 * @param 		maxHitsPerQry			Maximum hits per query.
 * @param[out]	lookupTable				Lookup table.
 * @param[out]	refIdx					Reference index.
 * @param[out]	refPos					Reference sequence offset.
 * @param[out]	numDistinctTuples		Total number of distinct tuples.
 * @param[out]	numRepeatsPerTuple		Number of repeats per tuple.
 * @param		tupleIgnoreThreshold	Tuple that have number of repeats higher
 * than this value will be ignored.
 */
void lookupTable4Create_wrap(const char *refFile, int seedLen,
		int maxHitsPerQry, int **lookupTable, char **refIdx, int **refPos,
		int *numDistinctTuples, int **numRepeatsPerTuple,
		int tupleIgnoreThreshold)
{
	lookupTable4Create(refFile, seedLen, maxHitsPerQry, tupleIgnoreThreshold);
	*lookupTable = _lookupTable;
	*refIdx = _refIdx;
	*refPos = _refPos;
	*numDistinctTuples = _numDistinctTuples;
	*numRepeatsPerTuple = _numRepeatsPerTuple;
}


/**
 * Creates a lookup table.
 *
 * @param	refFile					Reference file.
 * @param	seedLen					Seed length.
 * @param 	maxHitsPerQry			Maximum hits per query.
 * @param	tupleIgnoreThreshold	Tuple that have number of repeats higher
 * than this value will be ignored.
 */
void lookupTable4Create(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold)
{
	_seedLen = seedLen;
	_maxHitsPerQry = maxHitsPerQry;
	_tupleIgnoreThreshold = tupleIgnoreThreshold;
	_numDistinctTuples = (int) pow((float) DNA_ALPHABET_SIZE, (int) _seedLen);
	_lookupTable = (int *) calloc(_numDistinctTuples, sizeof(int));
	int *numRepeatsPerTuple = (int *) calloc(_numDistinctTuples, sizeof(int));
	_numRepeatsPerTuple = (int *) calloc(_numDistinctTuples, sizeof(int));
	_seed = (char *) calloc((_seedLen + 1), sizeof(char));

	/* First pass through reference to find the number of repeats for each
	 * distinct tuple. */
	time_t startTime, endTime;
	double diffTime;
	time(&startTime);
	fprintf(stderr, "   First pass...");
	FILE *filePtr = fopen(refFile, "r");
	int numTotalTuples = 0, numIterations, i, j, offset = 0;
	int hashVal, lineLength, numBases = 0, refIdx = -1;
	char line[MAX_LINE_LENGTH];

	char tmpRefFile[MAX_FILE_NAME_LENGTH];
	sprintf(tmpRefFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, TEMP_REF_FILE);
	FILE *tmpRefFilePtr = fopen(tmpRefFile, "w+");
	while (fgets(line + offset, MAX_LINE_LENGTH, filePtr) != NULL)
	{
		lineLength = strlen(line);
		if (line[lineLength - 1] == '\n')
		{
			line[lineLength - 1] = '\0';
			--lineLength;
		}

		/* This is an empty line. */
		if (lineLength == 0)
			offset = 0;
		/* This line contains reference ID. */
		else if (line[offset] == '>')
		{
			offset = 0;
			numBases = 0;
			++refIdx;
		}
		/* This is a line containing sequence. */
		else
		{
			numIterations = lineLength - _seedLen + 1;
			/* Consider non-overlapping tuples only. */
			for (i = 0; i < numIterations; i += _seedLen)
			{
				hashVal = getHash(line + i, _seedLen);
				++_numRepeatsPerTuple[hashVal];
				if (_numRepeatsPerTuple[hashVal] <= _tupleIgnoreThreshold)
				{
					++numRepeatsPerTuple[hashVal];
					++numTotalTuples;
					fprintf(tmpRefFilePtr, "%d\t%d\t%d\n", hashVal, refIdx,
							numBases + i);
				}
				else if (_numRepeatsPerTuple[hashVal]
				                            == numRepeatsPerTuple[hashVal] + 1)
				{
					numRepeatsPerTuple[hashVal] = 0;
					numTotalTuples -= _tupleIgnoreThreshold;
				}
			}
			numBases += lineLength;

			/* Copy the last few bases to the beginning of 'line' array */
			offset = lineLength - i;
			for (j = 0; j < offset; ++j)
				line[j] = line[i + j];
			numBases -= offset;
		}
	}
	_refPos = (int *) calloc((numTotalTuples + 1), sizeof(int));
	_refPos[0] = -1; /* First element set to -1 so that
	tuples that do not exist in the reference or that have number of repeats
	greater than a threshold value can point to this element. */
	_refIdx = (char *) calloc((numTotalTuples + 1), sizeof(int));
	_refIdx[0] = -1;
	time(&endTime);
	diffTime = difftime(endTime, startTime);
	fprintf(stderr, "done. (Time = %.2lf secs)\n", diffTime);

	/* Set values in the lookup table. */
	time(&startTime);
	fprintf(stderr, "   Set values in lookup table...");
	int start = 1;
	for (i = 0; i < _numDistinctTuples; ++i)
	{
		if (_numRepeatsPerTuple[i] > 0
				&& _numRepeatsPerTuple[i] <= tupleIgnoreThreshold)
		{
			_lookupTable[i] = start;
			start += _numRepeatsPerTuple[i];
		}
	}
	time(&endTime);
	diffTime = difftime(endTime, startTime);
	fprintf(stderr, "done. (Time = %.2lf secs)\n", diffTime);

	/* Set reference offsets. */
	time(&startTime);
	fprintf(stderr, "   Setting reference positions...");
	rewind(tmpRefFilePtr);
	int index, refPos;
	while (fgets(line, MAX_LINE_LENGTH, tmpRefFilePtr) != NULL)
	{
		sscanf(line, "%d\t%d\t%d", &hashVal, &refIdx, &refPos);
		if (_numRepeatsPerTuple[hashVal] <= _tupleIgnoreThreshold)
		{
			index = _lookupTable[hashVal] + numRepeatsPerTuple[hashVal] - 1;
			_refIdx[index] = (char) refIdx;
			_refPos[index] = refPos;
			--numRepeatsPerTuple[hashVal];
		}
	}
	time(&endTime);
	diffTime = difftime(endTime, startTime);
	fprintf(stderr, "done. (Time = %.2lf secs)\n", diffTime);

	fprintf(stderr, "   Freeing resources...");
	fclose(tmpRefFilePtr);
	remove(tmpRefFile);
	fclose(filePtr);
	free(numRepeatsPerTuple);
	fprintf(stderr, "done.\n");

	/* Create mapped-hits data structures. */
	mapHitsCreate(_seedLen);
//	mapHits5Create(_seedLen);
}


/**
 * Releases memory occupied by data structures in this file.
 */
void lookupTable4Delete()
{
	free(_lookupTable);
	_lookupTable = NULL;
	free(_refPos);
	_refPos = NULL;
	free(_refIdx);
	_refIdx = NULL;
	_seedLen = 0;
	_maxHitsPerQry = 0;
	_numDistinctTuples = 0;
	free(_seed);
	_seed = NULL;
	free(_numRepeatsPerTuple);
	_numRepeatsPerTuple = NULL;

	/* Delete mapped-hits data structures. */
	mapHitsDelete();
//	mapHits5Delete();
}


/**
 * Searches the query sequence in the reference and returns the best-matching
 * reference coordinates.
 *
 * @param	query 		Query object.
 * @param 	seqLength 	Sequence length.
 * @param 	isRevComp 	A value of '1' indicates that the
 * reverse complement of the query should be used for search; '0' indicates
 * otherwise.
 * @param	refIdx		Best-matching reference indexes.
 * @param	shift		Best-matching shifts.
 * @param	refPos		Best-matching reference positions.
 * @return 	Number of hits in @a refIdx.
 */
int lookupTable4MapQry(const Query *query, uint qryLen, uint isRevComp,
		char *refIdx, int *shift, int *refPos)
{
	static int numQryTuples, i, j, hash, startIdx, endIdx, numHits;
	static char *seq;

	if (isRevComp == 1)
		seq = query->revCompSeq;
	else
		seq = query->seq;

	/* Reset mapped-hits. */
	mapHitsReset();
//	mapHits5Reset();

	/* Step 1: Break the query into tuples. */
	numQryTuples = qryLen - _seedLen + 1;
	for (i = 0; i < numQryTuples; ++i)
	{
		/* Step 2: Calculate the hash for each tuple. */
		/* Step 3: Find the reference coordinates of the hash in the hash
		 * table. */
		/* Step 4: Aggregate the reference coordinate positions. */
//		strncpy(_seed, seq + i, _seedLen);
//		hash = getHash(_seed, _seedLen);
		hash = getHash(seq + i, _seedLen);
		if (_numRepeatsPerTuple[hash] > _tupleIgnoreThreshold)
			continue;

		/* Add the tuples to the list. */
		startIdx = _lookupTable[hash];
		endIdx = startIdx + _numRepeatsPerTuple[hash] - 1;
		for (j = startIdx; j <= endIdx; ++j)
			mapHitsAddHit(_refIdx[j], _refPos[j] - i, _refPos[j]);
//			mapHits5AddHit(_refIdx[j], _refPos[j]);
	}
	numHits = mapHitsGetBestHits(_maxHitsPerQry, refIdx, shift, refPos);
//	numHits = mapHits5GetBestHits(_maxHitsPerQry, refIdx, refPos);

	return numHits;
}

