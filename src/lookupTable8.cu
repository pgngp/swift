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

/**
 * @file
 *
 * @section DESCRIPTION
 *
 * Implements CPU-based filtering; uses linked list.
 */

#include "common.h"
#include "lookupTable8.h"
#include "preprocess.h"
#include "mapHits6.h"
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include "array.h"


#define	REF_POS_BITS2	25
#define	REF_POS_MASK2	33554431 /* Binary: 1111111111111111111111111 */
#define	TOTAL_NUM_CHRS	25

static int *_lookupTable = NULL;
static char *_refIdx = NULL;
static int *_refPos = NULL;
static int *_refPos2 = NULL;
static int _seedLen = 0;
static int _maxHitsPerQry = 0;
static int _numDistinctTuples = 0;
static int *_numRepeatsPerTuple = NULL;
static char *_seed = NULL;
static int _tupleIgnoreThreshold = 0;
static int _maxRefTuplesPerQry = 0;
static int _numTotalTuples = 0;
static int *_numActualRepeatsPerTuple = NULL;
static short _bgstClustSize[MAX_NUM_REFS];
static Hit *_list[MAX_NUM_REFS];


/**
 * Creates a lookup table.
 *
 * @param		refFile					Reference file.
 * @param		seedLen					Seed length.
 * @param 		maxHitsPerQry			Maximum hits per query.
 * @param		tupleIgnoreThreshold	Tuple that have number of repeats higher
 * than this value will be ignored.
 */
void lookupTable8Create(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold)
{
	_seedLen = seedLen;
	_maxHitsPerQry = maxHitsPerQry;
	_tupleIgnoreThreshold = tupleIgnoreThreshold;
	_numDistinctTuples = (int) pow((float) DNA_ALPHABET_SIZE, (int) _seedLen);
	_lookupTable = (int *) calloc(_numDistinctTuples, sizeof(int));
	_numActualRepeatsPerTuple = (int *) calloc(_numDistinctTuples, sizeof(int));
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
	int numIterations, i, j, offset = 0;
	int hashVal, lineLength, numBases = 0, refIdx = -1;
	char line[MAX_LINE_LENGTH];

	char tmpRefFile[MAX_FILE_NAME_LENGTH];
	sprintf(tmpRefFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, TEMP_REF_FILE);
	FILE *tmpRefFilePtr = fopen(tmpRefFile, "w");
	_numTotalTuples = 0;
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
					++_numActualRepeatsPerTuple[hashVal];
					++numRepeatsPerTuple[hashVal];
					++_numTotalTuples;
					fprintf(tmpRefFilePtr, "%d\t%d\t%d\n", hashVal, refIdx,
							numBases + i);
				}
				else if (_numRepeatsPerTuple[hashVal]
				                            == _numActualRepeatsPerTuple[hashVal] + 1)
				{
					_numActualRepeatsPerTuple[hashVal] = 0;
					numRepeatsPerTuple[hashVal] = 0;
					_numTotalTuples -= _tupleIgnoreThreshold;
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
	_refPos = (int *) calloc((_numTotalTuples + 1), sizeof(int));
	_refPos[0] = -1; /* First element set to -1 so that
	tuples that do not exist in the reference or that have number of repeats
	greater than a threshold value can point to this element. */
	_refPos2 = (int *) calloc((_numTotalTuples + 1), sizeof(int));
	_refPos2[0] = -1;
	_refIdx = (char *) calloc((_numTotalTuples + 1), sizeof(int));
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
				&& _numRepeatsPerTuple[i] <= _tupleIgnoreThreshold)
		{
			_lookupTable[i] = start;
			start += _numRepeatsPerTuple[i];
		}
	}
	time(&endTime);
	diffTime = difftime(endTime, startTime);
	fprintf(stderr, "done. (Time = %.2lf secs)\n", diffTime);
	fclose(tmpRefFilePtr);

	/* Set reference offsets. */
	time(&startTime);
	fprintf(stderr, "   Setting reference positions...");
	tmpRefFilePtr = fopen(tmpRefFile, "r");
	int index, refPos, tmp;
	while (fgets(line, MAX_LINE_LENGTH, tmpRefFilePtr) != NULL)
	{
		sscanf(line, "%d\t%d\t%d", &hashVal, &refIdx, &refPos);
		if (_numRepeatsPerTuple[hashVal] <= _tupleIgnoreThreshold)
		{
			index = _lookupTable[hashVal] + numRepeatsPerTuple[hashVal] - 1;
			_refIdx[index] = (char) refIdx;
			_refPos[index] = refPos;
			tmp = refIdx;
			tmp = tmp << REF_POS_BITS2;
			tmp += (refPos / _seedLen);
			_refPos2[index] = tmp;
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
}


/**
 * Reset key variables in this file.
 */
void lookupTable8Reset()
{

}


/**
 * Releases memory occupied by data structures in this file.
 */
void lookupTable8Delete()
{
	free(_lookupTable);
	_lookupTable = NULL;
	free(_refPos);
	_refPos = NULL;
	free(_refPos2);
	_refPos2 = NULL;
	free(_refIdx);
	_refIdx = NULL;
	_seedLen = 0;
	_maxHitsPerQry = 0;
	_numDistinctTuples = 0;
	free(_seed);
	_seed = NULL;
	free(_numRepeatsPerTuple);
	_numRepeatsPerTuple = NULL;
	free(_numActualRepeatsPerTuple);
	_numActualRepeatsPerTuple = NULL;
}


/**
 * Creates a new 'hit' instance.
 *
 * @param shift			Shift (reference position - query position) value.
 * @param refPos		Reference position.
 * @param clusterSize	Cluster size. Hits that have the same reference index
 * and shift belong to the same cluster.
 * @return Pointer to a new instance of 'hit'.
 */
static Hit *createHit(int shift, int refPos, short clusterSize)
{
	Hit *hit = (Hit *) malloc(sizeof(Hit));
	hit->shift = shift;
	hit->refPos = refPos;
	hit->clusterSize = clusterSize;
	hit->next = NULL;
	return hit;
}


/**
 * Searches the query sequence in the reference and returns the best-matching
 * reference coordinates.
 *
 * @param	query 				Query object.
 * @param 	seqLength 			Sequence length.
 * @param	refIdx_bestHits		Best-matching reference indexes.
 * @param	refPos_bestHits		Best-matching reference positions.
 * @return 	Number of hits in @a refIdx.
 */
int lookupTable8MapQry(char *qrySeq, int qryLen, char *refIdx_bestHits,
		int *refPos_bestHits)
{
	/* Initialize the list. */
	int i;
	for (i = 0; i < MAX_NUM_REFS; ++i)
	{
		_list[i] = NULL;
		_bgstClustSize[i] = 0;
	}

	/* Add new hits to proper clusters. */
	int numQryTuples = qryLen - _seedLen + 1;
	int j, hash, startIdx, endIdx;
	Hit *hit, *tmpHit, *prevHit;
	char refIdx;
	for (i = 0; i < numQryTuples; ++i)
	{
		hash = getHash(qrySeq + i, _seedLen);
		if (_numActualRepeatsPerTuple[hash] == 0)
			continue;

		startIdx = _lookupTable[hash];
		endIdx = startIdx + _numActualRepeatsPerTuple[hash] - 1;
		for (j = startIdx; j <= endIdx; ++j)
		{
			refIdx = _refIdx[j];
			hit = createHit( _refPos[j] - i, _refPos[j], 1);
			if (_list[refIdx] == NULL)
			{
				_list[refIdx] = hit;
				_bgstClustSize[refIdx] = 1;
			}
			else
			{
				tmpHit = _list[refIdx];
				prevHit = NULL;
				while (tmpHit != NULL && hit->shift > tmpHit->shift)
				{
					prevHit = tmpHit;
					tmpHit = tmpHit->next;
				}
				if (tmpHit != NULL)
				{
					if (tmpHit->shift == hit->shift)
					{
						tmpHit->clusterSize++;
						_bgstClustSize[refIdx] = max(
								_bgstClustSize[refIdx], tmpHit->clusterSize);
						free(hit);
					}
					else if (prevHit == NULL)
					{
						hit->next = tmpHit;
						_list[refIdx] = hit;
					}
					else
					{
						prevHit->next = hit;
						hit->next = tmpHit;
					}
				}
				else
					prevHit->next = hit;
			}
		}
	}

	/* Find the biggest cluster size. */
	short bgstClustSize = 0;
	for (i = 0; i < MAX_NUM_REFS; ++i)
		bgstClustSize = max(bgstClustSize, _bgstClustSize[i]);

	/* Find references that have biggest clusters. */
	int numBgstClust = 0;
	for (i = 0; i < MAX_NUM_REFS; ++i)
	{
		if (_bgstClustSize[i] == bgstClustSize)
		{
			tmpHit = _list[i];
			while (tmpHit != NULL)
			{
				if (tmpHit->clusterSize == bgstClustSize)
					++numBgstClust;
				tmpHit = tmpHit->next;
			}
		}
	}

	/* Fetch biggest clusters. */
	char *refIdxBestMtchs = (char *) calloc(numBgstClust, sizeof(char));
	int *refPosBestMtchs = (int *) calloc(numBgstClust, sizeof(int));
	j = 0;
	for (i = 0; i < MAX_NUM_REFS; ++i)
	{
		if (_bgstClustSize[i] == bgstClustSize)
		{
			tmpHit = _list[i];
			while (tmpHit != NULL)
			{
				if (tmpHit->clusterSize == bgstClustSize)
				{
					refIdxBestMtchs[j] = i;
					refPosBestMtchs[j] = tmpHit->refPos;
					++j;
				}
				tmpHit = tmpHit->next;
			}
		}
	}

	/* Assign biggest clusters to output array. */
	for (i = 0; i < _maxHitsPerQry; ++i)
		refIdx_bestHits[i] = -1;
	if (numBgstClust <= _maxHitsPerQry)
	{
		for (i = 0; i < numBgstClust; ++i)
		{
			refIdx_bestHits[i] = refIdxBestMtchs[i];
			refPos_bestHits[i] = refPosBestMtchs[i];
		}
	}
	else
	{
		int randNumArr[_maxHitsPerQry];
		arrGetRandomNums2(_maxHitsPerQry, 0, numBgstClust - 1, randNumArr);
		for (i = 0; i < _maxHitsPerQry; ++i)
		{
			refIdx_bestHits[i] = refIdxBestMtchs[randNumArr[i]];
			refPos_bestHits[i] = refPosBestMtchs[randNumArr[i]];
		}
		numBgstClust = _maxHitsPerQry;
	}

	/* Delete unused objects. */
	for (i = 0; i < MAX_NUM_REFS; ++i)
	{
		tmpHit = _list[i];
		while (tmpHit != NULL)
		{
			prevHit = tmpHit;
			tmpHit = tmpHit->next;
			free(prevHit);
		}
		_list[i] = NULL;
	}
	free(refIdxBestMtchs);
	free(refPosBestMtchs);

	return numBgstClust;
}
