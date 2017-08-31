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
 * Implements CPU-based filtering; uses arrays.
 */

#include "common.h"
#include "lookupTable9.h"
#include "preprocess.h"
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include "array.h"


#define	REF_POS_BITS2	25
#define	REF_POS_MASK2	33554431 /* Binary: 1111111111111111111111111 */
#define	TOTAL_NUM_CHRS	25
#define	MAX_NUM_SHIFTS	249255000

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
static int **_shift1 = NULL;
static short **_shift2 = NULL;
static short *_maxClustSizes = NULL;


/**
 * Creates a lookup table.
 *
 * @param		refFile					Reference file.
 * @param		seedLen					Seed length.
 * @param 		maxHitsPerQry			Maximum hits per query.
 * @param		tupleIgnoreThreshold	Tuple that have number of repeats higher
 * than this value will be ignored.
 */
void lookupTable9Create(const char *refFile, int seedLen, int maxHitsPerQry,
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

	_shift1 = (int **) calloc(MAX_NUM_SHIFTS, sizeof(int *));
	_shift2 = (short **) calloc(MAX_NUM_SHIFTS, sizeof(short *));
	_maxClustSizes = (short *) calloc(MAX_NUM_SHIFTS, sizeof(short));
	for (i = 0; i < MAX_NUM_SHIFTS; ++i)
	{
		_shift1[i] = (int *) calloc(MAX_NUM_REFS, sizeof(int));
		_shift2[i] = (short *) calloc(MAX_NUM_REFS, sizeof(short));
	}

	fprintf(stderr, "   Freeing resources...");
	fclose(tmpRefFilePtr);
	remove(tmpRefFile);
	fclose(filePtr);
	free(numRepeatsPerTuple);
	fprintf(stderr, "done.\n");
}


/**
 * This is a temporary way to create the lookup table.
 *
 * @param keysFile	File containing keys.
 * @param valsFile	File containing values and number of repeats per tuple.
 * @param maxHits	Max number of hits per second.
 * @param seedLen	Seed length.
 * @param tupleIgnoreThres	Threshold value used to ignore reference tuples.
 */
void lookupTable9Create(const char *keysFile, const char *valsFile, int maxHits,
		int seedLen, int tupleIgnoreThres)
{
	FILE *keysFilePtr = fopen(keysFile, "r");
	char line[200];
	int numKeys = 0;
	while (fgets(line, 200, keysFilePtr) != NULL)
		++numKeys;
	rewind(keysFilePtr);

	_lookupTable = (int *) calloc((numKeys + 1), sizeof(int));
	_numActualRepeatsPerTuple = (int *) calloc(numKeys, sizeof(int));
	int i = 0;
	while (fgets(line, 200, keysFilePtr) != NULL)
	{
		sscanf(line, "%d\t%d", &_lookupTable[i], &_numActualRepeatsPerTuple[i]);
		++i;
	}
	fclose(keysFilePtr);

	FILE *valsFilePtr = fopen(valsFile, "r");
	int numVals = 0;
	while (fgets(line, 200, valsFilePtr) != NULL)
		++numVals;
	rewind(valsFilePtr);

	_refIdx = (char *) calloc((numVals + 1), sizeof(char));
	_refPos = (int *) calloc((numVals + 1), sizeof(int));
	_refPos2 = (int *) calloc((numVals + 1), sizeof(int));
	i = 0;
	while (fgets(line, 200, valsFilePtr) != NULL)
	{
		sscanf(line, "%d", &_refPos2[i]);
		_refIdx[i] = (char) (_refPos2[i] >> REF_POS_BITS2);
		_refPos[i] = (_refPos2[i] & REF_POS_MASK2) * seedLen;
		++i;
	}
	fclose(valsFilePtr);

	_maxHitsPerQry = maxHits;
	int numQryTuples = MAX_QRY_SEQ_LENGTH - seedLen;
	_maxRefTuplesPerQry = numQryTuples * tupleIgnoreThres;
	_seedLen = seedLen;
	_tupleIgnoreThreshold = tupleIgnoreThres;
	_numDistinctTuples = (int) pow((float) DNA_ALPHABET_SIZE, (int) _seedLen);
	_numTotalTuples = numVals;

	_shift1 = (int **) calloc(MAX_NUM_SHIFTS, sizeof(int *));
	_shift2 = (short **) calloc(MAX_NUM_SHIFTS, sizeof(short *));
	_maxClustSizes = (short *) calloc(MAX_NUM_SHIFTS, sizeof(short));
	for (i = 0; i < MAX_NUM_SHIFTS; ++i)
	{
		_shift1[i] = (int *) calloc(MAX_NUM_REFS, sizeof(int));
		_shift2[i] = (short *) calloc(MAX_NUM_REFS, sizeof(short));
	}
}


/**
 * Reset key variables in this file.
 */
void lookupTable9Reset()
{

}


/**
 * Releases memory occupied by data structures in this file.
 */
void lookupTable9Delete()
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
	int i;
	for (i = 0; i < MAX_NUM_SHIFTS; ++i)
	{
		free(_shift1[i]);
		free(_shift2[i]);
	}
	free(_shift1);
	free(_shift2);
	free(_maxClustSizes);
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
int lookupTable9MapQry(char *qrySeq, int qryLen, char *refIdx_bestHits,
		int *refPos_bestHits)
{
	char refIdx;
	short bgstClustSize = 0;
	int i, j, k, hash, startIdx, endIdx, shift, *refPosArr;
	int numQryTuples = qryLen - _seedLen + 1;
	for (i = 0; i < numQryTuples; ++i)
	{
		hash = getHash(qrySeq + i, _seedLen);
		if (_numActualRepeatsPerTuple[hash] == 0)
			continue;

		startIdx = _lookupTable[hash];
		endIdx = startIdx + _numActualRepeatsPerTuple[hash] - 1;
		for (j = startIdx; j <= endIdx; ++j)
		{
			shift = _refPos[j] - i + MAX_QRY_SEQ_LENGTH;
//			if (_shift1[shift] == NULL)
//			{
//				refPosArr = (int *) calloc(MAX_NUM_REFS, sizeof(int));
//				for (k = 0; k < MAX_NUM_REFS; ++k)
//					refPosArr[k] = INT_MAX;
//				_shift1[shift] = refPosArr;
////				_shift2[shift] = (short *) calloc(MAX_NUM_REFS, sizeof(short));
//			}
			refIdx = _refIdx[j];
			_shift1[shift][refIdx] = min(_shift1[shift][refIdx], _refPos[j]);
			++(_shift2[shift][refIdx]);
			_maxClustSizes[shift] = max(_maxClustSizes[shift],
					_shift2[shift][refIdx]);
			bgstClustSize = max(bgstClustSize, _maxClustSizes[shift]);
		}
	}

	/* Find number of biggest clusters. */
	int numBgstClust = 0;
//	for (i = 0; i < MAX_NUM_SHIFTS; ++i)
//	{
//		if (maxClustSizes[i] == bgstClustSize)
//		{
//			for (j = 0; j < MAX_NUM_REFS; ++j)
//			{
//				if (shift2[i][j] == bgstClustSize)
//					++numBgstClust;
//			}
//		}
//	}

//	/* Fetch biggest clusters. */
//	char refIdxBestMtchs[numBgstClust];
//	int refPosBestMtchs[numBgstClust];
//	k = 0;
//	for (i = 0; i < MAX_NUM_SHIFTS; ++i)
//	{
//		if (maxClustSizes[i] == bgstClustSize)
//		{
//			for (j = 0; j < MAX_NUM_REFS; ++j)
//			{
//				if (shift2[i][j] == bgstClustSize)
//				{
//					refIdxBestMtchs[k] = j;
//					refPosBestMtchs[k] = shift1[i][j];
//					++k;
//				}
//			}
//		}
//	}
//
//	/* Assign biggest clusters to output array. */
//	for (i = 0; i < _maxHitsPerQry; ++i)
//		refIdx_bestHits[i] = -1;
//	if (numBgstClust <= _maxHitsPerQry)
//	{
//		for (i = 0; i < numBgstClust; ++i)
//		{
//			refIdx_bestHits[i] = refIdxBestMtchs[i];
//			refPos_bestHits[i] = refPosBestMtchs[i];
//		}
//	}
//	else
//	{
//		int randNumArr[_maxHitsPerQry];
//		arrGetRandomNums2(_maxHitsPerQry, 0, numBgstClust - 1, randNumArr);
//		for (i = 0; i < _maxHitsPerQry; ++i)
//		{
//			refIdx_bestHits[i] = refIdxBestMtchs[randNumArr[i]];
//			refPos_bestHits[i] = refPosBestMtchs[randNumArr[i]];
//		}
//		numBgstClust = _maxHitsPerQry;
//	}

	/* Reset data structures. */
	for (i = 0; i < MAX_NUM_SHIFTS; ++i)
	{
		if (_maxClustSizes[i] > 0)
		{
			for (j = 0; j < MAX_NUM_REFS; ++j)
				_shift1[i][j] = 0;
		}
		_maxClustSizes[i] = 0;
	}

	return numBgstClust;
}
