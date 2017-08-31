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
 * Implements CPU-based filtering using seed and extend. Uses bigger seed
 * lengths for tuples that have higher-than-threshold number of repeats.
 */

#include "common.h"
#include "lookupTable11.h"
#include "preprocess.h"
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include "array.h"
#include <time.h>


#define	REF_POS_BITS2	25
#define	REF_POS_MASK2	33554431 /* Binary: 1111111111111111111111111 */
#define	NUM_BGST_CLUST	120
#define	INS_DEL_MARGIN	5
#define	ADJ_REF_TUPLE_DIST	12	/* Adjacent reference tuple distance. */
#define	MAX_NUM_LONG_HASHES	3000
#define	TUPLE_IGNORE_THRES2	75000

static int *_lookupTable = NULL;
static char *_refIdx = NULL;
static int *_refPos = NULL;
static char *_refIdx2 = NULL;
static int *_refPos2 = NULL;
static int _seedLen = 0;
static int _maxHitsPerQry = 0;
static int _numDistinctTuples = 0;
static int *_numRepeatsPerTuple = NULL;
static char *_seed = NULL;
static int _tupleIgnoreThreshold = 0;
static int _numTotalTuples = 0;
static int *_numActualRepeatsPerTuple = NULL;
static char *_refIdxArr = NULL;
static int *_refPosArr = NULL;
static int *_refPosHashMap = NULL;
static long long _maxRefPosComposite = 0;
static long *_longHashes = NULL;


/**
 * Creates a lookup table.
 *
 * @param		refFile					Reference file.
 * @param		seedLen					Seed length.
 * @param 		maxHitsPerQry			Maximum hits per query.
 * @param		tupleIgnoreThreshold	Tuple that have number of repeats higher
 * than this value will be ignored.
 */
void lookupTable12Create(const char *refFile, int seedLen, int maxHitsPerQry,
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
	_longHashes = (long *) calloc(MAX_NUM_LONG_HASHES, sizeof(long));

	/* First pass through reference to find the number of repeats for each
	 * distinct tuple. */
	time_t startTime, endTime;
	double diffTime;
	time(&startTime);
	fprintf(stderr, "   First pass...");
	FILE *filePtr = fopen(refFile, "r");
	int numIterations, i, j, offset = 0;
	int hashVal, prevHashVal, lineLength, numBases = 0, refIdx = -1;
	char line[MAX_LINE_LENGTH];

	char tmpRefFile[MAX_FILE_NAME_LENGTH];
	sprintf(tmpRefFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, TEMP_REF_FILE);
	FILE *tmpRefFilePtr = fopen(tmpRefFile, "w");
	setbuf(tmpRefFilePtr, NULL);
	_numTotalTuples = 0;
	long long refPosComposite = 0;
	int numLongHashPos = 0;
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
			prevHashVal = -1;
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
					refPosComposite = refIdx;
					refPosComposite = refPosComposite << REF_POS_BITS2;
					refPosComposite += ((numBases + i) / _seedLen);
					_maxRefPosComposite = max(_maxRefPosComposite, refPosComposite);
				}
				else if (_numRepeatsPerTuple[hashVal]
				                            == _numActualRepeatsPerTuple[hashVal] + 1)
				{
					_numActualRepeatsPerTuple[hashVal] = 0;
					numRepeatsPerTuple[hashVal] = 0;
					_numTotalTuples -= _tupleIgnoreThreshold;
					fprintf(tmpRefFilePtr, "%d\t%d\t%d\n", hashVal, refIdx,
							numBases + i);
					refPosComposite = refIdx;
					refPosComposite = refPosComposite << REF_POS_BITS2;
					refPosComposite += ((numBases + i) / _seedLen);
					_maxRefPosComposite = max(_maxRefPosComposite, refPosComposite);
					++numLongHashPos;
				}
				else if (_numRepeatsPerTuple[hashVal] <= TUPLE_IGNORE_THRES2)
				{
					fprintf(tmpRefFilePtr, "%d\t%d\t%d\n", hashVal, refIdx,
							numBases + i);
					refPosComposite = refIdx;
					refPosComposite = refPosComposite << REF_POS_BITS2;
					refPosComposite += ((numBases + i) / _seedLen);
					_maxRefPosComposite = max(_maxRefPosComposite, refPosComposite);
					++numLongHashPos;
				}
				prevHashVal = hashVal;
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
	_refIdx = (char *) calloc((_numTotalTuples + 1), sizeof(int));
	_refIdx[0] = -1;
	_refIdx2 = (char *) calloc(numLongHashPos, sizeof(char));
	_refPos2 = (int *) calloc(numLongHashPos, sizeof(int));
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
	_refPosHashMap = (int *) calloc(_maxRefPosComposite, sizeof(int));
	tmpRefFilePtr = fopen(tmpRefFile, "r");
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

			refPosComposite = refIdx;
			refPosComposite = refPosComposite << REF_POS_BITS2;
			refPosComposite += (refPos / _seedLen);
			_refPosHashMap[refPosComposite] = hashVal;
		}
		else if (_numRepeatsPerTuple[hashVal] <= TUPLE_IGNORE_THRES2)
		{
			refPosComposite = refIdx;
			refPosComposite = refPosComposite << REF_POS_BITS2;
			refPosComposite += (refPos / _seedLen);
			_refPosHashMap[refPosComposite] = hashVal;
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

	_refIdxArr = (char *) malloc(NUM_BGST_CLUST * sizeof(char));
	_refPosArr = (int *) malloc(NUM_BGST_CLUST * sizeof(int));

	/* Print hash table to a file. */
//	FILE *keysFilePtr = fopen("/home/pgupta/data/gpusw/keys5000.txt", "w");
//	for (i = 0; i < _numDistinctTuples; ++i)
//		fprintf(keysFilePtr, "%d\t%d\n", _lookupTable[i], _numActualRepeatsPerTuple[i]);
//	fclose(keysFilePtr);
//
//	FILE *valsFilePtr = fopen("/home/pgupta/data/gpusw/vals5000.txt", "w");
//	for (i = 0; i <= _numTotalTuples; ++i)
//		fprintf(valsFilePtr, "%d\t%d\n", _refIdx[i], _refPos[i]);
//	fclose(valsFilePtr);
//
//	FILE *keysFilePtr2 = fopen("/home/pgupta/data/gpusw/keys5000_2.txt", "w");
//	for (i = 0; i < _maxRefPosComposite; ++i)
//		fprintf(keysFilePtr2, "%d\n", _refPosHashMap[i]);
//	fclose(keysFilePtr2);
}


/**
 * Releases memory occupied by data structures in this file.
 */
void lookupTable12Delete()
{
	free(_lookupTable);
	_lookupTable = NULL;
	free(_refPos);
	_refPos = NULL;
	free(_refIdx);
	_refIdx = NULL;
	free(_refIdx2);
	_refIdx2 = NULL;
	free(_refPos2);
	_refPos2 = NULL;
	_seedLen = 0;
	_maxHitsPerQry = 0;
	_numDistinctTuples = 0;
	free(_seed);
	_seed = NULL;
	free(_numRepeatsPerTuple);
	_numRepeatsPerTuple = NULL;
	free(_numActualRepeatsPerTuple);
	_numActualRepeatsPerTuple = NULL;
	free(_refIdxArr);
	_refIdxArr = NULL;
	free(_refPosArr);
	_refPosArr = NULL;
	free(_refPosHashMap);
	_refPosHashMap = NULL;
//	remove(_tmpRefFile2);
}
