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
 * Implements CPU-based filtering using seed and extend.
 */

#include "common.h"
#include "lookupTable11.h"
#include "preprocess.h"
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include "array.h"


#define	REF_POS_BITS2	25
#define	REF_POS_MASK2	33554431 /* Binary: 1111111111111111111111111 */
#define	NUM_BGST_CLUST	120
#define	INS_DEL_MARGIN	5
#define	ADJ_REF_TUPLE_DIST	12	/* Adjacent reference tuple distance. */
#define	TUPLE_IGNORE_THRES2	75000

static int *_lookupTable = NULL;
static char *_refIdx = NULL;
static int *_refPos = NULL;
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
static char _overflowFile[MAX_FILE_NAME_LENGTH];
static FILE *_overflowFilePtr = NULL;
static char _line[MAX_LINE_LENGTH];
static int *_prevHashes = NULL;
static int _numPrevHashes = 0;


/**
 * Creates a lookup table.
 *
 * @param		refFile					Reference file.
 * @param		seedLen					Seed length.
 * @param 		maxHitsPerQry			Maximum hits per query.
 * @param		tupleIgnoreThreshold	Tuple that have number of repeats higher
 * than this value will be ignored.
 */
void lookupTable11Create(const char *refFile, int seedLen, int maxHitsPerQry,
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
	setbuf(tmpRefFilePtr, NULL);
	_numTotalTuples = 0;
	long long refPosComposite = 0;
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
//			for (i = 0; i < numIterations; i += ADJ_REF_TUPLE_DIST)
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

	_prevHashes = (int *) calloc(MAX_QRY_SEQ_LENGTH, sizeof(int));

	/* Print hash table to a file. */
//	FILE *keysFilePtr = fopen("/home/pgupta/data/gpusw/keys500_12M11.txt", "w");
//	for (i = 0; i < _numDistinctTuples; ++i)
//		fprintf(keysFilePtr, "%d\t%d\n", _lookupTable[i], _numActualRepeatsPerTuple[i]);
//	fclose(keysFilePtr);
//
//	FILE *valsFilePtr = fopen("/home/pgupta/data/gpusw/vals500_12M11.txt", "w");
//	for (i = 0; i <= _numTotalTuples; ++i)
//		fprintf(valsFilePtr, "%d\t%d\n", _refIdx[i], _refPos[i]);
//	fclose(valsFilePtr);
}


/**
 * This is a temporary way to create the lookup table.
 *
 * @param keysFile	File containing keys.
 * @param valsFile	File containing values and number of repeats per tuple.
 * @param keysFile2	File containing the mapping of reference positions to hashes.
 * @param maxHits	Max number of hits per second.
 * @param seedLen	Seed length.
 * @param tupleIgnoreThres	Threshold value used to ignore reference tuples.
 */
void lookupTable11Create(const char *keysFile, const char *valsFile,
		const char *keysFile2, int maxHits, int seedLen, int tupleIgnoreThres)
{
	time_t startTime, endTime;
	double diffTime = 0.0;
	time(&startTime);
	FILE *keysFilePtr = fopen(keysFile, "r");
	char line[200];
	int numKeys = 0;
	while (fgets(line, 200, keysFilePtr) != NULL)
		++numKeys;
	rewind(keysFilePtr);

	_lookupTable = (int *) calloc((numKeys + 1), sizeof(int));
	_numRepeatsPerTuple = (int *) calloc(numKeys, sizeof(int));
	int i = 0;
	while (fgets(line, 200, keysFilePtr) != NULL)
	{
		sscanf(line, "%d\t%d", &_lookupTable[i], &_numRepeatsPerTuple[i]);
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
	i = 0;
	long long refPosComposite = 0;
	while (fgets(line, 200, valsFilePtr) != NULL)
	{
		sscanf(line, "%d\t%d", &_refIdx[i], &_refPos[i]);
		refPosComposite = _refIdx[i];
		refPosComposite = refPosComposite << REF_POS_BITS2;
		refPosComposite += (_refPos[i] / seedLen);
		_maxRefPosComposite = max(_maxRefPosComposite, refPosComposite);
		++i;
	}
	fclose(valsFilePtr);

	_refPosHashMap = (int *) calloc(_maxRefPosComposite, sizeof(int));
	FILE *keysFilePtr2 = fopen(keysFile2, "r");
	i = 0;
	while (fgets(line, 200, keysFilePtr2) != NULL)
	{
		sscanf(line, "%d", &_refPosHashMap[i]);
		++i;
	}
	fclose(keysFilePtr2);

	_maxHitsPerQry = maxHits;
	_seedLen = seedLen;
	_tupleIgnoreThreshold = tupleIgnoreThres;
	_numDistinctTuples = (int) pow((float) DNA_ALPHABET_SIZE, (int) _seedLen);
	_numTotalTuples = numVals;
	_refIdxArr = (char *) malloc(NUM_BGST_CLUST * sizeof(char));
	_refPosArr = (int *) malloc(NUM_BGST_CLUST * sizeof(int));

	_prevHashes = (int *) calloc(MAX_QRY_SEQ_LENGTH, sizeof(int));

	sprintf(_overflowFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "overflow.txt");
	_overflowFilePtr = fopen(_overflowFile, "r");

	time(&endTime);
	diffTime = difftime(endTime, startTime);
	fprintf(stderr, "(Time = %.2lf)...", diffTime);

//	int xx = 3357627;
//	fprintf(stderr, "repeats[%d] = %d\n", xx, _numRepeatsPerTuple[xx]);
//	xx = 1118515;
//	fprintf(stderr, "repeats[%d] = %d\n", xx, _numRepeatsPerTuple[xx]);
//	xx = 3346705;
//	fprintf(stderr, "repeats[%d] = %d\n", xx, _numRepeatsPerTuple[xx]);
//	xx = 4453171;
//	fprintf(stderr, "repeats[%d] = %d\n", xx, _numRepeatsPerTuple[xx]);
//	xx = 384224;
//	fprintf(stderr, "repeats[%d] = %d\n", xx, _numRepeatsPerTuple[xx]);
//	xx = 8717395;
//	fprintf(stderr, "repeats[%d] = %d\n", xx, _numRepeatsPerTuple[xx]);
//	xx = 15728443;
//	fprintf(stderr, "repeats[%d] = %d\n", xx, _numRepeatsPerTuple[xx]);
//	xx = 11677695;
//	fprintf(stderr, "repeats[%d] = %d\n", xx, _numRepeatsPerTuple[xx]);
}


/**
 * Creates a lookup table.
 *
 * @param		refFile					Reference file.
 * @param		seedLen					Seed length.
 * @param 		maxHitsPerQry			Maximum hits per query.
 * @param		tupleIgnoreThreshold	Tuple that have number of repeats higher
 * than this value will be ignored.
 */
void lookupTable11Create2(const char *refFile, int seedLen, int maxHitsPerQry,
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
	int hashVal, prevHashVal, lineLength, numBases = 0, refIdx = -1;
	char line[MAX_LINE_LENGTH];

	char tmpRefFile[MAX_FILE_NAME_LENGTH];
	sprintf(tmpRefFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, TEMP_REF_FILE);
	FILE *tmpRefFilePtr = fopen(tmpRefFile, "w");
	setbuf(tmpRefFilePtr, NULL);
	_numTotalTuples = 0;
	long long refPosComposite = 0;
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
				++numRepeatsPerTuple[hashVal];
				if (_numRepeatsPerTuple[hashVal] <= _tupleIgnoreThreshold)
				{
					++_numActualRepeatsPerTuple[hashVal];
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
					_numTotalTuples -= _tupleIgnoreThreshold;
					fprintf(tmpRefFilePtr, "%d\t%d\t%d\n", hashVal, refIdx,
							numBases + i);
					refPosComposite = refIdx;
					refPosComposite = refPosComposite << REF_POS_BITS2;
					refPosComposite += ((numBases + i) / _seedLen);
					_maxRefPosComposite = max(_maxRefPosComposite, refPosComposite);
				}
				else if (_numRepeatsPerTuple[hashVal] <= TUPLE_IGNORE_THRES2)
				{
					fprintf(tmpRefFilePtr, "%d\t%d\t%d\n", hashVal, refIdx,
							numBases + i);
					refPosComposite = refIdx;
					refPosComposite = refPosComposite << REF_POS_BITS2;
					refPosComposite += ((numBases + i) / _seedLen);
					_maxRefPosComposite = max(_maxRefPosComposite, refPosComposite);
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
	_refIdx = (char *) calloc((_numTotalTuples + 1), sizeof(int));
	_refIdx[0] = -1;
	time(&endTime);
	diffTime = difftime(endTime, startTime);
	fprintf(stderr, "done. (Time = %.2lf secs)\n", diffTime);

	/* Set values in the lookup table. */
	time(&startTime);
	fprintf(stderr, "   Set values in lookup table...");
	int start = 1, start2 = 0;
	for (i = 0; i < _numDistinctTuples; ++i)
	{
		if (_numRepeatsPerTuple[i] > 0
				&& _numRepeatsPerTuple[i] <= _tupleIgnoreThreshold)
		{
			_lookupTable[i] = start;
			start += _numRepeatsPerTuple[i];
		}
		else if (_numRepeatsPerTuple[i] <= TUPLE_IGNORE_THRES2)
		{
			_lookupTable[i] = start2;
			start2 += _numRepeatsPerTuple[i];
		}
	}
	time(&endTime);
	diffTime = difftime(endTime, startTime);
	fprintf(stderr, "done. (Time = %.2lf secs)\n", diffTime);
	fclose(tmpRefFilePtr);

	/* Set reference offsets. */
	time(&startTime);
	fprintf(stderr, "   Setting reference positions...");
	int numOverflowRepeats = 0;
	for (i = 0; i < _numDistinctTuples; ++i)
	{
		if (_numRepeatsPerTuple[i] > _tupleIgnoreThreshold
				&& _numRepeatsPerTuple[i] <= TUPLE_IGNORE_THRES2)
			numOverflowRepeats += _numRepeatsPerTuple[i];
	}
	char *refIdx2 = (char *) calloc(numOverflowRepeats + 1, sizeof(char));
	refIdx2[0] = -1;
	int *refPos2 = (int *) calloc(numOverflowRepeats + 1, sizeof(int));
	refPos2[0] = -1;
	_refPosHashMap = (int *) calloc(_maxRefPosComposite, sizeof(int));
	tmpRefFilePtr = fopen(tmpRefFile, "r");
	int index, refPos;
	char unhash[_seedLen + 1];
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
			index = _lookupTable[hashVal] + numRepeatsPerTuple[hashVal] - 1;
			refIdx2[index] = (char) refIdx;
			refPos2[index] = refPos;
			--numRepeatsPerTuple[hashVal];

			refPosComposite = refIdx;
			refPosComposite = refPosComposite << REF_POS_BITS2;
			refPosComposite += (refPos / _seedLen);
			_refPosHashMap[refPosComposite] = hashVal;
		}
	}
	time(&endTime);
	diffTime = difftime(endTime, startTime);
	fprintf(stderr, "done. (Time = %.2lf secs)\n", diffTime);

	/* Overflow reference coordinates. */
	sprintf(_overflowFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, "overflow.txt");
	_overflowFilePtr = fopen(_overflowFile, "w");
	int startIdx, endIdx;
	for (i = 0; i < _numDistinctTuples; ++i)
	{
		if (_numRepeatsPerTuple[i] > _tupleIgnoreThreshold
				&& _numRepeatsPerTuple[i] <= TUPLE_IGNORE_THRES2)
		{
			startIdx = _lookupTable[i];
			endIdx = startIdx + _numRepeatsPerTuple[i] - 1;
			_lookupTable[i] = ftell(_overflowFilePtr);
			for (j = endIdx; j >= startIdx; --j)
				fprintf(_overflowFilePtr, "%d\t%d\n", refIdx2[j], refPos2[j]);
		}
	}
	fclose(_overflowFilePtr);
	free(refIdx2);
	free(refPos2);

	/* Free resources. */
	fprintf(stderr, "   Freeing resources...");
	fclose(tmpRefFilePtr);
	remove(tmpRefFile);
	fclose(filePtr);
	free(numRepeatsPerTuple);
	fprintf(stderr, "done.\n");

	_refIdxArr = (char *) malloc(NUM_BGST_CLUST * sizeof(char));
	_refPosArr = (int *) malloc(NUM_BGST_CLUST * sizeof(int));

	_prevHashes = (int *) calloc(MAX_QRY_SEQ_LENGTH, sizeof(int));

	_overflowFilePtr = fopen(_overflowFile, "r");

	/* Print hash table to a file. */
	FILE *keysFilePtr = fopen("/home/pgupta/data/gpusw/keys1000.txt", "w");
	for (i = 0; i < _numDistinctTuples; ++i)
		fprintf(keysFilePtr, "%d\t%d\n", _lookupTable[i], _numRepeatsPerTuple[i]);
	fclose(keysFilePtr);

	FILE *valsFilePtr = fopen("/home/pgupta/data/gpusw/vals1000.txt", "w");
	for (i = 0; i <= _numTotalTuples; ++i)
		fprintf(valsFilePtr, "%d\t%d\n", _refIdx[i], _refPos[i]);
	fclose(valsFilePtr);

	FILE *keysFilePtr2 = fopen("/home/pgupta/data/gpusw/keys1000_2.txt", "w");
	for (i = 0; i < _maxRefPosComposite; ++i)
		fprintf(keysFilePtr2, "%d\n", _refPosHashMap[i]);
	fclose(keysFilePtr2);
}


/**
 * Reset key variables in this file.
 */
void lookupTable11Reset()
{

}


/**
 * Releases memory occupied by data structures in this file.
 */
void lookupTable11Delete()
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
	free(_numActualRepeatsPerTuple);
	_numActualRepeatsPerTuple = NULL;
	free(_refIdxArr);
	_refIdxArr = NULL;
	free(_refPosArr);
	_refPosArr = NULL;
	free(_refPosHashMap);
	_refPosHashMap = NULL;
	fclose(_overflowFilePtr);
//	remove(_overflowFile);
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
int lookupTable11MapQry(char *qrySeq, int qryLen, char *refIdx_bestHits,
		int *refPos_bestHits)
{
	int i, j, hash, startIdx, endIdx, bgstClustSize = 0, clustSize;
	int isContinous, isContinous2, fileOffset, tmp, refIdx, refPos;
	int numBgstClust = 0, numBgstClust2 = 0;
	int numQryTuples = qryLen - _seedLen + 1;
	numQryTuples = min(25, numQryTuples);
	for (i = 0; i < numQryTuples; ++i)
	{
		hash = getHash(qrySeq + i, _seedLen);
		if (_numRepeatsPerTuple[hash] > _tupleIgnoreThreshold)
			continue;
		startIdx = _lookupTable[hash];
		endIdx = startIdx + _numRepeatsPerTuple[hash] - 1;
		for (j = endIdx; j >= startIdx; --j)
		{
			clustSize = getClustSize(qrySeq + i, qryLen, _refIdx[j],
					_refPos[j], &isContinous);
			if (bgstClustSize <= clustSize)
			{
				if (bgstClustSize < clustSize)
				{
					bgstClustSize = clustSize;
					_refIdxArr[0] = _refIdx[j];
					_refPosArr[0] = _refPos[j];
					numBgstClust = 1;
					numBgstClust2 = 1;
					isContinous2 = -1;
				}
				else if (numBgstClust < NUM_BGST_CLUST)
				{
					if (isContinous2 < isContinous)
					{
						bgstClustSize = clustSize;
						_refIdxArr[0] = _refIdx[j];
						_refPosArr[0] = _refPos[j];
						numBgstClust = 1;
						numBgstClust2 = 1;
						isContinous2 = isContinous;
					}
					else if (isContinous2 == isContinous)
					{
						_refIdxArr[numBgstClust] = _refIdx[j];
						_refPosArr[numBgstClust] = _refPos[j];
						++numBgstClust;
						++numBgstClust2;
					}
				}
				else
					++numBgstClust;
			}
		}
	}

	if (numBgstClust2 <= _maxHitsPerQry)
	{
		for (i = 0; i < _maxHitsPerQry; ++i)
		{
			if (i < numBgstClust2)
			{
				refIdx_bestHits[i] = _refIdxArr[i];
				refPos_bestHits[i] = _refPosArr[i];
			}
			else
				refIdx_bestHits[i] = -1;
		}
		numBgstClust = numBgstClust2;
	}
	else if (numBgstClust2 >= numBgstClust)
	{
		int *randNumArr = (int *) calloc(_maxHitsPerQry, sizeof(int));
		arrGetRandomNums2(_maxHitsPerQry, 0, numBgstClust2 - 1, randNumArr);
		for (i = 0; i < _maxHitsPerQry; ++i)
		{
			refIdx_bestHits[i] = _refIdxArr[randNumArr[i]];
			refPos_bestHits[i] = _refPosArr[randNumArr[i]];
		}
		free(randNumArr);
		numBgstClust = _maxHitsPerQry;
	}
	else
		numBgstClust = 0;

	return numBgstClust;
}


/**
 * Finds and returns the number of reference hits for tuples in a query.
 *
 * @param qrySeq	Query sequence
 * @param qryLen	Query sequence length.
 * @param refIdx	Reference index.
 * @param refPos	Reference position.
 * @param[out]	isContinous	A value of 0 indicates there is a mismatching
 * tuple; 1 indicates there are no mismatching tuples.
 * @return Number of reference hits for tuples in a query.
 */
static int getClustSize(char *qrySeq, int qryLen, char refIdx, int refPos,
		int *isContinous)
{
	int clustSize = 0, hash, offset = 0;
	int numQryTuples =  qryLen - _seedLen;
	long refPosComposite = refIdx;
	refPosComposite = refPosComposite << REF_POS_BITS2;
	refPosComposite += (refPos / _seedLen);
	int refPosTmp = 0, refIdxTmp;
	int i, haveEarlierTuplesMatched = 0;
	*isContinous = 1;
	while (offset < numQryTuples)
	{
		hash = getHash(qrySeq + offset, _seedLen);
		refPosTmp = (refPosComposite & 33554431) * _seedLen;
		refIdxTmp = refPosComposite >> REF_POS_BITS2;
		if (_maxRefPosComposite > refPosComposite
				&& _refPosHashMap[refPosComposite] == hash)
		{
			++clustSize;
			offset += _seedLen;
			haveEarlierTuplesMatched = 1;
		}
		else if (haveEarlierTuplesMatched == 0)
			++offset;
		else
		{
			offset += _seedLen;
			*isContinous = 0;
		}
		++refPosComposite;
	}

	return clustSize;
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
int lookupTable11MapQry2(char *qrySeq, int qryLen, char *refIdx_bestHits,
		int *refPos_bestHits)
{
	int i, j, hash, startIdx, endIdx, bgstClustSize = 0, clustSize;
	int isContinous, isContinous2, fileOffset, tmp, refIdx, refPos;
	int numBgstClust = 0, numBgstClust2 = 0;
	int numQryTuples = qryLen - _seedLen + 1;
	numQryTuples = min(25, numQryTuples);
	_numPrevHashes = 0;

	for (i = 0; i < numQryTuples; ++i)
	{
		hash = getHash(qrySeq + i, _seedLen);
		if (_numRepeatsPerTuple[hash] > _tupleIgnoreThreshold)
		{
			_prevHashes[_numPrevHashes] = hash;
			++_numPrevHashes;
			continue;
		}
		startIdx = _lookupTable[hash];
		endIdx = startIdx + _numRepeatsPerTuple[hash] - 1;
		for (j = endIdx; j >= startIdx; --j)
		{
//			if (_refIdx[j] == 0 && _refPos[j] >= 1705000 && _refPos[j] <= 1706000)
//			{
//				fprintf(stderr, "refIdx = %d, refPos = %d\n", _refIdx[j], _refPos[j]);
				clustSize = getClustSize2(qrySeq + i, qryLen, _refIdx[j],
						_refPos[j], &isContinous);
				if (bgstClustSize <= clustSize)
				{
					if (bgstClustSize < clustSize)
					{
						bgstClustSize = clustSize;
						_refIdxArr[0] = _refIdx[j];
						_refPosArr[0] = _refPos[j];
						numBgstClust = 1;
						numBgstClust2 = 1;
						isContinous2 = isContinous;
					}
					else if (numBgstClust < NUM_BGST_CLUST)
					{
						if (isContinous2 < isContinous)
						{
							bgstClustSize = clustSize;
							_refIdxArr[0] = _refIdx[j];
							_refPosArr[0] = _refPos[j];
							numBgstClust = 1;
							numBgstClust2 = 1;
							isContinous2 = isContinous;
						}
						else if (isContinous2 == isContinous)
						{
							_refIdxArr[numBgstClust] = _refIdx[j];
							_refPosArr[numBgstClust] = _refPos[j];
							++numBgstClust;
							++numBgstClust2;
						}
					}
					else
						++numBgstClust;
				}
//			}
		}
	}

//	fprintf(stderr, "*bgstClustSize = %d, numBgstClust = %d\n",
//			bgstClustSize, numBgstClust);

	i = 0;
	while (bgstClustSize < 7)
	{
		int fileOffset = 0, refIdx, refPos;
		while (TRUE)
		{
			hash = getHash(qrySeq + i, _seedLen);
			if (_numRepeatsPerTuple[hash] > _tupleIgnoreThreshold
				&& _numRepeatsPerTuple[hash] <= TUPLE_IGNORE_THRES2)
			{
				fileOffset = _lookupTable[hash];
				break;
			}
			++i;
		}
		if (i >= 20)
			break;
		int k = _numRepeatsPerTuple[hash];

		if (fseek(_overflowFilePtr, fileOffset, SEEK_SET) != 0)
		{
			fprintf(stderr, "Error: fseek returned non-zero value.\n");
			exit(0);
		}
		while (fgets(_line, MAX_LINE_LENGTH, _overflowFilePtr))
		{
			sscanf(_line, "%d\t%d", &refIdx, &refPos);
//			if (refIdx == 0 && refPos >= 1705000 && refPos <= 1706000)
//			{
//				fprintf(stderr, "refIdx = %d, refPos = %d\n", refIdx, refPos);
				clustSize = getClustSize2(qrySeq + i, qryLen, refIdx, refPos,
						&isContinous);
				if (bgstClustSize <= clustSize)
				{
					if (bgstClustSize < clustSize)
					{
						bgstClustSize = clustSize;
						_refIdxArr[0] = refIdx;
						_refPosArr[0] = refPos;
						numBgstClust = 1;
						numBgstClust2 = 1;
						isContinous2 = isContinous;
					}
					else if (numBgstClust < NUM_BGST_CLUST)
					{
						if (isContinous2 < isContinous)
						{
							bgstClustSize = clustSize;
							_refIdxArr[0] = refIdx;
							_refPosArr[0] = refPos;
							numBgstClust = 1;
							numBgstClust2 = 1;
							isContinous2 = isContinous;
						}
						else if (isContinous2 == isContinous)
						{
							_refIdxArr[numBgstClust] = refIdx;
							_refPosArr[numBgstClust] = refPos;
							++numBgstClust;
							++numBgstClust2;
						}
					}
					else
						++numBgstClust;
				}
//			}

			--k;
			if (k == 0)
				break;
		}
		++i;
	}

//	fprintf(stderr, "bgstClustSize = %d, numBgstClust = %d\n", bgstClustSize,
//			numBgstClust);
//	for (i = 0; i < numBgstClust2; ++i)
//	{
//		fprintf(stderr, "   _refIdx[i] = %d, _refPos[i] = %d\n", _refIdxArr[i],
//				_refPosArr[i]);
//	}

	if (numBgstClust2 <= _maxHitsPerQry)
	{
		for (i = 0; i < _maxHitsPerQry; ++i)
		{
			if (i < numBgstClust2)
			{
				refIdx_bestHits[i] = _refIdxArr[i];
				refPos_bestHits[i] = _refPosArr[i];
			}
			else
				refIdx_bestHits[i] = -1;
		}
		numBgstClust = numBgstClust2;
	}
	else if (numBgstClust2 >= numBgstClust)
	{
		int *randNumArr = (int *) calloc(_maxHitsPerQry, sizeof(int));
		arrGetRandomNums2(_maxHitsPerQry, 0, numBgstClust2 - 1, randNumArr);
		for (i = 0; i < _maxHitsPerQry; ++i)
		{
			refIdx_bestHits[i] = _refIdxArr[randNumArr[i]];
			refPos_bestHits[i] = _refPosArr[randNumArr[i]];
		}
		free(randNumArr);
		numBgstClust = _maxHitsPerQry;
	}
	else
		numBgstClust = 0;

	return numBgstClust;
}


/**
 * Finds and returns the number of reference hits for tuples in a query.
 *
 * @param qrySeq	Query sequence
 * @param qryLen	Query sequence length.
 * @param refIdx	Reference index.
 * @param refPos	Reference position.
 * @param[out]	isContinous	A value of 0 indicates there is a mismatching
 * tuple; 1 indicates there are no mismatching tuples.
 * @return Number of reference hits for tuples in a query.
 */
static int getClustSize2(char *qrySeq, int qryLen, char refIdx, int refPos,
		int *isContinous)
{
	int clustSize = 0, clustSize2 = 0, hash, offset = 0;
	int numQryTuples =  qryLen - _seedLen;
	long refPosComposite = refIdx, refPosComposite2;
	refPosComposite = refPosComposite << REF_POS_BITS2;
	refPosComposite += (refPos / _seedLen);
	int refPosTmp = 0, refIdxTmp;
	int i, j, haveEarlierTuplesMatched = 0;
	*isContinous = 1;
	while (offset < numQryTuples)
	{
		hash = getHash(qrySeq + offset, _seedLen);
		refPosTmp = (refPosComposite & 33554431) * _seedLen;
		refIdxTmp = refPosComposite >> REF_POS_BITS2;
		if (_maxRefPosComposite > refPosComposite
				&& _refPosHashMap[refPosComposite] == hash)
		{
			++clustSize;
			offset += _seedLen;
			haveEarlierTuplesMatched = 1;

//			refPosComposite2 = refPosComposite - 1;
//			for (j = _numPrevHashes - 1; j >= 0; j = j - _seedLen)
//			{
//				if (refPosComposite2 >= 0
//						&& _refPosHashMap[refPosComposite2] == _prevHashes[j])
//					++clustSize2;
//				else if (refPosComposite2 < 0)
//					break;
//				--refPosComposite2;
//			}
//			if (clustSize2 > 0)
//			{
//				clustSize += clustSize2;
//				_numPrevHashes = 0;
//				clustSize2 = 0;
//			}
		}
		else if (haveEarlierTuplesMatched == 0)
			++offset;
		else
		{
			offset += _seedLen;
			*isContinous = 0;
		}
		++refPosComposite;
	}

	return clustSize;
}
