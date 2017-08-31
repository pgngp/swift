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
 * Implements CPU-based filtering using bins.
 */

#include "common.h"
#include "lookupTable10.h"
#include "preprocess.h"
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include "array.h"
#include <time.h>

#define	REF_POS_BITS2	25
#define	REF_POS_MASK2	33554431 /* Binary: 1111111111111111111111111 */
#define	TOTAL_NUM_CHRS	25
#define	CHR_SIZE_MARGIN	200

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
static int _chrSizes[MAX_NUM_REFS];
static int _numClustPerChr[MAX_NUM_REFS];
static int _numClust = 0;
static uchar *_clustSizes = NULL;
static int _refClustMap[MAX_NUM_REFS];
static char *_clustRefIdxMap = NULL;
static int *_clustRefPosMap = NULL;
static int _netClustSpan;
static int *_bestMatches = NULL;


/**
 * Creates clusters.
 */
static void createClust()
{
	_chrSizes[0] = CHR1_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[1] = CHR2_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[2] = CHR3_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[3] = CHR4_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[4] = CHR5_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[5] = CHR6_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[6] = CHR7_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[7] = CHR8_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[8] = CHR9_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[9] = CHR10_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[10] = CHR11_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[11] = CHR12_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[12] = CHR13_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[13] = CHR14_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[14] = CHR15_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[15] = CHR16_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[16] = CHR17_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[17] = CHR18_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[18] = CHR19_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[19] = CHR20_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[20] = CHR21_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[21] = CHR22_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[22] = CHRX_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[23] = CHRY_MAX_SIZE + CHR_SIZE_MARGIN;
	_chrSizes[24] = CHRMT_MAX_SIZE + CHR_SIZE_MARGIN;

	_netClustSpan = CLUST_SPAN - CLUST_OVERLAP;
	int i;
	for (i = 0; i < MAX_NUM_REFS; ++i)
	{
		_numClustPerChr[i] = (int) ceil(((float) _chrSizes[i]) / _netClustSpan);
		_refClustMap[i] = _numClust;
		_numClust += _numClustPerChr[i];
	}
	fprintf(stderr, "Number of clusters = %d\n", _numClust);
	_clustSizes = (uchar *) calloc(_numClust, sizeof(uchar));
	_clustRefIdxMap = (char *) calloc(_numClust, sizeof(char));
	_clustRefPosMap = (int *) calloc(_numClust, sizeof(int));
	int refPos = 0 - _netClustSpan;
	int j, k = 0, numClustPerChr;
	for (i = 0; i < MAX_NUM_REFS; ++i)
	{
		numClustPerChr = _numClustPerChr[i];
		for (j = 0; j < numClustPerChr; ++j)
		{
			_clustRefIdxMap[k] = i;
			refPos += _netClustSpan;
			_clustRefPosMap[k] = refPos;
			++k;
		}
		refPos = 0 - _netClustSpan;
	}
	_bestMatches = (int *) calloc(100000, sizeof(int));

//	for (i = 0; i < _numClust; ++i)
//		fprintf(stderr, "i = %d, refIdx = %d, refPos = %d\n", i,
//				_clustRefIdxMap[i], _clustRefPosMap[i]);
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
void lookupTable10Create(const char *refFile, int seedLen, int maxHitsPerQry,
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

	createClust();

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
void lookupTable10Create(const char *keysFile, const char *valsFile, int maxHits,
		int seedLen, int tupleIgnoreThres)
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

	createClust();
	time(&endTime);
	diffTime = difftime(endTime, startTime);
	fprintf(stderr, "(Time = %.2lf)...", diffTime);
}


/**
 * Reset key variables in this file.
 */
void lookupTable10Reset()
{

}


/**
 * Releases memory occupied by data structures in this file.
 */
void lookupTable10Delete()
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
	free(_clustSizes);
	_clustSizes = NULL;
	free(_clustRefIdxMap);
	_clustRefIdxMap = NULL;
	free(_clustRefPosMap);
	_clustRefPosMap = NULL;
	free(_bestMatches);
	_bestMatches = NULL;
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
int lookupTable10MapQry(char *qrySeq, int qryLen, char *refIdx_bestHits,
		int *refPos_bestHits, int *refSize_bestHits)
{
#define	NUM_BGST_CLUST	5
	int i, j, hash, startIdx, endIdx, clustIdx, numBgstClust = 0;
	int numQryTuples = qryLen - _seedLen + 1;
	short bgstClustSize = 0;
	int bgstClustIdxArr[NUM_BGST_CLUST], numBgstClust2 = 0;
	short *clustSizes = (short *) calloc(_numClust, sizeof(short));
	for (i = 0; i < numQryTuples; ++i)
	{
		hash = getHash(qrySeq + i, _seedLen);
		if (_numActualRepeatsPerTuple[hash] == 0)
			continue;

		startIdx = _lookupTable[hash];
		endIdx = startIdx + _numActualRepeatsPerTuple[hash] - 1;
		for (j = startIdx; j <= endIdx; ++j)
		{
			clustIdx = _refClustMap[_refIdx[j]] + (_refPos[j] / _netClustSpan);
			++clustSizes[clustIdx];
			if (bgstClustSize <= clustSizes[clustIdx])
			{
				if (bgstClustSize < clustSizes[clustIdx])
				{
					bgstClustSize = clustSizes[clustIdx];
					bgstClustIdxArr[0] = clustIdx;
					numBgstClust = 1;
					numBgstClust2 = 1;
				}
				else if (numBgstClust < NUM_BGST_CLUST)
				{
					bgstClustIdxArr[numBgstClust] = clustIdx;
					++numBgstClust;
					++numBgstClust2;
				}
				else
					++numBgstClust;
			}
		}
	}

//	for (i = 0; i < _numClust; ++i)
//		fprintf(stderr, "%d\n", clustSizes[i]);

	for (i = 0; i < _maxHitsPerQry; ++i)
		refIdx_bestHits[i] = -1;

	if (numBgstClust2 == 0)
	{
		free(clustSizes);
		return 0;
	}
	else if (numBgstClust2 <= _maxHitsPerQry)
	{
		for (i = 0; i < numBgstClust; ++i)
		{
			refIdx_bestHits[i] = _clustRefIdxMap[bgstClustIdxArr[i]];
			refPos_bestHits[i] = _clustRefPosMap[bgstClustIdxArr[i]];
			refSize_bestHits[i] = _chrSizes[refIdx_bestHits[i]] - CHR_SIZE_MARGIN;
		}
	}
	else if (numBgstClust2 >= numBgstClust)
	{
		int *randNumArr = (int *) calloc(_maxHitsPerQry, sizeof(int));
		arrGetRandomNums2(_maxHitsPerQry, 0, numBgstClust2 - 1, randNumArr);
		for (i = 0; i < _maxHitsPerQry; ++i)
		{
			refIdx_bestHits[i] = _clustRefIdxMap[bgstClustIdxArr[randNumArr[i]]];
			refPos_bestHits[i] = _clustRefPosMap[bgstClustIdxArr[randNumArr[i]]];
			refSize_bestHits[i] = _chrSizes[refIdx_bestHits[i]] - CHR_SIZE_MARGIN;
		}
		free(randNumArr);
		numBgstClust = _maxHitsPerQry;
	}
	else
	{
		j = 0;
		for (i = 0; i < _numClust; ++i)
		{
			if (clustSizes[i] == bgstClustSize)
			{
				_bestMatches[j] = i;
				++j;
			}
			if (j >= 100000)
			{
				free(clustSizes);
				return 0;
			}
		}
		numBgstClust = j;
		int *randNumArr = (int *) calloc(_maxHitsPerQry, sizeof(int));
		arrGetRandomNums2(_maxHitsPerQry, 0, numBgstClust - 1, randNumArr);
		for (i = 0; i < _maxHitsPerQry; ++i)
		{
			refIdx_bestHits[i] = _clustRefIdxMap[_bestMatches[randNumArr[i]]];
			refPos_bestHits[i] = _clustRefPosMap[_bestMatches[randNumArr[i]]];
			refSize_bestHits[i] = _chrSizes[refIdx_bestHits[i]] - CHR_SIZE_MARGIN;
		}
		free(randNumArr);
		numBgstClust = _maxHitsPerQry;
	}
	free(clustSizes);

	return numBgstClust;
}
