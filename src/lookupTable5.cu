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
 * Contains CPU-based filtering; based on SSAHA
 */

#include "common.h"
#include "lookupTable5.h"
#include "preprocess.h"
#include "mapHits6.h"
#include <stdio.h>
#include <limits.h>
#include <time.h>


#define	REF_POS_BITS2	25
#define	REF_POS_MASK2	33554431 /* Binary: 1111111111111111111111111 */


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
__constant__ int hashCodes_gpu[NUM_ASCII_CHARS];
__constant__ int powerVals_gpu[MAX_SEED_LENGTH];
extern __shared__ int arr_shr[];
__shared__ uint randNumSeed;


/**
 * This is a wrapper function that wraps @a lookupTable5Create function.
 * This function has been added so that @a lookupTable5Create can be
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
void lookupTable5Create_wrap(const char *refFile, int seedLen,
		int maxHitsPerQry, int **lookupTable, char **refIdx, int **refPos,
		int *numDistinctTuples, int **numRepeatsPerTuple,
		int tupleIgnoreThreshold)
{
	int numTotalTuples;
	lookupTable5Create(refFile, seedLen, maxHitsPerQry, tupleIgnoreThreshold,
			&numTotalTuples);
	*lookupTable = _lookupTable;
	*refIdx = _refIdx;
	*refPos = _refPos;
	*numDistinctTuples = _numDistinctTuples;
	*numRepeatsPerTuple = _numRepeatsPerTuple;
}


/**
 * Creates a lookup table.
 *
 * @param		refFile					Reference file.
 * @param		seedLen					Seed length.
 * @param 		maxHitsPerQry			Maximum hits per query.
 * @param		tupleIgnoreThreshold	Tuple that have number of repeats higher
 * than this value will be ignored.
 * @param[out]	totalTuples				Total number of reference tuples that
 * will be used.
 */
void lookupTable5Create(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold, int *totalTuples)
{
	_seedLen = seedLen;
	_maxHitsPerQry = maxHitsPerQry;
	_tupleIgnoreThreshold = tupleIgnoreThreshold;
	_numDistinctTuples = (int) pow((float) DNA_ALPHABET_SIZE, (int) _seedLen);
	_lookupTable = (int *) calloc(_numDistinctTuples, sizeof(int));
	int *numRepeatsPerTuple = (int *) calloc(_numDistinctTuples, sizeof(int));
	_numRepeatsPerTuple = (int *) calloc(_numDistinctTuples, sizeof(int));
	_numActualRepeatsPerTuple = (int *) calloc(_numDistinctTuples, sizeof(int));
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
//			for (i = 0; i < numIterations; i += 7)
//			for (i = 0; i < numIterations; i += 3)
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
	*totalTuples = _numTotalTuples;
	_refPos = (int *) calloc((_numTotalTuples + 1), sizeof(int));
	_refPos[0] = -1; /* First element set to -1 so that
	tuples that do not exist in the reference or that have number of repeats
	greater than a threshold value can point to this element. */
	_refIdx = (char *) calloc((_numTotalTuples + 1), sizeof(int));
	_refIdx[0] = -1;
	time(&endTime);
	diffTime = difftime(endTime, startTime);
	fprintf(stderr, "done. (Time = %.2lf secs)\n", diffTime);
	fclose(tmpRefFilePtr);

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
//	rewind(tmpRefFilePtr);
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

	int numQryTuples = MAX_QRY_SEQ_LENGTH - seedLen;
	_maxRefTuplesPerQry = numQryTuples * _tupleIgnoreThreshold;

//	FILE *someFile = fopen("/home/pgupta/data/gpusw/tupleRepeats.txt", "w");
//	char unhash[_seedLen + 1];
//	for (i = 0; i < _numDistinctTuples; ++i)
//	{
//		getUnhash(i, unhash, _seedLen);
//		fprintf(someFile, "%s\t%d\n", unhash, _numRepeatsPerTuple[i]);
//	}
//	fclose(someFile);
}


/**
 * Creates a lookup table.
 *
 * @param		refFile					Reference file.
 * @param		seedLen					Seed length.
 * @param 		maxHitsPerQry			Maximum hits per query.
 * @param		tupleIgnoreThreshold	Tuple that have number of repeats higher
 * than this value will be ignored.
 * @param[out]	totalTuples				Total number of reference tuples that
 * will be used.
 */
void lookupTable5Create2(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold, int *totalTuples)
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
	FILE *tmpRefFilePtr = fopen(tmpRefFile, "w+");
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
	*totalTuples = _numTotalTuples;
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

	int numQryTuples = MAX_QRY_SEQ_LENGTH - seedLen;
	_maxRefTuplesPerQry = numQryTuples * _tupleIgnoreThreshold;
}


/**
 * Reset key variables in this file.
 */
void lookupTable5Reset()
{

}


/**
 * Releases memory occupied by data structures in this file.
 */
void lookupTable5Delete()
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
 * Copies hash table from CPU to GPU.
 *
 * @param[out]	keys		Keys array on the GPU.
 * @param[out]	numKeys		Number of elements in @a keys array.
 * @param[out]	values		Values array on the GPU.
 * @param[out]	numValues	Number of elements in @a values array.
 * @param[out]	numRepeatsPerTuple	Number of tuples per hash.
 */
void lookupTable5CpyHashTableToGPU(int **keys, int *numKeys, int **values,
		int *numValues, int **numRepeatsPerTuple)
{
	cudaMalloc((void **) keys, _numDistinctTuples * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(*keys, _lookupTable, _numDistinctTuples * sizeof(int),
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()
	*numKeys = _numDistinctTuples;

	cudaMalloc((void **) values, (_numTotalTuples + 1) * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(*values, _refPos2, (_numTotalTuples + 1) * sizeof(int),
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()
	*numValues = _numTotalTuples + 1;

	cudaMalloc((void **) numRepeatsPerTuple, _numDistinctTuples * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(*numRepeatsPerTuple, _numActualRepeatsPerTuple,
			_numDistinctTuples * sizeof(int), cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()
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
void lookupTable5Create(const char *keysFile, const char *valsFile, int maxHits,
		int seedLen, int tupleIgnoreThres)
{
	FILE *keysFilePtr = fopen(keysFile, "r");
	char line[200];
	int numKeys = 0;
	while (fgets(line, 200, keysFilePtr) != NULL)
		++numKeys;
//	rewind(keysFilePtr);
	fclose(keysFilePtr);

	keysFilePtr = fopen(keysFile, "r");
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
//	rewind(valsFilePtr);
	fclose(valsFilePtr);

	valsFilePtr = fopen(valsFile, "r");
	_refIdx = (char *) calloc((numVals + 1), sizeof(char));
	_refPos = (int *) calloc((numVals + 1), sizeof(int));
	_refPos2 = (int *) calloc((numVals + 1), sizeof(int));
	i = 0;
	while (fgets(line, 200, valsFilePtr) != NULL)
	{
		sscanf(line, "%*d\t%d\t%d", &_refIdx[i], &_refPos[i]);
//		_refIdx[i] = (char) (_refPos2[i] >> REF_POS_BITS2);
//		_refPos[i] = (_refPos2[i] & REF_POS_MASK2) * seedLen;
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
 * @param	refIdx_bestHits		Best-matching reference indexes.
 * @param	shift_bestHits		Best-matching shifts.
 * @param	refPos_bestHits		Best-matching reference positions.
 * @return 	Number of hits in @a refIdx.
 */
int lookupTable5MapQry(const Query *query, uint qryLen, uint isRevComp,
		char *refIdx_bestHits, int *shift_bestHits, int *refPos_bestHits)
{
	int numQryTuples, i, j, k = 0, hash, startIdx, endIdx, numHits;
	char *seq;

	if (isRevComp == 1)
		seq = query->revCompSeq;
	else
		seq = query->seq;

	char *refIdx = (char *) calloc(_maxRefTuplesPerQry, sizeof(char));
	int *shift = (int *) calloc(_maxRefTuplesPerQry, sizeof(int));
	int *refPos = (int *) calloc(_maxRefTuplesPerQry, sizeof(int));

	/* Step 1: Break the query into tuples. */
	numQryTuples = qryLen - _seedLen + 1;
	for (i = 0; i < numQryTuples; ++i)
	{
		/* Step 2: Calculate the hash for each tuple. */
		/* Step 3: Find the reference coordinates of the hash in the hash
		 * table. */
		/* Step 4: Aggregate the reference coordinate positions. */
		hash = getHash(seq + i, _seedLen);
		if (_numActualRepeatsPerTuple[hash] == 0)
			continue;

		/* Add the tuples to the list. */
		startIdx = _lookupTable[hash];
		endIdx = startIdx + _numActualRepeatsPerTuple[hash] - 1;
		for (j = startIdx; j <= endIdx; ++j)
		{
			refIdx[k] = _refIdx[j];
			shift[k] = _refPos[j] - i;
			refPos[k] = _refPos[j];
			++k;
		}
	}
	numHits = mapHits6GetBestHits(_maxHitsPerQry, refIdx, shift, refPos, k,
			refIdx_bestHits, shift_bestHits, refPos_bestHits);

	free(refIdx);
	free(shift);
	free(refPos);

	return numHits;
}


/**
 * Searches the query sequence in the reference and returns the best-matching
 * reference coordinates.
 *
 * @param	query 				Query object.
 * @param 	seqLength 			Sequence length.
 * @param	refIdx_bestHits		Best-matching reference indexes.
 * @param	shift_bestHits		Best-matching shifts.
 * @param	refPos_bestHits		Best-matching reference positions.
 * @return 	Number of hits in @a refIdx.
 */
int lookupTable5MapQry2(char *qrySeq, int qryLen, char *refIdx_bestHits,
		int *shift_bestHits, int *refPos_bestHits)
{
	int numQryTuples, i, j, k = 0, hash, startIdx, endIdx, numHits;
	char *refIdx = (char *) calloc(_maxRefTuplesPerQry, sizeof(char));
	int *shift = (int *) calloc(_maxRefTuplesPerQry, sizeof(int));
	int *refPos = (int *) calloc(_maxRefTuplesPerQry, sizeof(int));

	/* Step 1: Break the query into tuples. */
	numQryTuples = qryLen - _seedLen + 1;
	for (i = 0; i < numQryTuples; ++i)
	{
		/* Step 2: Calculate the hash for each tuple. */
		/* Step 3: Find the reference coordinates of the hash in the hash
		 * table. */
		/* Step 4: Aggregate the reference coordinate positions. */
		hash = getHash(qrySeq + i, _seedLen);
		if (_numRepeatsPerTuple[hash] > _tupleIgnoreThreshold)
			continue;

		/* Add the tuples to the list. */
		startIdx = _lookupTable[hash];
		endIdx = startIdx + _numRepeatsPerTuple[hash] - 1;
		for (j = startIdx; j <= endIdx; ++j)
		{
			refIdx[k] = _refIdx[j];
			shift[k] = _refPos[j] - i;
			refPos[k] = _refPos[j];
			++k;
		}
	}
	numHits = mapHits6GetBestHits(_maxHitsPerQry, refIdx, shift, refPos, k,
			refIdx_bestHits, shift_bestHits, refPos_bestHits);

	free(refIdx);
	free(shift);
	free(refPos);

	return numHits;
}


/**
 * Copies constant memory from CPU to GPU.
 */
void lookupTable5CpyConstMemToGPU()
{
	int hashCodes[NUM_ASCII_CHARS];
	int i;
	for (i = 0; i < NUM_ASCII_CHARS; ++i)
		hashCodes[i] = 0;
	hashCodes['A'] = CODE_A;
	hashCodes['a'] = CODE_A;
	hashCodes['C'] = CODE_C;
	hashCodes['c'] = CODE_C;
	hashCodes['G'] = CODE_G;
	hashCodes['g'] = CODE_G;
	hashCodes['T'] = CODE_T;
	hashCodes['t'] = CODE_T;
	hashCodes['N'] = CODE_N;
	hashCodes['n'] = CODE_N;
	cudaMemcpyToSymbol(hashCodes_gpu, hashCodes, sizeof(hashCodes), 0,
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	int powerVals[MAX_SEED_LENGTH];
	for (i = 0; i < MAX_SEED_LENGTH; ++i)
		powerVals[i] = (int) pow((float) DNA_ALPHABET_SIZE, i);
	cudaMemcpyToSymbol(powerVals_gpu, powerVals, sizeof(powerVals), 0,
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()
}


/**
 * Searches the query sequence in the reference on the GPU.
 *
 * @param		keys	Keys of the lookup table.
 * @param		values	Values of the lookup table.
 * @param		numRptsPerTuple	Number of repeats per tuple.
 * @param		qrs		Query sequences.
 * @param		qryLen	Length of each query.
 * @param[out]	refIdx	Indices of the best matching reference sequences.
 * @param[out]	refPos	Positions of the best matching reference sequences.
 * @param		maxHits	Max number of hits.
 * @param		randNum	Random number.
 * @param		arrSize	Size of the arrays in shared memory.
 *
 * @note		Make sure the number of threads launched per block is less than
 * or equal to the number of elements in the arrays in the shared memory.
 */
__global__ void lookupTable5MapQry_gpu(int *keys, int *values,
		int *numRptsPerTuple, char *qrs, uchar *qryLen, char *refIdx,
		int *refPos, int maxHits, int seedLen, int randNum, int arrSize)
{
	int blockId = (blockIdx.y * gridDim.x) + blockIdx.x;
	short threadId = (threadIdx.y * blockDim.x) + threadIdx.x;
	int *randNums_shr = (int *) arr_shr;
	int *shift_shr = (int *) &randNums_shr[maxHits];
	int *refPos_shr = (int *) &shift_shr[arrSize];
	int *clusterSize_shr = (int *) &refPos_shr[arrSize];
	int *hash_shr = (int *) clusterSize_shr;
	char *refIdx_shr = (char *) &clusterSize_shr[arrSize];
	short numQryTuples = qryLen[blockId] - seedLen + 1;
	int i;
	randNumSeed = blockId + randNum;

	/* Initialize shared memory. */
	short totalThreads = blockDim.x * blockDim.y * blockDim.z;
	short binSize = (short) ceil(((float) (arrSize)) / totalThreads);
	short idx = threadId * binSize;
	for (i = 0; i < binSize; ++i)
	{
		if (idx < arrSize)
		{
			refIdx_shr[idx] = CHAR_MAX;
			shift_shr[idx] = -1;
			refPos_shr[idx] = -1;
			clusterSize_shr[idx] = -1;
		}
		++idx;
	}
	__syncthreads();

	/* Calculate query tuple hash. */
	if (threadId < numQryTuples)
	{
		int hash = 0;
		int idx = (blockId * MAX_QRY_SEQ_LENGTH) + threadId;
		for (i = 0; i < seedLen; ++i)
			hash += powerVals_gpu[i] * hashCodes_gpu[qrs[idx + i]];
		hash_shr[threadId] = hash;
	}
	__syncthreads();

	/* Fetch query tuple hits from the global memory. */
	if (threadId == 0)
	{
		short numRefTuples = 0;
		int begin, end, hash, j, value;
		for (i = 0; i < numQryTuples; ++i)
		{
			hash = hash_shr[i];
			begin = keys[hash];
			end = begin + numRptsPerTuple[hash] - 1;
			for (j = begin; j <= end; ++j)
			{
				value = values[j];
				refIdx_shr[numRefTuples] = (char) (value >> REF_POS_BITS2);
				refPos_shr[numRefTuples] = (value & REF_POS_MASK2) * seedLen;
				shift_shr[numRefTuples] = refPos_shr[numRefTuples] - i;
				++numRefTuples;

				/* If tuples cannot be fit into the shared memory, return. */
				if (numRefTuples >= arrSize)
					return;
			}
		}
	}
	__syncthreads();

	/* Sort query tuples (parallel merge-sort). */
	char tmpChar, isSorted = FALSE;
	int tmpInt;
	short shrMemBegin, shrMemEnd;
	short numThreads = totalThreads / 2;
	while (numThreads >= 1)
	{
		if (threadId < numThreads)
		{
			binSize = (short) ceil(((float) (arrSize)) / numThreads);
			shrMemBegin = threadId * binSize;
			shrMemEnd = shrMemBegin + binSize - 1;
			isSorted = FALSE;
			while (isSorted == FALSE)
			{
				isSorted = TRUE;
				for (i = shrMemBegin + 1; i <= shrMemEnd; ++i)
				{
					if ((i < arrSize)
							&& ((refIdx_shr[i] < refIdx_shr[i - 1])
							|| ((refIdx_shr[i] == refIdx_shr[i - 1])
									&& (shift_shr[i] < shift_shr[i - 1]))
							|| ((refIdx_shr[i] == refIdx_shr[i - 1])
									&& (shift_shr[i] == shift_shr[i - 1]))
									&& (refPos_shr[i] < refPos_shr[i - 1])))
					{
						tmpChar = refIdx_shr[i];
						refIdx_shr[i] = refIdx_shr[i - 1];
						refIdx_shr[i - 1] = tmpChar;

						tmpInt = shift_shr[i];
						shift_shr[i] = shift_shr[i - 1];
						shift_shr[i - 1] = tmpInt;

						tmpInt = refPos_shr[i];
						refPos_shr[i] = refPos_shr[i - 1];
						refPos_shr[i - 1] = tmpInt;
						isSorted = FALSE;
					} /* end-if */
				} /* end-for */
			} /* end-while */
		} /* end-if */
		numThreads = numThreads / 2;
		__syncthreads();
	} /* end-while */

	/* Create clusters, find biggest clusters, and choose randomly among
	 * biggest clusters. */
	if (threadId == 0)
	{
		/* Create clusters. */
		short biggestClusterSize = 1;
		int numClusters = 0;
		clusterSize_shr[numClusters] = 1;
		++numClusters;
		for (i = 1; i < arrSize; ++i)
		{
			if (refIdx_shr[i] != CHAR_MAX)
				break;
			else if ((refIdx_shr[i] == refIdx_shr[i - 1])
					&& (shift_shr[i] == shift_shr[i - 1]))
				++clusterSize_shr[numClusters - 1];
			else
			{
				biggestClusterSize = max(biggestClusterSize,
						clusterSize_shr[numClusters - 1]);
				clusterSize_shr[numClusters] = 1;
				++numClusters;
			}
		}
		biggestClusterSize = max(biggestClusterSize,
				clusterSize_shr[numClusters - 1]);

		/* Find biggest clusters. */
		char numBiggestHits = 0;
		int arrIdx = 0;
		for (i = 0; i < numClusters; ++i)
		{
			if (biggestClusterSize == clusterSize_shr[i])
			{
				shift_shr[numBiggestHits] = arrIdx;
				++numBiggestHits;
			}
			arrIdx += clusterSize_shr[i];
		}

		/* Choose randomly among biggest clusters if number of biggest clusters
		 * is greater than the maximum allowed number. */
		int globalMemIdx = blockId * maxHits;
		if (numBiggestHits <= maxHits)
		{
			for (i = 0; i < numBiggestHits; ++i)
			{
				refIdx[globalMemIdx + i] = refIdx_shr[shift_shr[i]];
				refPos[globalMemIdx + i] = refPos_shr[shift_shr[i]];
			}
		}
		else
		{
			/* Randomly choose the required number of hits from among all the
			 * hits. */
			arrGetRandomNums_gpu(maxHits, 0, numBiggestHits - 1, randNums_shr);
			for (i = 0; i < maxHits; ++i)
			{
				refIdx[globalMemIdx + i] = refIdx_shr[shift_shr[randNums_shr[i]]];
				refPos[globalMemIdx + i] = refPos_shr[shift_shr[randNums_shr[i]]];
			}
		}
	}
}


/**
 * Searches the query sequence in the reference on the GPU.
 *
 * @param		keys	Keys of the lookup table.
 * @param		values	Values of the lookup table.
 * @param		numRptsPerTuple	Number of repeats per tuple.
 * @param		qrs		Query sequences.
 * @param		qryLen	Length of each query.
 * @param		maxQrySeqLen	Max query sequence length.
 * @param[out]	refIdx	Indices of the best matching reference sequences.
 * @param[out]	refPos	Positions of the best matching reference sequences.
 * @param		maxHits	Max number of hits.
 * @param		randNum	Random number.
 * @param		arrSize	Size of the arrays in shared memory.
 *
 * @note		Make sure the number of threads launched per block is less than
 * or equal to the number of elements in the arrays in the shared memory.
 */
__global__ void lookupTable5MapQry2_gpu(int *keys, int *values,
		int *numRptsPerTuple, char *qrs, uchar *qryLen, short maxQrySeqLen,
		char *refIdx, int *refPos, int maxHits, int seedLen, int randNum,
		int arrSize)
{
	int blockId = (blockIdx.y * gridDim.x) + blockIdx.x;
	short threadId = (threadIdx.y * blockDim.x) + threadIdx.x;
	int *randNums_shr = (int *) arr_shr;
	int *shift_shr = (int *) &randNums_shr[maxHits];
	int *bgstHits_shr = (int *) shift_shr;
	int *refPos_shr = (int *) &shift_shr[arrSize];
	int *clusterSize_shr = (int *) &refPos_shr[arrSize];
	int *hash_shr = (int *) clusterSize_shr;
	char *refIdx_shr = (char *) &clusterSize_shr[arrSize];
	short numQryTuples = qryLen[blockId] - seedLen + 1;
	randNumSeed = blockId + randNum;

	/* Initialize shared memory. */
//	initializeShrMem_gpu(refIdx_shr, shift_shr, refPos_shr, clusterSize_shr,
//			arrSize, threadId);

	/* Calculate query tuple hash. */
	short totalThreads = blockDim.x * blockDim.y * blockDim.z;
	short binSize = (short) ceil(((float) numQryTuples) / totalThreads);
	short i, threadId2 = threadId * binSize;
	int idx;
	for (i = 0; i < binSize; ++i)
	{
		if (threadId2 < numQryTuples)
		{
			idx = (blockId * maxQrySeqLen) + threadId2;
			hash_shr[threadId2] = getHash_gpu(qrs + idx, seedLen);
		}
		++threadId2;
	}
	__syncthreads();

	/* Fetch query tuple hits from the global memory. */
	if (threadId == 0)
	{
		cpyHitsFromGlobalToShr_gpu(refIdx_shr, shift_shr, refPos_shr, hash_shr,
				arrSize, keys, values, numRptsPerTuple, numQryTuples, seedLen);
	}
	__syncthreads();

	/* Sort query tuples (parallel merge-sort). */
	sort_gpu(refIdx_shr, shift_shr, refPos_shr, arrSize, threadId);

//	/* Create clusters, find biggest clusters, and choose randomly among
//	 * biggest clusters. */
//	if (threadId == 0)
//	{
//		short biggestClusterSize;
//		short numClusters = createClusters_gpu(refIdx_shr, shift_shr,
//				clusterSize_shr, arrSize, &biggestClusterSize);
//		char numBiggestHits = findBiggestClusters_gpu(numClusters,
//				biggestClusterSize, clusterSize_shr, bgstHits_shr);
//		assignResults_gpu(blockId, maxHits, numBiggestHits, refIdx_shr, refIdx,
//				refPos_shr, refPos, bgstHits_shr, randNums_shr);
//	}
}


/**
 * This is a wrapper function that wraps @a cpyHitsFromGlobalToShr_gpu
 * function. It has been added so that @a cpyHitsFromGlobalToShr_gpu can be
 * unit-tested.
 *
 * @param[out]	refIdx	Reference indices.
 * @param[out]	shift	Shifts (reference position - query position).
 * @param[out]	refPos	Reference positions.
 * @param		hashes	Array containing tuple hash values for a query.
 * @param		arrSize	Number of elements in @a refIdx.
 * @param		keys	Array containing the indices to @a values.
 * @param		values	Array containing reference positions.
 * @param		numRptsPerTuple	Array containing number of repeats per reference
 * tuple.
 * @param		numQryTuples	Number of query tuples for a query.
 * @param		seedLen	Seed length.
 */
__global__ void cpyHitsFromGlobalToShr_gpu_wrap(char *refIdx, int *shift,
		int *refPos, int *hashes, int arrSize, int *keys, int *values,
		int *numRptsPerTuple, short numQryTuples, int seedLen)
{
	cpyHitsFromGlobalToShr_gpu(refIdx, shift, refPos, hashes, arrSize, keys,
			values, numRptsPerTuple, numQryTuples, seedLen);
}


/**
 * Copies tuple hits from global memory to shared memory.
 *
 * @param[out]	refIdx	Reference indices.
 * @param[out]	shift	Shifts (reference position - query position).
 * @param[out]	refPos	Reference positions.
 * @param		hashes	Array containing tuple hash values for a query.
 * @param		arrSize	Number of elements in @a refIdx.
 * @param		keys	Array containing the indices to @a values.
 * @param		values	Array containing reference positions.
 * @param		numRptsPerTuple	Array containing number of repeats per reference
 * tuple.
 * @param		numQryTuples	Number of query tuples for a query.
 * @param		seedLen	Seed length.
 */
__device__ void cpyHitsFromGlobalToShr_gpu(char *refIdx, int *shift,
		int *refPos, int *hashes, int arrSize, int *keys, int *values,
		int *numRptsPerTuple, short numQryTuples, int seedLen)
{
	short numRefTuples = 0, i;
	int begin, end, hash, j, value;
	for (i = 0; i < numQryTuples; ++i)
	{
		hash = hashes[i];
		begin = keys[hash];
		end = begin + numRptsPerTuple[hash] - 1;
		for (j = begin; j <= end; ++j)
		{
			value = values[j];
			refIdx[numRefTuples] = (char) (value >> REF_POS_BITS2);
			refPos[numRefTuples] = (value & REF_POS_MASK2) * seedLen;
			shift[numRefTuples] = refPos[numRefTuples] - i;
			++numRefTuples;

			/* If tuples cannot be fit into the shared memory, return. */
			if (numRefTuples >= arrSize)
				return;
		}
	}
}


/**
 * This is a wrapper function that wraps @a initializeShrMem_gpu function. It
 * has been added so that @a initializeShrMem_gpu can be unit-tested.
 *
 * @param[out]	refIdx		Reference indices.
 * @param[out]	shift		Shifts (reference position - query position).
 * @param[out]	refPos		Reference positions.
 * @param[out]	clusterSize	Cluster sizes.
 * @param		arrSize		Number of elements in @a refIdx.
 */
__global__ void intializeShrMem_gpu_wrap(char *refIdx, int *shift, int *refPos,
		int *clusterSize, int arrSize)
{
	short threadId = (threadIdx.y * blockDim.x) + threadIdx.x;
	initializeShrMem_gpu(refIdx, shift, refPos, clusterSize, arrSize, threadId);
}


/**
 * Initialize shared memory.
 *
 * @param[out]	refIdx		Reference indices.
 * @param[out]	shift		Shifts (reference position - query position).
 * @param[out]	refPos		Reference positions.
 * @param[out]	clusterSize	Cluster sizes.
 * @param		arrSize		Number of elements in @a refIdx.
 * @param		threadId	Thread ID.
 */
__device__ void initializeShrMem_gpu(char *refIdx, int *shift, int *refPos,
		int *clusterSize, int arrSize, short threadId)
{
	short totalThreads = blockDim.x * blockDim.y * blockDim.z;
	short binSize = (short) ceil(((float) (arrSize)) / totalThreads);
	short idx = threadId * binSize;
	short i;
	for (i = 0; i < binSize; ++i)
	{
		if (idx < arrSize)
		{
			refIdx[idx] = CHAR_MAX;
			shift[idx] = -1;
			refPos[idx] = -1;
			clusterSize[idx] = -1;
		}
		++idx;
	}
	__syncthreads();
}


/**
 * This is a wrapper function that wraps @a assignResults_gpu. It has been
 * added so that @a assignResults_gpu can be unit-tested.
 *
 * @param	maxHits		Max number of allowed hits.
 * @param	numBgstHits	Number of biggest hits.
 * @param	refIdx		Reference indices.
 * @param	refIdx_global	Reference index array on the global memory.
 * @param	refPos		Reference positions.
 * @param	refPos_global	Reference positions on the global memory.
 * @param	bgstClust	Array containing the indices of biggest clusters.
 * @param	randNum		Random number.
 */
__global__ void assignResults_gpu_wrap(char maxHits, char numBgstHits,
		char *refIdx, char *refIdx_global, int *refPos, int *refPos_global,
		int *bgstClust, int randNum)
{
	randNumSeed = randNum;
	__shared__ int randNums[10];
	int blockId = (blockIdx.y * gridDim.x) + blockIdx.x;
	assignResults_gpu(blockId, maxHits, numBgstHits, refIdx, refIdx_global,
			refPos, refPos_global, bgstClust, randNums);
}


/**
 * Assigns results to global memory.
 *
 * If the number of biggest clusters is greater than the maximum allowed
 * biggest clusters, then this function will randomly choose among all
 * the biggest clusters.
 *
 * @param	blockId		Current block ID.
 * @param	maxHits		Max number of allowed hits.
 * @param	numBgstHits	Number of biggest hits.
 * @param	refIdx		Reference indices.
 * @param	refIdx_global	Reference index array on the global memory.
 * @param	refPos		Reference positions.
 * @param	refPos_global	Reference positions on the global memory.
 * @param	bgstClust	Array containing the indices of biggest clusters.
 * @param	randNums	Array that can be used for storing random numbers.
 */
__device__ void assignResults_gpu(int blockId, char maxHits, char numBgstHits,
		char *refIdx, char *refIdx_global, int *refPos, int *refPos_global,
		int *bgstClust, int *randNums)
{
	int globalMemIdx = blockId * maxHits;
	char i;
	if (numBgstHits <= maxHits)
	{
		for (i = 0; i < numBgstHits; ++i)
		{
			refIdx_global[globalMemIdx + i] = refIdx[bgstClust[i]];
			refPos_global[globalMemIdx + i] = refPos[bgstClust[i]];
		}
	}
	else
	{
		/* Randomly choose the required number of hits from among all the
		 * hits. */
		arrGetRandomNums_gpu(maxHits, 0, numBgstHits - 1, randNums);
		for (i = 0; i < maxHits; ++i)
		{
			refIdx_global[globalMemIdx + i] = refIdx[bgstClust[randNums[i]]];
			refPos_global[globalMemIdx + i] = refPos[bgstClust[randNums[i]]];
		}
	}
}


/**
 * This is a wrapper function that wraps @a findBiggestClusters_gpu function.
 * It has been added so that @a findBiggestClusters_gpu can be unit-tested.
 *
 * @param		numClust		Number of clusters.
 * @param		bgstClustSize	Biggest cluster size.
 * @param		clusterSize		Array containing cluster sizes.
 * @param[out]	bgstClust		Array containing biggest clusters.
 */
__global__ void findBiggestClusters_gpu_wrap(short numClust,
		short bgstClustSize, int *clusterSize, int *bgstClust,
		char *numBgstClust)
{
	*numBgstClust = findBiggestClusters_gpu(numClust, bgstClustSize,
			clusterSize, bgstClust);
}


/**
 * Find biggest clusters.
 *
 * @param		numClust		Number of clusters.
 * @param		bgstClustSize	Biggest cluster size.
 * @param		clusterSize		Array containing cluster sizes.
 * @param[out]	bgstClust		Array containing biggest clusters.
 * @return		Number of biggest clusters.
 */
__device__ char findBiggestClusters_gpu(short numClust, short bgstClustSize,
		int *clusterSize, int *bgstClust)
{
	char numBiggestHits = 0;
	short arrIdx = 0, i;
	for (i = 0; i < numClust; ++i)
	{
		if (bgstClustSize == clusterSize[i])
		{
			bgstClust[numBiggestHits] = arrIdx;
			++numBiggestHits;
		}
		arrIdx += clusterSize[i];
	}
	return numBiggestHits;
}


/**
 * This is a wrapper function that wraps @a createClusters_gpu function. It
 * has been added so that @a createClusters_gpu can be unit-tested.
 *
 * @param		refIdx			Reference indices.
 * @param		shift			Shifts (reference position - query position).
 * @param[out]	clusterSize		Cluster sizes.
 * @param		arrSize			Size of @a refIdx array.
 * @param[out]	bgstClustSize	Biggest cluster size.
 * @param[out]	numClusters		Number of clusters.
 */
__global__ void createClusters_gpu_wrap(char *refIdx, int *shift,
		int *clusterSize, int arrSize, short *bgstClustSize, short *numClusters)
{
	*numClusters = createClusters_gpu(refIdx, shift, clusterSize,
			arrSize, bgstClustSize);
}


/**
 * Creates clusters and returns the number of clusters.
 *
 * @note This function assumes that the input arrays are already sorted, first
 * by reference index and then by shift.
 *
 * @param		refIdx			Reference indices.
 * @param		shift			Shifts (reference position - query position).
 * @param[out]	clusterSize		Cluster sizes.
 * @param		arrSize			Size of @a refIdx array.
 * @param[out]	bgstClustSize	Biggest cluster size.
 * @return		Number of clusters.
 */
__device__ short createClusters_gpu(char *refIdx, int *shift, int *clusterSize,
		int arrSize, short *bgstClustSize)
{
	short numClusters = 0, i;
	clusterSize[numClusters] = 1;
	++numClusters;
	*bgstClustSize = 1;
	for (i = 1; i < arrSize; ++i)
	{
		if (refIdx[i] == CHAR_MAX)
			continue;
		else if ((refIdx[i] == refIdx[i - 1]) && (shift[i] == shift[i - 1]))
			++clusterSize[numClusters - 1];
		else
		{
			*bgstClustSize = max(*bgstClustSize, clusterSize[numClusters - 1]);
			clusterSize[numClusters] = 1;
			++numClusters;
		}
	}
	*bgstClustSize = max(*bgstClustSize, clusterSize[numClusters - 1]);
	return numClusters;
}


/**
 * This is a wrapper function that wraps @a sort_gpu function. It has been
 * added so that @a sort_gpu can be unit-tested.
 *
 * @param[in,out]	refIdx			Reference indices.
 * @param[in,out]	shift			Shifts (reference position - query position).
 * @param[in,out]	refPos			Reference positions.
 * @param			arrSize			Size of @a refIdx array.
 */
__global__ void sort_gpu_wrap(char *refIdx, int *shift, int *refPos, int arrSize)
{
	short threadId = (threadIdx.y * blockDim.x) + threadIdx.x;
	sort_gpu(refIdx, shift, refPos, arrSize, threadId);
}


/**
 * Sorts the given arrays in-place using parallel merge-sort. It first sorts
 * using the reference index, then shift, and finally by reference position.
 *
 * @param[in,out]	refIdx			Reference indices.
 * @param[in,out]	shift			Shifts (reference position - query position).
 * @param[in,out]	refPos			Reference positions.
 * @param			arrSize			Size of @a refIdx array.
 * @param			threadId		Current thread ID.
 */
//__device__ void sort_gpu(char *refIdx, int *shift, int *refPos, int arrSize,
//		short threadId)
//{
//	__shared__ char isSorted;
//	short numThreads = blockDim.x * blockDim.y * blockDim.z;
//	short binSize = (short) ceil(((float) arrSize) / numThreads);
//	if (binSize < 2)
//	{
//		numThreads = numThreads / 2;
//		binSize = (short) ceil(((float) arrSize) / numThreads);
//	}
//
//	if (threadId == 0)
//		isSorted = FALSE;
//	__syncthreads();
//
//	short i, begin, end;
//	char tmpChar;
//	int tmpInt, tmp = 0;
//	while (isSorted == FALSE)
//	{
//		if (threadId == 0)
//			isSorted = TRUE;
//		__syncthreads();
//
//		begin = (threadId * binSize) + (tmp % 2);
//		end = i + binSize - 1;
//		for (i = begin + 1; (i <= end) && (i < arrSize); ++i)
//		{
//			if ((refIdx[i] < refIdx[i - 1])
//					|| ((refIdx[i] == refIdx[i - 1])
//							&& (shift[i] < shift[i - 1]))
//					|| ((refIdx[i] == refIdx[i - 1])
//							&& (shift[i] == shift[i - 1])
//							&& (refPos[i] < refPos[i - 1]))
//							)
//			{
//				tmpChar = refIdx[i];
//				refIdx[i] = refIdx[i - 1];
//				refIdx[i - 1] = tmpChar;
//
//				tmpInt = shift[i];
//				shift[i] = shift[i - 1];
//				shift[i - 1] = tmpInt;
//
//				tmpInt = refPos[i];
//				refPos[i] = refPos[i - 1];
//				refPos[i - 1] = tmpInt;
//				isSorted = FALSE;
//			}
//		}
//		++tmp;
//		__syncthreads();
//	}
//}


/**
 * Sorts the given arrays in-place using parallel merge-sort. It first sorts
 * using the reference index, then shift, and finally by reference position.
 *
 * @param[in,out]	refPos			Reference positions.
 * @param			arrSize			Size of @a refIdx array.
 * @param			threadId		Current thread ID.
 */
__device__ void sort_gpu(char *refIdx, int *shift, int *refPos, int arrSize,
		short threadId)
{
	short numThreads = blockDim.x * blockDim.y * blockDim.z;
	short binSize = (short) ceil(((float) arrSize) / numThreads);
	if (binSize < 2)
	{
		numThreads = numThreads / 2;
		binSize = (short) ceil(((float) arrSize) / numThreads);
	}
	int numIter = (int) ceil(((float) arrSize) / binSize);

	char isSorted, tmpChar;
	short i, begin, end;
	int tmp = 0, tmpInt;
	while (numIter >= 0)
	{
		begin = (threadId * binSize) + (tmp % 2);
		end = begin + binSize - 1;
		if (end >= arrSize)
			end = arrSize - 1;
		isSorted = FALSE;
		while (isSorted == FALSE)
		{
			isSorted = TRUE;
			for (i = begin + 1; i <= end; ++i)
			{
				if ((refIdx[i] < refIdx[i - 1])
						|| ((refIdx[i] == refIdx[i - 1])
								&& (shift[i] < shift[i - 1]))
						|| ((refIdx[i] == refIdx[i - 1])
								&& (shift[i] == shift[i - 1])
								&& (refPos[i] < refPos[i - 1])))
				{
					tmpChar = refIdx[i];
					refIdx[i] = refIdx[i - 1];
					refIdx[i - 1] = tmpChar;

					tmpInt = shift[i];
					shift[i] = shift[i - 1];
					shift[i - 1] = tmpInt;

					tmpInt = refPos[i];
					refPos[i] = refPos[i - 1];
					refPos[i - 1] = tmpInt;
					isSorted = FALSE;
				}
			}
		}
		++tmp;
		--numIter;
		__syncthreads();
	}
}


///**
// * Sorts the given arrays in-place using parallel merge-sort. It first sorts
// * using the reference index, then shift, and finally by reference position.
// *
// * @param[in,out]	refIdx			Reference indices.
// * @param[in,out]	shift			Shifts (reference position - query position).
// * @param[in,out]	refPos			Reference positions.
// * @param			arrSize			Size of @a refIdx array.
// * @param			threadId		Current thread ID.
// */
//__device__ void sort_gpu(char *refIdx, int *shift, int *refPos, int arrSize,
//		short threadId)
//{
//	char tmpChar, isSorted;
//	short begin, end, binSize, i;
//	short numThreads = blockDim.x * blockDim.y * blockDim.z;
//	int tmpInt;
//	while (numThreads >= 1)
//	{
//		if (threadId < numThreads)
//		{
//			binSize = (short) ceil(((float) arrSize) / numThreads);
//			begin = threadId * binSize;
//			end = begin + binSize - 1;
//			isSorted = FALSE;
//			while (isSorted == FALSE)
//			{
//				isSorted = TRUE;
//				for (i = begin + 1; i <= end; ++i)
//				{
//					if (i >= arrSize)
//						break;
//					else if ((refIdx[i] < refIdx[i - 1])
//							|| ((refIdx[i] == refIdx[i - 1])
//									&& (shift[i] < shift[i - 1]))
//							|| ((refIdx[i] == refIdx[i - 1])
//									&& (shift[i] == shift[i - 1]))
//									&& (refPos[i] < refPos[i - 1]))
//					{
//						tmpChar = refIdx[i];
//						refIdx[i] = refIdx[i - 1];
//						refIdx[i - 1] = tmpChar;
//
//						tmpInt = shift[i];
//						shift[i] = shift[i - 1];
//						shift[i - 1] = tmpInt;
//
//						tmpInt = refPos[i];
//						refPos[i] = refPos[i - 1];
//						refPos[i - 1] = tmpInt;
//						isSorted = FALSE;
//					} /* end-if */
//				} /* end-for */
//			} /* end-while */
//		} /* end-if */
//		numThreads = numThreads / 2;
//		__syncthreads();
//	} /* end-while */
//}


/**
 * This is a wrapper function that wraps @a getHash_gpu function. It has been
 * added so that @a getHash_gpu can be unit-tested.
 *
 * @param		str		String for which hash is to be calculated.
 * @param		len		Length of @a str.
 * @param[out]	hash	Hash value of @a str.
 */
__global__ void getHash_gpu_wrap(char *str, int len, int *hash)
{
	*hash = getHash_gpu(str, len);
}


/**
 * Returns the hash value of the given string.
 *
 * @param	str	String for which hash is to be calculated.
 * @param	len	Length of @a str.
 * @return	Hash value.
 */
__device__ int getHash_gpu(char *str, int len)
{
	int hash = 0, i;
	for (i = 0; i < len; ++i)
		hash += powerVals_gpu[i] * hashCodes_gpu[str[i]];
	return hash;
}


/**
 * Returns the index of the given number if it is already present in the
 * given array; otherwise, returns -1.
 *
 * @param	arr		Array in which the number will be searched.
 * @param	arrSize	Number of elements in the array. This number must be
 * equal to the number of elements in @a arr, otherwise it may result in invalid
 * behavior.
 * @param	num		The number to be searched in the array.
 * @return 	Index of the searched number in the array; otherwise, -1.
 */
__device__ int arrSearch_gpu(int *arr, int arrSize, int num)
{
	int i = 0;
	for (i = 0; i < arrSize; ++i)
	{
		if (arr[i] == num)
			return i;
	}
	return -1;
}


/**
 * Fetches the given number of random numbers between the given limits.
 *
 * @param		n		Number of random numbers to be created. This number
 * should not be greater than the range specified by @a lowerLimit and
 * @a upperLimit.
 * @param		lLimit	The lower limit of the random numbers.
 * @param		uLimit	The upper limit of the random numbers.
 * @param[out] 	arr 	Array in which the random numbers will be stored.
 * The size of the array should be atleast @a n.
 */
__device__ void arrGetRandomNums_gpu(int n, int lLimit, int uLimit, int *arr)
{
	int range = uLimit - lLimit + 1;
	int i = 0;
	int randNum;
	while (i < n)
	{
		randNum = (getRandNum_gpu() % range) + lLimit;
		if (arrSearch_gpu(arr, i, randNum) == -1)
		{
			arr[i] = randNum;
			++i;
		}
	}
}


/**
 * Returns a pseudo-random number.
 *
 * @return	Random number.
 *
 * @note This algorithm has been taken from "The C Programming Language"
 * by Kernighan and Ritchie.
 */
__device__ int getRandNum_gpu()
{
	randNumSeed = (randNumSeed * 1103515245) + 12345;
	return ((randNumSeed / 65536) % 32768);
}


/**
 * Calculates the value of base raised to the n-th power.
 *
 * @param 	base	Base value.
 * @param 	n		Exponent value.
 * @return			The calculated value.
 */
__device__ int pow_gpu(int base, int n)
{
	int p = 1;
	while (n > 0)
	{
		p = p * base;
		--n;
	}
	return p;
}

