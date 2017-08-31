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
 * Contains GPU-based filtering.
 */

#include "common.h"
#include "lookupTable7.h"
#include "preprocess.h"
#include "mapHits6.h"
#include <stdio.h>
#include <limits.h>
#include <time.h>

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
__constant__ int hashCodes_gpu[NUM_ASCII_CHARS];
__constant__ int powerVals_gpu[MAX_SEED_LENGTH];
extern __shared__ int arr_shr[];
__shared__ uint randNumSeed;


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
void lookupTable7Create2(const char *refFile, int seedLen, int maxHitsPerQry,
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

	int numQryTuples = MAX_QRY_SEQ_LENGTH - seedLen;
	_maxRefTuplesPerQry = numQryTuples * _tupleIgnoreThreshold;
	fprintf(stderr, "Max reference tuples per query = %d\n", _maxRefTuplesPerQry);

//	FILE *keysFilePtr = fopen("/home/pgupta/data/gpusw/repeats.txt", "w");
//	char seed[seedLen + 1];
//	for (i = 0; i <= _numDistinctTuples; ++i)
//	{
//		getUnhash(i, seed, seedLen);
//		fprintf(keysFilePtr, "%s\t%d\t%d\n", seed, i, _numRepeatsPerTuple[i]);
//	}
//	fclose(keysFilePtr);

//	FILE *valsFilePtr = fopen("/home/pgupta/data/gpusw/vals100.txt", "w");
//	for (i = 0; i <= _numTotalTuples; ++i)
//		fprintf(valsFilePtr, "%d\t%d\t%d\n", _refPos2[i], _refIdx[i], _refPos[i]);
//	fclose(valsFilePtr);
//
//	FILE *keysFilePtr = fopen("/home/pgupta/data/gpusw/keys100.txt", "w");
//	for (i = 0; i < _numDistinctTuples; ++i)
//		fprintf(keysFilePtr, "%d\t%d\n", _lookupTable[i],
//				_numActualRepeatsPerTuple[i]);
//	fclose(keysFilePtr);
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
void lookupTable7Create3(const char *keysFile, const char *valsFile,
		int maxHits, int seedLen, int tupleIgnoreThres)
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
}


/**
 * Reset key variables in this file.
 */
void lookupTable7Reset()
{

}


/**
 * Releases memory occupied by data structures in this file.
 */
void lookupTable7Delete()
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
void lookupTable7CpyHashTableToGPU(int **keys, int *numKeys, int **values,
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
 * Copies hash table from CPU to GPU.
 *
 * @param[out]	keys		Keys array on the GPU.
 * @param[out]	numKeys		Number of elements in @a keys array.
 * @param[out]	values		Values array on the GPU.
 * @param[out]	numValues	Number of elements in @a values array.
 * @param[out]	numRepeatsPerTuple	Number of tuples per hash.
 */
void lookupTable7FetchHashTable(int **keys, int *numKeys, int **values,
		int *numValues, int **numRepeatsPerTuple)
{
	*keys = _lookupTable;
	*numKeys = _numDistinctTuples;

	*values = _refPos2;
	*numValues = _numTotalTuples + 1;

	*numRepeatsPerTuple = _numActualRepeatsPerTuple;
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
int lookupTable7MapQry2(char *qrySeq, int qryLen, char *refIdx_bestHits,
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
 * Copies constant memory from CPU to GPU.
 */
void lookupTable7CpyConstMemToGPU()
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
 * Copies constant memory from CPU to GPU.
 */
static int hashCodes[NUM_ASCII_CHARS];
static int powerVals[MAX_SEED_LENGTH];
void lookupTable7CreateConstMem()
{
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

	for (i = 0; i < MAX_SEED_LENGTH; ++i)
		powerVals[i] = (int) pow((float) DNA_ALPHABET_SIZE, i);
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
 * @param		maxNumChr	Max number of chromosomes.
 *
 * @note		Make sure the number of threads launched per block is less than
 * or equal to the number of elements in the arrays in the shared memory.
 */
__global__ void lookupTable7MapQry2_gpu2(int *keys, int *values,
		int *numRptsPerTuple, char *qrs, uchar *qryLen, short maxQrySeqLen,
		char *refIdx, int *refPos, int maxHits, int seedLen, int randNum,
		int arrSize, int maxNumChr)
{
	int *shift_shr = (int *) arr_shr;
	int *refPos_shr = (int *) &shift_shr[arrSize];
	int *values_shr = (int *) &refPos_shr[arrSize];
	int *randNums_shr = (int *) &values_shr[arrSize];
	short *clusterSize_shr = (short *) &randNums_shr[maxHits];
	short *hitsPerChr_shr = (short *) &clusterSize_shr[arrSize];
	short *hitsStartIdx_shr = (short *) &hitsPerChr_shr[maxNumChr];
	char *refIdx_shr = (char *) &hitsStartIdx_shr[maxNumChr];

	int *distinctShift_shr = (int *) values_shr;
	int *bestMatches_shr = (int *) values_shr;
	short *qryTupleIdx_shr = (short *) clusterSize_shr;
	int blockId = (blockIdx.y * gridDim.x) + blockIdx.x;
	short threadId = (threadIdx.y * blockDim.x) + threadIdx.x;

	/* Fetch values from global memory and split it into reference index,
	 * shift and position. */
	__shared__ int numHits;
	if (threadId == 0)
	{
		int globalMemIdx = blockId * maxHits;
		short k;
		for (k = 0; k < maxHits; ++k)
		{
			refIdx[globalMemIdx + k] = -1;
			refPos[globalMemIdx + k] = -1;
		}

		/* Put values in shared memory. */
		int idx = blockId * maxQrySeqLen;
		short numQryTuples = qryLen[blockId] - seedLen + 1;
		int hash, i, j, begin, end, value;
		numHits = 0;
		for (i = 0; i < numQryTuples; ++i)
		{
			hash = getHash_gpu7(qrs + idx + i, seedLen);
			begin = keys[hash];
			end = begin + numRptsPerTuple[hash] - 1;

			for (j = begin; j <= end; ++j)
			{
				value = values[j];
				if (value < 0)
					continue;
				values_shr[numHits] = value;
				qryTupleIdx_shr[numHits] = i;
				++numHits;

				if (numHits >= arrSize)
					break;
			}
			if (numHits >= arrSize)
				break;
		}

		/* Split values into reference index, shift, and reference position. */
		if (numHits > 0 && numHits <= arrSize)
		{
			for (i = 0; i < maxNumChr; ++i)
				hitsPerChr_shr[i] = 0;
			char rIdx;
			for (i = 0; i < numHits; ++i)
			{
				rIdx = (char) (values_shr[i] >> REF_POS_BITS2);
				++hitsPerChr_shr[rIdx];
			}

			hitsStartIdx_shr[0] = 0;
			for (i = 1; i < maxNumChr; ++i)
				hitsStartIdx_shr[i] = hitsStartIdx_shr[i - 1]
				                                       + hitsPerChr_shr[i - 1];

			for (i = 0; i < arrSize; ++i)
			{
				refIdx_shr[i] = -1;
				shift_shr[i] = -1;
				refPos_shr[i] = -1;
			}

			short startIdx;
			for (i = 0; i < numHits; ++i)
			{
				rIdx = (char) (values_shr[i] >> REF_POS_BITS2);
				startIdx = hitsStartIdx_shr[rIdx];
				refIdx_shr[startIdx] = rIdx;
				refPos_shr[startIdx] = (values_shr[i] & REF_POS_MASK2) * seedLen;
				shift_shr[startIdx] = refPos_shr[startIdx] -
						((int) qryTupleIdx_shr[i]) + maxQrySeqLen;
				++hitsStartIdx_shr[rIdx];
			}

			for (i = 0; i < arrSize; ++i)
			{
				distinctShift_shr[i] = 0;
				clusterSize_shr[i] = 0;
			}

			hitsStartIdx_shr[0] = 0;
			for (i = 1; i < maxNumChr; ++i)
				hitsStartIdx_shr[i] = hitsStartIdx_shr[i - 1]
				                                       + hitsPerChr_shr[i - 1];
		}
	}
	__syncthreads();

	/* Create clusters. */
	if (threadId < maxNumChr && numHits > 0 && numHits <= arrSize)
	{
		short start = hitsStartIdx_shr[threadId];
		short end = start + hitsPerChr_shr[threadId] - 1;
		createClusters_gpu7_2(shift_shr, refPos_shr, start, end, clusterSize_shr,
				distinctShift_shr);
	}
	__syncthreads();

	/* Find the biggest cluster and pick the first hit of the biggest
	 * cluster as the result. */
	if (threadId == 0 && numHits > 0 && numHits <= arrSize)
	{
		/* Find biggest cluster size. */
		short bgstClusterSize = 0, i;
		for (i = 0; i < arrSize; ++i)
			bgstClusterSize = max(bgstClusterSize, clusterSize_shr[i]);

		/* Find biggest clusters. */
		short numBgstHits = 0;
		for (i = 0; i < arrSize; ++i)
		{
			if (clusterSize_shr[i] == bgstClusterSize)
			{
				bestMatches_shr[numBgstHits] = shift_shr[i];
				++numBgstHits;
			}
		}

		/* Assign the first hit of the biggest clusters to global memory. */
		int globalMemIdx = blockId * maxHits;
		if (numBgstHits <= maxHits)
		{
			for (i = 0; i < numBgstHits; ++i)
			{
				refIdx[globalMemIdx + i] = refIdx_shr[bestMatches_shr[i]];
				refPos[globalMemIdx + i] = refPos_shr[bestMatches_shr[i]];
			}
		}
		else
		{
			/* Randomly choose the required number of hits from among all the
			 * hits. */
			arrGetRandomNums_gpu7(maxHits, 0, numBgstHits - 1, randNums_shr);
			for (i = 0; i < maxHits; ++i)
			{
				refIdx[globalMemIdx + i] =
						refIdx_shr[bestMatches_shr[randNums_shr[i]]];
				refPos[globalMemIdx + i] =
						refPos_shr[bestMatches_shr[randNums_shr[i]]];
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
 * @param		maxNumChr	Max number of chromosomes.
 * @param		blockId		Block ID.
 *
 * @note		Make sure the number of threads launched per block is less than
 * or equal to the number of elements in the arrays in the shared memory.
 */
void lookupTable7MapQry2_cpu2(int *keys, int *values, int *numRptsPerTuple,
		char *qrs, uchar *qryLen, short maxQrySeqLen, char *refIdx,
		int *refPos, int maxHits, int seedLen, int randNum, int arrSize,
		int maxNumChr, int blockId)
{
	int *shift_shr = (int *) calloc(arrSize, sizeof(int));
	int *refPos_shr = (int *) calloc(arrSize, sizeof(int));
	int *values_shr = (int *) calloc(arrSize, sizeof(int));
	int *randNums_shr = (int *) calloc(maxHits, sizeof(int));
	short *clusterSize_shr = (short *) calloc(arrSize, sizeof(short));
	short *hitsPerChr_shr = (short *) calloc(maxNumChr, sizeof(short));
	short *hitsStartIdx_shr = (short *) calloc(maxNumChr, sizeof(short));
	char *refIdx_shr = (char *) calloc(arrSize, sizeof(char));

	int *distinctShift_shr = (int *) calloc(arrSize, sizeof(int));
	int *bestMatches_shr = (int *) calloc(arrSize, sizeof(int));
	short *qryTupleIdx_shr = (short *) calloc(arrSize, sizeof(short));

	int numHits;
	int globalMemIdx = blockId * maxHits;
	short k;
	for (k = 0; k < maxHits; ++k)
	{
		refIdx[globalMemIdx + k] = -1;
		refPos[globalMemIdx + k] = -1;
	}

	int idx = blockId * maxQrySeqLen;
	short numQryTuples = qryLen[blockId] - seedLen + 1;
	int hash, i, j, begin, end, value;
	numHits = 0;
	for (i = 0; i < numQryTuples; ++i)
	{
		hash = getHash_cpu7(qrs + idx + i, seedLen);
		begin = keys[hash];
		end = begin + numRptsPerTuple[hash] - 1;

		for (j = begin; j <= end; ++j)
		{
			value = values[j];
			if (value < 0)
				continue;
			values_shr[numHits] = value;
			qryTupleIdx_shr[numHits] = i;
			++numHits;

			if (numHits >= arrSize)
				break;
		}
		if (numHits >= arrSize)
			break;
	}

	if (numHits > 0 && numHits <= arrSize)
	{
		for (i = 0; i < maxNumChr; ++i)
			hitsPerChr_shr[i] = 0;
		char rIdx;
		for (i = 0; i < numHits; ++i)
		{
			rIdx = (char) (values_shr[i] >> REF_POS_BITS2);
			++hitsPerChr_shr[rIdx];
		}

		hitsStartIdx_shr[0] = 0;
		for (i = 1; i < maxNumChr; ++i)
			hitsStartIdx_shr[i] = hitsStartIdx_shr[i - 1]
												   + hitsPerChr_shr[i - 1];

		for (i = 0; i < arrSize; ++i)
		{
			refIdx_shr[i] = -1;
			shift_shr[i] = -1;
			refPos_shr[i] = -1;
		}

		short startIdx;
		for (i = 0; i < numHits; ++i)
		{
			rIdx = (char) (values_shr[i] >> REF_POS_BITS2);
			startIdx = hitsStartIdx_shr[rIdx];
			refIdx_shr[startIdx] = rIdx;
			refPos_shr[startIdx] = (values_shr[i] & REF_POS_MASK2) * seedLen;
			shift_shr[startIdx] = refPos_shr[startIdx] -
					((int) qryTupleIdx_shr[i]) + maxQrySeqLen;
			++hitsStartIdx_shr[rIdx];
		}

		for (i = 0; i < arrSize; ++i)
		{
			distinctShift_shr[i] = 0;
			clusterSize_shr[i] = 0;
		}

		hitsStartIdx_shr[0] = 0;
		for (i = 1; i < maxNumChr; ++i)
			hitsStartIdx_shr[i] = hitsStartIdx_shr[i - 1] + hitsPerChr_shr[i - 1];
	}

	if (numHits > 0 && numHits <= arrSize)
	{
		for (i = 0; i < maxNumChr; ++i)
		{
			short start = hitsStartIdx_shr[i];
			short end = start + hitsPerChr_shr[i] - 1;
			createClusters_cpu7_2(shift_shr, refPos_shr, start, end,
					clusterSize_shr, distinctShift_shr);
		}
	}

	if (numHits > 0 && numHits <= arrSize)
	{
		/* Find biggest cluster size. */
		short bgstClusterSize = 0, i;
		for (i = 0; i < arrSize; ++i)
			bgstClusterSize = max(bgstClusterSize, clusterSize_shr[i]);

		/* Find biggest clusters. */
		short numBgstHits = 0;
		for (i = 0; i < arrSize; ++i)
		{
			if (clusterSize_shr[i] == bgstClusterSize)
			{
				bestMatches_shr[numBgstHits] = shift_shr[i];
				++numBgstHits;
			}
		}

		/* Assign the first hit of the biggest clusters to global memory. */
		int globalMemIdx = blockId * maxHits;
		if (numBgstHits <= maxHits)
		{
			for (i = 0; i < numBgstHits; ++i)
			{
				refIdx[globalMemIdx + i] = refIdx_shr[bestMatches_shr[i]];
				refPos[globalMemIdx + i] = refPos_shr[bestMatches_shr[i]];
			}
		}
		else
		{
			/* Randomly choose the required number of hits from among all the
			 * hits. */
			arrGetRandomNums_cpu7(maxHits, 0, numBgstHits - 1, randNums_shr,
					randNum);
			for (i = 0; i < maxHits; ++i)
			{
				refIdx[globalMemIdx + i] =
						refIdx_shr[bestMatches_shr[randNums_shr[i]]];
				refPos[globalMemIdx + i] =
						refPos_shr[bestMatches_shr[randNums_shr[i]]];
			}
		}
	}

	free(shift_shr);
	free(refPos_shr);
	free(values_shr);
	free(randNums_shr);
	free(clusterSize_shr);
	free(hitsPerChr_shr);
	free(hitsStartIdx_shr);
	free(refIdx_shr);
	free(distinctShift_shr);
	free(bestMatches_shr);
	free(qryTupleIdx_shr);
}


/**
 * This is a wrapper function that wraps @a createClusters_gpu7_2 function. It
 * has been added so that @a createClusters_gpu7_2 can be unit-tested.
 *
 * @param[in,out]	shift		Shifts (reference position - query position).
 * @param			refPos		Reference positions.
 * @param			start		The starting index of the array.
 * @param			end			The ending index of the array.
 * @param[out]		clusterSize	Array where cluster sizes will be stored.
 * @param[out]		distinctShift	Array where distinct shifts will be stored.
 * @param			arrSize		Number of elements in the given arrays.
 */
__global__ void createClusters_gpu7_2_wrap(int *shift, int *refPos, short start,
		short end, short *clusterSize, int *distinctShift, int arrSize)
{
	int i;
	for (i = 0; i < arrSize; ++i)
	{
		clusterSize[i] = 0;
		distinctShift[i] = -1;
	}
	createClusters_gpu7_2(shift, refPos, start, end, clusterSize, distinctShift);
}


/**
 * Creates clusters.
 *
 * @param[in,out]	shift		Shifts (reference position - query position).
 * @param			refPos		Reference positions.
 * @param			start		The starting index of the array.
 * @param			end			The ending index of the array.
 * @param[out]		clusterSize	Array where cluster sizes will be stored.
 * @param			distinctShift	Array where distinct shifts will be stored.
 */
__device__ void createClusters_gpu7_2(int *shift, int *refPos, short start,
		short end, short *clusterSize, int *distinctShift)
{
	char isFound = FALSE;
	short numShifts = 0, i, j;
	for (i = start; i <= end; ++i)
	{
		isFound = FALSE;
		for (j = 0; j < numShifts; ++j)
		{
			if (shift[i] == distinctShift[start + j])
			{
				isFound = TRUE;
				++clusterSize[start + j];
				if (refPos[i] < refPos[shift[start + j]])
					shift[start + j] = i;
				break;
			}
		}
		if (isFound == FALSE)
		{
			clusterSize[start + numShifts] = 1;
			distinctShift[start + numShifts] = shift[i];
			shift[start + numShifts] = i;
			++numShifts;
		}
	}
}


/**
 * Creates clusters.
 *
 * @param[in,out]	shift		Shifts (reference position - query position).
 * @param			refPos		Reference positions.
 * @param			start		The starting index of the array.
 * @param			end			The ending index of the array.
 * @param[out]		clusterSize	Array where cluster sizes will be stored.
 * @param			distinctShift	Array where distinct shifts will be stored.
 */
void createClusters_cpu7_2(int *shift, int *refPos, short start, short end,
		short *clusterSize, int *distinctShift)
{
	char isFound = FALSE;
	short numShifts = 0, i, j;
	for (i = start; i <= end; ++i)
	{
		isFound = FALSE;
		for (j = 0; j < numShifts; ++j)
		{
			if (shift[i] == distinctShift[start + j])
			{
				isFound = TRUE;
				++clusterSize[start + j];
				if (refPos[i] < refPos[shift[start + j]])
					shift[start + j] = i;
				break;
			}
		}
		if (isFound == FALSE)
		{
			clusterSize[start + numShifts] = 1;
			distinctShift[start + numShifts] = shift[i];
			shift[start + numShifts] = i;
			++numShifts;
		}
	}
}


/**
 * Performs QuickSort on the given arrays, first by @a refIdxArr, then
 * by @a shiftArr, and finally by @a posArr.
 *
 * @param[in,out]	refIdxArr	Reference index array.
 * @param[in,out]	shiftArr	Shift array.
 * @param[in,out]	posArr		Position array.
 * @param			arrSize		Number of elements in each of @a refIdxArr,
 * @a shiftArr, and @a posArr.
 * @param			begin		Array for holding starting indices.
 * @param			end			Array for holding ending indices.
 *
 * @note This code has been adapted from the code written by Darel Rex Finley
 * posted on the site http://alienryderflex.com/quicksort/.
 */
__device__ int quickSort_gpu7(char *refIdxArr, int *shiftArr, int *posArr,
		int arrSize, short *begin, short *end)
{
	int pivotShift, pivotPos, i, left, right;
	char pivotRefIdx;

	begin[0] = 0;
	end[0] = arrSize;
	i = 0;
	while (i >= 0)
	{
		left = begin[i];
		right = end[i] - 1;
		if (left < right)
		{
			pivotRefIdx = refIdxArr[left];
			pivotShift = shiftArr[left];
			pivotPos = posArr[left];

			if (i == arrSize - 1)
				return 0;

			while (left < right)
			{
				while ((refIdxArr[right] > pivotRefIdx
						|| (refIdxArr[right] == pivotRefIdx
								&& shiftArr[right] > pivotShift)
						|| (refIdxArr[right] == pivotRefIdx
								&& shiftArr[right] == pivotShift
								&& posArr[right] >= pivotPos))
						&& left < right)
					--right;

				if (left < right)
				{
					refIdxArr[left] = refIdxArr[right];
					shiftArr[left] = shiftArr[right];
					posArr[left] = posArr[right];
					++left;
				}

				while ((refIdxArr[left] < pivotRefIdx
						|| (refIdxArr[left] == pivotRefIdx
								&& shiftArr[left] < pivotShift)
						|| (refIdxArr[left] == pivotRefIdx
								&& shiftArr[left] == pivotShift
								&& posArr[left] <= pivotPos))
						&& left < right)
					++left;

				if (left < right)
				{
					refIdxArr[right] = refIdxArr[left];
					shiftArr[right] = shiftArr[left];
					posArr[right] = posArr[left];
					--right;
				}
			}
			refIdxArr[left] = pivotRefIdx;
			shiftArr[left] = pivotShift;
			posArr[left] = pivotPos;
			begin[i + 1] = left + 1;
			end[i + 1] = end[i];
			end[i] = left;
			++i;
		}
		else
			--i;
	}
	return 1;
}


/**
 * This is a wrapper function that wraps @a findBiggestClusters_gpu7 function.
 * It has been added so that @a findBiggestClusters_gpu7 can be unit-tested.
 *
 * @param		numClust		Number of clusters.
 * @param		bgstClustSize	Biggest cluster size.
 * @param		clusterSize		Array containing cluster sizes.
 * @param[out]	bgstClust		Array containing biggest clusters.
 */
__global__ void findBiggestClusters_gpu_wrap7(short numClust,
		short bgstClustSize, int *clusterSize, short *bgstClust,
		char *numBgstClust)
{
	*numBgstClust = findBiggestClusters_gpu7(numClust, bgstClustSize,
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
__device__ char findBiggestClusters_gpu7(short numClust, short bgstClustSize,
		int *clusterSize, short *bgstClust)
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
 * This is a wrapper function that wraps @a sort_gpu7 function. It has been
 * added so that @a sort_gpu7 can be unit-tested.
 *
 * @param[in,out]	refPos			Reference positions.
 * @param			arrSize			Size of @a refIdx array.
 */
__global__ void sort_gpu_wrap7(long *refPos, int arrSize)
{
	short threadId = (threadIdx.y * blockDim.x) + threadIdx.x;
	sort_gpu7(refPos, arrSize, threadId);
}


/**
 * Sorts the given arrays in-place using parallel merge-sort. It first sorts
 * using the reference index, then shift, and finally by reference position.
 *
 * @param[in,out]	refPos			Reference positions.
 * @param			arrSize			Size of @a refIdx array.
 * @param			threadId		Current thread ID.
 */
__device__ void sort_gpu7(long *refPos, int arrSize, short threadId)
{
	short numThreads = blockDim.x * blockDim.y * blockDim.z;
	short binSize = (short) ceil(((float) arrSize) / numThreads);
	if (binSize < 2)
	{
		numThreads = numThreads / 2;
		binSize = (short) ceil(((float) arrSize) / numThreads);
	}
	int numBins = (int) ceil(((float) arrSize) / binSize);
	char isSorted;
	short i, begin, end;
	int tmp = 0;
	long tmpLong;
	while (numBins >= 0)
	{
		begin = (threadId * binSize) + (tmp % 2);
		end = min(begin + binSize - 1, arrSize - 1);
		isSorted = FALSE;
		while (isSorted == FALSE)
		{
			isSorted = TRUE;
			for (i = begin + 1; i <= end; ++i)
			{
				if (refPos[i] < refPos[i - 1])
				{
					tmpLong = refPos[i];
					refPos[i] = refPos[i - 1];
					refPos[i - 1] = tmpLong;
					isSorted = FALSE;
				}
			}
		}
		++tmp;
		--numBins;
		__syncthreads();
	}
}


/**
 * This is a wrapper function that wraps @a getHash_gpu7 function. It has been
 * added so that @a getHash_gpu7 can be unit-tested.
 *
 * @param		str		String for which hash is to be calculated.
 * @param		len		Length of @a str.
 * @param[out]	hash	Hash value of @a str.
 */
__global__ void getHash_gpu_wrap7(char *str, int len, int *hash)
{
	*hash = getHash_gpu7(str, len);
}


/**
 * Returns the hash value of the given string.
 *
 * @param	str	String for which hash is to be calculated.
 * @param	len	Length of @a str.
 * @return	Hash value.
 */
__device__ int getHash_gpu7(char *str, int len)
{
	int hash = 0, i;
	for (i = 0; i < len; ++i)
		hash += powerVals_gpu[i] * hashCodes_gpu[str[i]];
	return hash;
}


/**
 * Returns the hash value of the given string.
 *
 * @param	str	String for which hash is to be calculated.
 * @param	len	Length of @a str.
 * @return	Hash value.
 */
int getHash_cpu7(char *str, int len)
{
	int hash = 0, i;
	for (i = 0; i < len; ++i)
		hash += powerVals[i] * hashCodes[str[i]];
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
__device__ int arrSearch_gpu7(int *arr, int arrSize, int num)
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
int arrSearch_cpu7(int *arr, int arrSize, int num)
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
__device__ void arrGetRandomNums_gpu7(int n, int lLimit, int uLimit, int *arr)
{
	int range = uLimit - lLimit + 1;
	int i = 0;
	int randNum;
	while (i < n)
	{
		randNum = (getRandNum_gpu7() % range) + lLimit;
		if (arrSearch_gpu7(arr, i, randNum) == -1)
		{
			arr[i] = randNum;
			++i;
		}
	}
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
void arrGetRandomNums_cpu7(int n, int lLimit, int uLimit, int *arr,
		uint randNumSeed)
{
	int range = uLimit - lLimit + 1;
	int i = 0;
	int randNum;
	while (i < n)
	{
		randNum = (getRandNum_cpu7(randNumSeed) % range) + lLimit;
		if (arrSearch_cpu7(arr, i, randNum) == -1)
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
__device__ int getRandNum_gpu7()
{
	randNumSeed = (randNumSeed * 1103515245) + 12345;
	return ((randNumSeed / 65536) % 32768);
}


/**
 * Returns a pseudo-random number.
 *
 * @return	Random number.
 *
 * @note This algorithm has been taken from "The C Programming Language"
 * by Kernighan and Ritchie.
 */
int getRandNum_cpu7(uint randNumSeed)
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
__device__ int pow_gpu7(int base, int n)
{
	int p = 1;
	while (n > 0)
	{
		p = p * base;
		--n;
	}
	return p;
}

