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
 * Implements GPU-based filtering using seed-and-extend.
 */


#include "common.h"
#include "lookupTable13.h"
#include "preprocess.h"
#include <stdio.h>

#define	REF_POS_BITS2	25
#define	REF_POS_MASK2	33554431
#define	DUMMY_VAL		32

static int *_lookupTable = NULL;
static char *_refIdx = NULL;
static int *_refPos = NULL;
static int *_refPos2 = NULL;
static int _seedLen = 0;
static int _numDistinctTuples = 0;
static int *_numRepeatsPerTuple = NULL;
static int _tupleIgnoreThreshold = 0;
static int _numTotalTuples = 0;
static int *_numActualRepeatsPerTuple = NULL;
static int *_refPosHashMap = NULL;
static long long _maxRefPosComposite = 0;
extern __shared__ int arr_shr[];
__shared__ uint randNumSeed;
__constant__ int hashCodes_gpu[NUM_ASCII_CHARS];
__constant__ int powerVals_gpu[MAX_SEED_LENGTH];


/**
 * Creates a lookup table.
 *
 * @param		refFile					Reference file.
 * @param		seedLen					Seed length.
 * @param		tupleIgnoreThreshold	Tuple that have number of repeats higher
 * than this value will be ignored.
 * @param[out]	totalTuples				Total number of reference tuples that
 * will be used.
 */
void lookupTable13Create(const char *refFile, int seedLen,
		int tupleIgnoreThreshold, int *totalTuples)
{
	_seedLen = seedLen;
	_tupleIgnoreThreshold = tupleIgnoreThreshold;
	_numDistinctTuples = (int) pow((float) DNA_ALPHABET_SIZE, (int) _seedLen);
	_lookupTable = (int *) calloc(_numDistinctTuples, sizeof(int));
	_numActualRepeatsPerTuple = (int *) calloc(_numDistinctTuples, sizeof(int));
	int *numRepeatsPerTuple = (int *) calloc(_numDistinctTuples, sizeof(int));
	_numRepeatsPerTuple = (int *) calloc(_numDistinctTuples, sizeof(int));

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
	_refPosHashMap = (int *) calloc(_maxRefPosComposite, sizeof(int));
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

			refPosComposite = refIdx;
			refPosComposite = refPosComposite << REF_POS_BITS2;
			refPosComposite += (refPos / _seedLen);
			_refPosHashMap[refPosComposite] = hashVal;
		}
		else
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

	int numQryTuples = MAX_QRY_SEQ_LENGTH - seedLen;
	int maxRefTuplesPerQry = numQryTuples * _tupleIgnoreThreshold;
	fprintf(stderr, "Max reference tuples per query = %d\n", maxRefTuplesPerQry);

//	FILE *valsFilePtr = fopen("/home/pgupta/data/gpusw/vals1000.txt", "w");
//	for (i = 0; i <= _numTotalTuples; ++i)
//		fprintf(valsFilePtr, "%d\t%d\t%d\n", _refPos2[i], _refIdx[i], _refPos[i]);
//	fclose(valsFilePtr);
//
//	FILE *keysFilePtr = fopen("/home/pgupta/data/gpusw/keys1000.txt", "w");
//	for (i = 0; i < _numDistinctTuples; ++i)
//		fprintf(keysFilePtr, "%d\t%d\n", _lookupTable[i],
//				_numActualRepeatsPerTuple[i]);
//	fclose(keysFilePtr);

//	FILE *refPosHashMapFilePtr = fopen("/home/pgupta/data/gpusw/refPosHashMap1000.txt", "w");
//	for (i = 0; i < _maxRefPosComposite; ++i)
//		fprintf(refPosHashMapFilePtr, "%d\n", _refPosHashMap[i]);
//	fclose(refPosHashMapFilePtr);
}


/**
 * Creates lookup table using data stored in files.
 *
 * @param keysFile			File containing keys.
 * @param valsFile			File containing values.
 * @param refPosHashMapFile	File containing reference position hash map.
 * @param maxHits			Max number of hits per query.
 * @param seedLen			Seed length.
 * @param tupleIgnoreThres	Threshold value over which any reference tuple
 * will be ignored.
 * @param numTotalTuples	Total number of tuples.
 */
void lookupTable13Create2(const char *keysFile, const char *valsFile,
		const char *refPosHashMapFile, int maxHits, int seedLen,
		int tupleIgnoreThres, int *numTotalTuples)
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
	long long refPosComposite = 0;
	while (fgets(line, 200, valsFilePtr) != NULL)
	{
		sscanf(line, "%d\t%d\t%d", &_refPos2[i], &_refIdx[i], &_refPos[i]);
		refPosComposite = _refIdx[i];
		refPosComposite = refPosComposite << REF_POS_BITS2;
		refPosComposite += (_refPos[i] / seedLen);
		_maxRefPosComposite = max(_maxRefPosComposite, refPosComposite);
		++i;
	}
	fclose(valsFilePtr);

	_refPosHashMap = (int *) calloc(_maxRefPosComposite, sizeof(int));
	FILE *refPosHashMapFilePtr = fopen(refPosHashMapFile, "r");
	i = 0;
	while (fgets(line, 200, refPosHashMapFilePtr) != NULL)
	{
		sscanf(line, "%d", &_refPosHashMap[i]);
		++i;
	}
	fclose(refPosHashMapFilePtr);

	_seedLen = seedLen;
	_tupleIgnoreThreshold = tupleIgnoreThres;
	_numDistinctTuples = (int) pow((float) DNA_ALPHABET_SIZE, (int) _seedLen);
	_numTotalTuples = numVals;
	*numTotalTuples = numVals;

	time(&endTime);
	diffTime = difftime(endTime, startTime);
	fprintf(stderr, "(Time = %.2lf)...", diffTime);
}


/**
 * Copies hash table from CPU to GPU.
 *
 * @param[out]	keys				Keys array on the GPU.
 * @param[out]	numKeys				Number of elements in @a keys array.
 * @param[out]	values				Values array on the GPU.
 * @param[out]	numValues			Number of elements in @a values array.
 * @param[out]	numRepeatsPerTuple	Number of tuples per hash.
 * @param[out]	refPosHashMap		Hash table containing reference positions
 * as the keys and hash values as the values.
 * @param[out]	maxRefPosComposite	Max composite reference position.
 */
void lookupTable13CpyHashTableToGPU(int **keys, int *numKeys, int **values,
		int *numValues, int **numRepeatsPerTuple, int **refPosHashMap,
		int *maxRefPosComposite)
{
	int memKeys = _numDistinctTuples * sizeof(int);
	fprintf(stderr, "   Memory for keys = %d Bytes\n", memKeys);
	cudaMalloc((void **) keys, _numDistinctTuples * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(*keys, _lookupTable, _numDistinctTuples * sizeof(int),
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()
	*numKeys = _numDistinctTuples;

	int memVals = (_numTotalTuples + 1) * sizeof(int);
	fprintf(stderr, "   Memory for values = %d Bytes\n", memVals);
	cudaMalloc((void **) values, (_numTotalTuples + 1) * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(*values, _refPos2, (_numTotalTuples + 1) * sizeof(int),
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()
	*numValues = _numTotalTuples + 1;

	int memNumRptsPerTuple = _numDistinctTuples * sizeof(int);
	fprintf(stderr, "   Memory for numRepeatsPerTuple = %d Bytes\n",
			memNumRptsPerTuple);
	cudaMalloc((void **) numRepeatsPerTuple, _numDistinctTuples * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(*numRepeatsPerTuple, _numActualRepeatsPerTuple,
			_numDistinctTuples * sizeof(int), cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	*maxRefPosComposite = _maxRefPosComposite;
	long memRefPosHashMap = _maxRefPosComposite * sizeof(int);
	fprintf(stderr, "   Memory for memRefPosHashMap = %ld Bytes "
			"(_maxRefPosComposite = %ld)\n", memRefPosHashMap,
			_maxRefPosComposite);

	cudaMalloc((void **) refPosHashMap, _maxRefPosComposite * sizeof(int));
	PRINT_CUDA_ERROR()
	cudaMemcpy(*refPosHashMap, _refPosHashMap, _maxRefPosComposite * sizeof(int),
			cudaMemcpyHostToDevice);
	PRINT_CUDA_ERROR()

	long memTotal = memKeys + memVals + memNumRptsPerTuple + memRefPosHashMap;
	fprintf(stderr, "   Total required GPU memory = %ld Bytes\n", memTotal);
}


/**
 * Maps query sequences to reference sequences on the GPU.
 *
 * @param	keys				Keys array on the GPU.
 * @param	vals				Values array on the GPU.
 * @param	numTuplesPerHash	Array containing number of tuples for each
 * hash.
 * @param	refPosHashMap		Hash table containing reference position as
 * the key and reference tuple hash as value.
 * @param	qrs					Query sequences.
 * @param	qryLen				Length of each query sequence in @qrs array.
 * @param	qryDistance			Distance of each query sequence from the beginning
 * of the @qrs array.
 * @param[out]	refIdx			Output array that will store reference indexes.
 * @param[out]	refPos			Output array that will store reference positions.
 * @param	maxNumHits			Max number of allowed hits.
 * @param	seedLen				Seed length.
 * @param	randNum				Random number.
 * @param	idealArrSize		Size of shared memory in 4 byte blocks.
 */
__global__ void lookupTable13MapQry_gpu(int *keys, int *vals,
		int *numTuplesPerHash, int *refPosHashMap, char *qrs, uchar *qryLen,
		int *qryDistance, int maxQryLen, char *refIdx, int *refPos,
		int maxNumHits, int seedLen, int randNum, int idealArrSize)
{
	int threadId = threadIdx.x;
	int blockId = blockIdx.x;
	int *refPos_shr = (int *) arr_shr;
	int *numMatches_shr = (int *) &refPos_shr[idealArrSize];
	int *refPos2_shr = (int *) &numMatches_shr[idealArrSize];
	int *refPos3_shr = (int *) &refPos2_shr[idealArrSize];
	int *randNum_shr = (int *) &refPos3_shr[idealArrSize];
	int *numMatches3_shr = (int *) &randNum_shr[maxNumHits];
	char *qrySeq_shr = (char *) &numMatches3_shr[idealArrSize];
	__shared__ int maxNumMatches_shr;

	/* Initialize shared arrays. */
	int numIterations2 = (int) ceilf(((float) idealArrSize) / blockDim.x);
	int i, j;
	for (i = 0; i < numIterations2; ++i)
	{
		refPos_shr[(i * blockDim.x) + threadId] = -1;
		numMatches_shr[(i * blockDim.x) + threadId] = -1;
		refPos2_shr[(i * blockDim.x) + threadId] = -1;
		refPos3_shr[(i * blockDim.x) + threadId] = -1;
		numMatches3_shr[(i * blockDim.x) + threadId] = 0;
		qrySeq_shr[(i * blockDim.x) + threadId] = '\0';
	}
	__syncthreads();

	/* Copy the query sequence to shared memory. */
	qrs += blockId * maxQryLen;
	int qrySeqLen = qryLen[blockId];
	int numIterations = (int) ceilf(((float) qrySeqLen) / blockDim.x);
	for (i = 0; i < numIterations; ++i)
		qrySeq_shr[(i * blockDim.x) + threadId] = qrs[(i * blockDim.x) + threadId];
	__syncthreads();
	if (threadId == 0)
	{
		for (i = qrySeqLen; i < (numIterations * blockDim.x); ++i)
			qrySeq_shr[i] = '\0';
	}
	__syncthreads();

	/* Copy each query tuple's hits to shared memory and perform seed and
	 * extend. */
	int increment = -1, qryTupleHash, numHitsPerQryTuple;
	int val, key;
	int refTupleHash, qryTupleHash2, increment2, numMatches, numMismatches;
	int lastTupleIdx = (qrySeqLen / seedLen) * seedLen;
	char s1[MAX_SEED_LENGTH + 1];
	int idx, maxNumMatches = 0;
	__shared__ int count_shr;
	if (threadId == 0)
	{
		count_shr = 0;
		for (i = 0; i <= MAX_SEED_LENGTH; ++i)
			s1[i] = '\0';
	}
	__syncthreads();

	int numTuplesPerQry = qrySeqLen - seedLen + 1;
	while (increment <= numTuplesPerQry)
	{
		/* Copy each query tuple's hits to shared memory. */
		++increment;
		qryTupleHash = getHash13_gpu(qrySeq_shr + increment, seedLen);
		numHitsPerQryTuple = numTuplesPerHash[qryTupleHash];
		numIterations = (int) ceilf(((float) numHitsPerQryTuple) / blockDim.x);
		key = keys[qryTupleHash];
		for (i = 0; i < numIterations; ++i)
			refPos2_shr[(i * blockDim.x) + threadId]
			            = vals[key + (i * blockDim.x) + threadId];
		__syncthreads();
		if (threadId == 0)
		{
			for (i = numHitsPerQryTuple; i < (numIterations * blockDim.x); ++i)
				refPos2_shr[i] = -1;
		}
		__syncthreads();

		/* Seed and extend. */
		for (i = 0; i < numIterations; ++i)
		{
			val = refPos2_shr[(i * blockDim.x) + threadId];
			numMatches = 0;
			numMismatches = 0;
			increment2 = increment;
			if (val < 0)
			{
				numMatches_shr[count_shr + (i * blockDim.x) + threadId] = 0;
				continue;
			}
			while (numMismatches < 5 && increment2 < lastTupleIdx)
			{
				refTupleHash = refPosHashMap[val];
				qryTupleHash2 = getHash13_gpu(qrySeq_shr + increment2, seedLen);
				if (refTupleHash == qryTupleHash2)
					++numMatches;
				else
					++numMismatches;
				++val;
				increment2 += seedLen;
			}
			increment2 -= seedLen;
			if (numMismatches < 5 && (qrySeqLen - increment2) < seedLen
					&& (qrySeqLen - increment2) > 0)
			{
				getUnhash13_gpu(refTupleHash, s1, seedLen);
				refTupleHash = getHash13_gpu(s1, qrySeqLen - increment2);
				if (refTupleHash == qryTupleHash2)
					++numMatches;
			}
			idx = count_shr + (i * blockDim.x) + threadId;
			numMatches_shr[idx] = numMatches;
			refPos_shr[idx] = refPos2_shr[(i * blockDim.x) + threadId];
		}
		__syncthreads();

		int offset = count_shr + numHitsPerQryTuple;
		int gap = (idealArrSize - offset) % blockDim.x;
		if (threadId == 0)
		{
			for (i = offset; i < (offset + gap); ++i)
				numMatches_shr[i] = -1;
		}
		int numIterations3 = (idealArrSize - offset) / blockDim.x;
		for (i = 0; i < numIterations3; ++i)
			numMatches_shr[offset + (i * blockDim.x) + threadId] = -1;
		__syncthreads();

		/* Keep track of the biggest matches. */
		if (threadId == 0)
		{
			int k = count_shr;
			for (j = count_shr; j < offset; ++j)
			{
				if (maxNumMatches <= numMatches_shr[j])
				{
					if (maxNumMatches < numMatches_shr[j])
					{
						maxNumMatches = numMatches_shr[j];
						k = 0;
					}
					numMatches_shr[j] = -1;
					numMatches_shr[k] = maxNumMatches;
					refPos_shr[k] = refPos_shr[j];
					++k;
				}
			}
			if (k > 1000)
				k = 0;
			count_shr = k;
		}
		__syncthreads();
	}

	/* Keep track of the biggest matches. */
	if (threadId == 0)
	{
		maxNumMatches = 0;
		for (i = 0; i <= count_shr; ++i)
			maxNumMatches = (int) fmaxf((float) maxNumMatches,
					(float) numMatches_shr[i]);
		maxNumMatches_shr = maxNumMatches;
	}
	__syncthreads();

	for (i = 0; i < numIterations2; ++i)
	{
		if (numMatches_shr[(i * blockDim.x) + threadId] == maxNumMatches_shr)
			++numMatches3_shr[threadId];
	}
	__syncthreads();
	if (threadId == 0)
	{
		refPos2_shr[0] = 0;
		for (i = 1; i < blockDim.x; ++i)
			refPos2_shr[i] = refPos2_shr[i - 1] + numMatches3_shr[i - 1];
	}
	__syncthreads();
	int tmp = refPos2_shr[threadId];
	for (i = 0; i < numIterations2; ++i)
	{
		if (numMatches_shr[(i * blockDim.x) + threadId] == maxNumMatches_shr)
		{
			refPos3_shr[tmp] = refPos_shr[(i * blockDim.x) + threadId];
			++tmp;
		}
	}
	__syncthreads();

	/* Store the best matches in the global memory. */
	if (threadId == 0)
	{
		int numHits = refPos2_shr[blockDim.x - 1] + numMatches3_shr[blockDim.x - 1];
		if (numHits > maxNumHits)
		{
			arrGetRandomNums13_gpu(maxNumHits, 0, numHits - 1, randNum_shr);
			for (i = 0; i < maxNumHits; ++i)
			{
				refIdx[(blockId * maxNumHits) + i] =
						(char) (refPos3_shr[randNum_shr[i]] >> REF_POS_BITS2);
				refPos[(blockId * maxNumHits) + i] =
						(refPos3_shr[randNum_shr[i]] & REF_POS_MASK2) * seedLen;
			}
		}
		else
		{
			for (i = 0; i < maxNumHits; ++i)
			{
				refIdx[(blockId * maxNumHits) + i] =
						(char) (refPos3_shr[i] >> REF_POS_BITS2);
				refPos[(blockId * maxNumHits) + i] =
						(refPos3_shr[i] & REF_POS_MASK2) * seedLen;
			}
		}
	}
}


/**
 * Maps query sequences to reference sequences on the GPU.
 *
 * @param	keys				Keys array on the GPU.
 * @param	vals				Values array on the GPU.
 * @param	numTuplesPerHash	Array containing number of tuples for each
 * hash.
 * @param	refPosHashMap		Hash table containing reference position as
 * the key and reference tuple hash as value.
 * @param	qrs					Query sequences.
 * @param	qryLen				Length of each query sequence in @qrs array.
 * @param	qryDistance			Distance of each query sequence from the beginning
 * of the @qrs array.
 * @param[out]	refIdx			Output array that will store reference indexes.
 * @param[out]	refPos			Output array that will store reference positions.
 * @param	maxNumHits			Max number of allowed hits.
 * @param	seedLen				Seed length.
 * @param	randNum				Random number.
 * @param	idealArrSize		Size of shared memory in 4 byte blocks.
 */
__global__ void lookupTable13MapQry2_gpu(int *keys, int *vals,
		int *numTuplesPerHash, int *refPosHashMap, char *qrs, uchar *qryLen,
		int *qryDistance, int maxQryLen, char *refIdx, int *refPos,
		int maxNumHits, int seedLen, int randNum, int idealArrSize)
{
	int threadId = threadIdx.x;
	int blockId = blockIdx.x;
	int *refPos_shr = (int *) arr_shr;
	int *numMatches_shr = (int *) &refPos_shr[idealArrSize];
	int *refPos2_shr = (int *) &numMatches_shr[idealArrSize];
	int *refPos3_shr = (int *) &refPos2_shr[idealArrSize];
	int *randNum_shr = (int *) &refPos3_shr[idealArrSize];
	int *numMatches3_shr = (int *) &randNum_shr[maxNumHits];
	char *qrySeq_shr = (char *) &numMatches3_shr[idealArrSize];
	__shared__ int maxNumMatches_shr;

	/* Initialize shared arrays. */
	int numIterations2 = (int) ceilf(((float) idealArrSize) / blockDim.x);
	int i, j;
	for (i = 0; i < numIterations2; ++i)
	{
		refPos_shr[(i * blockDim.x) + threadId] = -1;
		numMatches_shr[(i * blockDim.x) + threadId] = -1;
		refPos2_shr[(i * blockDim.x) + threadId] = -1;
		refPos3_shr[(i * blockDim.x) + threadId] = -1;
		numMatches3_shr[(i * blockDim.x) + threadId] = 0;
		qrySeq_shr[(i * blockDim.x) + threadId] = '\0';
	}
	__syncthreads();

	/* Copy the query sequence to shared memory. */
	qrs += blockId * maxQryLen;
	int qrySeqLen = qryLen[blockId];
	int numIterations = (int) ceilf(((float) qrySeqLen) / blockDim.x);
	for (i = 0; i < numIterations; ++i)
		qrySeq_shr[(i * blockDim.x) + threadId] = qrs[(i * blockDim.x) + threadId];
	__syncthreads();
	if (threadId == 0)
	{
		for (i = qrySeqLen; i < (numIterations * blockDim.x); ++i)
			qrySeq_shr[i] = '\0';
	}
	__syncthreads();

	/* Copy each query tuple's hits to shared memory and perform seed and
	 * extend. */
	int increment = -1, qryTupleHash, numHitsPerQryTuple;
	int val, key;
	int refTupleHash, increment2, numMatches, numMismatches;
	int lastTupleIdx = (qrySeqLen / seedLen) * seedLen;
	int idx, maxNumMatches = 0, numMismatches2;
	__shared__ int count_shr;
	if (threadId == 0)
		count_shr = 0;
	__syncthreads();

	int numTuplesPerQry = qrySeqLen - seedLen + 1;
	while (increment <= numTuplesPerQry)
	{
		/* Copy each query tuple's hits to shared memory. */
		++increment;
		qryTupleHash = getHash13_gpu(qrySeq_shr + increment, seedLen);
		numHitsPerQryTuple = numTuplesPerHash[qryTupleHash];
		numIterations = (int) ceilf(((float) numHitsPerQryTuple) / blockDim.x);
		key = keys[qryTupleHash];
		for (i = 0; i < numIterations; ++i)
			refPos2_shr[(i * blockDim.x) + threadId]
			            = vals[key + (i * blockDim.x) + threadId];
		__syncthreads();
		if (threadId == 0)
		{
			for (i = numHitsPerQryTuple; i < (numIterations * blockDim.x); ++i)
				refPos2_shr[i] = -1;
		}
		__syncthreads();

		/* Seed and extend. */
		for (i = 0; i < numIterations; ++i)
		{
			val = refPos2_shr[(i * blockDim.x) + threadId];
			numMatches = 0;
			numMismatches = 0;
			numMismatches2 = 0;
			increment2 = increment;
			if (val < 0)
			{
				numMatches_shr[count_shr + (i * blockDim.x) + threadId] = 0;
				continue;
			}
			refTupleHash = refPosHashMap[val];
			while (increment2 < lastTupleIdx)
			{
				if (refTupleHash == getHash13_gpu(qrySeq_shr + increment2, seedLen))
				{
					++val;
					refTupleHash = refPosHashMap[val];
					increment2 += seedLen;
					++numMatches;
					numMismatches = 0;
					numMismatches2 = 0;
				}
				else if (numMismatches < (seedLen - 1))
				{
					++increment2;
					++numMismatches;
				}
				else if (numMismatches2 == 0)
				{
					++val;
					refTupleHash = refPosHashMap[val];
					increment2 -= numMismatches;
					numMismatches = 0;
					++numMismatches2;
				}
				else
					break;
			}
			idx = count_shr + (i * blockDim.x) + threadId;
			numMatches_shr[idx] = numMatches;
			refPos_shr[idx] = refPos2_shr[(i * blockDim.x) + threadId];
		}
		__syncthreads();

		int offset = count_shr + numHitsPerQryTuple;
		int gap = (idealArrSize - offset) % blockDim.x;
		if (threadId == 0)
		{
			for (i = offset; i < (offset + gap); ++i)
				numMatches_shr[i] = -1;
		}
		int numIterations3 = (idealArrSize - offset) / blockDim.x;
		for (i = 0; i < numIterations3; ++i)
			numMatches_shr[offset + (i * blockDim.x) + threadId] = -1;
		__syncthreads();

		/* Keep track of the biggest matches. */
		if (threadId == 0)
		{
			int k = count_shr;
			for (j = count_shr; j < offset; ++j)
			{
				if (maxNumMatches <= numMatches_shr[j])
				{
					if (maxNumMatches < numMatches_shr[j])
					{
						maxNumMatches = numMatches_shr[j];
						k = 0;
					}
					numMatches_shr[j] = -1;
					numMatches_shr[k] = maxNumMatches;
					refPos_shr[k] = refPos_shr[j];
					++k;
				}
			}
			if (k > 1000)
				k = 0;
			count_shr = k;
		}
		__syncthreads();
	}

	/* Keep track of the biggest matches. */
	if (threadId == 0)
	{
		maxNumMatches = 0;
		for (i = 0; i <= count_shr; ++i)
			maxNumMatches = (int) fmaxf((float) maxNumMatches,
					(float) numMatches_shr[i]);
		maxNumMatches_shr = maxNumMatches;
	}
	__syncthreads();

	for (i = 0; i < numIterations2; ++i)
	{
		if (numMatches_shr[(i * blockDim.x) + threadId] == maxNumMatches_shr)
			++numMatches3_shr[threadId];
	}
	__syncthreads();
	if (threadId == 0)
	{
		refPos2_shr[0] = 0;
		for (i = 1; i < blockDim.x; ++i)
			refPos2_shr[i] = refPos2_shr[i - 1] + numMatches3_shr[i - 1];
	}
	__syncthreads();
	int tmp = refPos2_shr[threadId];
	for (i = 0; i < numIterations2; ++i)
	{
		if (numMatches_shr[(i * blockDim.x) + threadId] == maxNumMatches_shr)
		{
			refPos3_shr[tmp] = refPos_shr[(i * blockDim.x) + threadId];
			++tmp;
		}
	}
	__syncthreads();

	/* Store the best matches in the global memory. */
	if (threadId == 0)
	{
		int numHits = refPos2_shr[blockDim.x - 1] + numMatches3_shr[blockDim.x - 1];
		if (numHits > maxNumHits)
		{
			arrGetRandomNums13_gpu(maxNumHits, 0, numHits - 1, randNum_shr);
			for (i = 0; i < maxNumHits; ++i)
			{
				refIdx[(blockId * maxNumHits) + i] =
						(char) (refPos3_shr[randNum_shr[i]] >> REF_POS_BITS2);
				refPos[(blockId * maxNumHits) + i] =
						(refPos3_shr[randNum_shr[i]] & REF_POS_MASK2) * seedLen;
			}
		}
		else
		{
			for (i = 0; i < maxNumHits; ++i)
			{
				refIdx[(blockId * maxNumHits) + i] =
						(char) (refPos3_shr[i] >> REF_POS_BITS2);
				refPos[(blockId * maxNumHits) + i] =
						(refPos3_shr[i] & REF_POS_MASK2) * seedLen;
			}
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
__device__ void arrGetRandomNums13_gpu(int n, int lLimit, int uLimit, int *arr)
{
	int range = uLimit - lLimit + 1;
	int i = 0;
	int randNum;
	while (i < n)
	{
		randNum = (getRandNum13_gpu() % range) + lLimit;
		if (arrSearch13_gpu(arr, i, randNum) == -1)
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
__device__ int getRandNum13_gpu()
{
	randNumSeed = (randNumSeed * 1103515245) + 12345;
	return ((randNumSeed / 65536) % 32768);
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
__device__ int arrSearch13_gpu(int *arr, int arrSize, int num)
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
 * Returns the hash value of the given string.
 *
 * @param	str	String for which hash is to be calculated.
 * @param	len	Length of @a str.
 * @return	Hash value.
 */
__device__ int getHash13_gpu(char *str, int len)
{
	int hash = 0, i;
	for (i = 0; i < len; ++i)
		hash += powerVals_gpu[i] * hashCodes_gpu[str[i]];
	return hash;
}


/**
 * Decodes the given hash.
 *
 * @param		hash	Hash value to be decoded.
 * @param[out]	s		Decoded sequence.
 * @param		seedLen	Seed length.
 */
__device__ void getUnhash13_gpu(int hash, char *s, int seedLen)
{
	int remainder, i = -1;
	while (hash > 0)
	{
		remainder = hash % 4;
		++i;
		if (remainder == 0)
			s[i] = 'A';
		else if (remainder == 1)
			s[i] = 'C';
		else if (remainder == 2)
			s[i] = 'G';
		else if (remainder == 3)
			s[i] = 'T';
		else
			s[i] = 'N';
		hash = hash / 4;
	}

	/* Pad with 'A's */
	while (i < (seedLen - 1))
	{
		++i;
		s[i] = 'A';
	}
	s[seedLen] = '\0';
}


/**
 * Returns the reference index and reference position by breaking down the
 * given composite value.
 *
 * @param		compositeVal	Composite value that is to be broken down.
 * @param[out]	refIdx			Reference index.
 * @param[out]	refPos			Reference position.
 * @param		seedLen			Seed length.
 */
__device__ void getRefIdxAndPos13_gpu(int compositeVal, int *refIdx, int *refPos,
		int seedLen)
{
	*refIdx = compositeVal >> REF_POS_BITS2;
	*refPos = (compositeVal & REF_POS_MASK2) * seedLen;
}


/**
 * Copies constant memory from CPU to GPU.
 */
void lookupTable13CpyConstMemToGPU()
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
