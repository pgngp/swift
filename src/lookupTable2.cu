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
#include "lookupTable2.h"
#include "preprocess.h"
#include <stdio.h>
#include "array.h"
#include <limits.h>
#include "refMap.h"
#include <time.h>


static int *_lookupTable = NULL;
static uint *_refPos = NULL;
static int _seedLen = 0;
static int _maxHitsPerQry = 0;
static int _numDistinctTuples = 0;
static int *_numRepeatsPerTuple = NULL;
static char *_seed = NULL;
static int _tupleIgnoreThreshold = 0;


/**
 * This is a wrapper function that wraps @a lookupTable2Create function.
 * This function has been added so that @a lookupTable2Create can be
 * unit-tested.
 *
 * @param		refFile					Reference file.
 * @param		seedLen					Seed length.
 * @param 		maxHitsPerQry			Maximum hits per query.
 * @param[out]	lookupTable				Lookup table.
 * @param[out]	refPos					Reference sequence offset.
 * @param[out]	numDistinctTuples		Total number of distinct tuples.
 * @param[out]	numRepeatsPerTuple		Number of repeats per tuple.
 * @param		tupleIgnoreThreshold	Tuple that have number of repeats higher
 * than this value will be ignored.
 */
void lookupTable2Create_wrap(const char *refFile, int seedLen,
		int maxHitsPerQry, int **lookupTable, uint **refPos,
		int *numDistinctTuples, int **numRepeatsPerTuple,
		int tupleIgnoreThreshold)
{
	lookupTable2Create(refFile, seedLen, maxHitsPerQry, tupleIgnoreThreshold);
	*lookupTable = _lookupTable;
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
void lookupTable2Create(const char *refFile, int seedLen, int maxHitsPerQry,
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
	int hashVal, lineLength, numBases = 0, numNewLines = 0;
	char line[MAX_LINE_LENGTH];

	char tmpRefFile[MAX_FILE_NAME_LENGTH];
	sprintf(tmpRefFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR, TEMP_REF_FILE);
	FILE *tmpRefFilePtr = fopen(tmpRefFile, "w+");
	uint refPos = 0;
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
			numNewLines = 0;
			refPos = ftell(filePtr);
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
					if (offset > 0)
					{
						fprintf(tmpRefFilePtr, "%d\t%u\n", hashVal,
								refPos + numBases + i + numNewLines - 1);
						offset = 0;
					}
					else
					{
						fprintf(tmpRefFilePtr, "%d\t%u\n", hashVal,
								refPos + numBases + i + numNewLines);
					}
				}
				else if (_numRepeatsPerTuple[hashVal]
				                            == numRepeatsPerTuple[hashVal] + 1)
				{
					numRepeatsPerTuple[hashVal] = 0;
					numTotalTuples -= _tupleIgnoreThreshold;
				}
			}
			numBases += lineLength;
			++numNewLines;

			/* Copy the last few bases to the beginning of 'line' array */
			offset = lineLength - i;
			for (j = 0; j < offset; ++j)
				line[j] = line[i + j];
			numBases -= offset;
		}
	}
	_refPos = (uint *) calloc((numTotalTuples + 1), sizeof(uint));
	_refPos[0] = UINT_MAX; /* First element set to max value so that
	tuples that do not exist in the reference or that have number of repeats
	greater than a threshold value can point to this element. */
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
	int index;
	while (fgets(line, MAX_LINE_LENGTH, tmpRefFilePtr) != NULL)
	{
		sscanf(line, "%d\t%u", &hashVal, &refPos);
		if (_numRepeatsPerTuple[hashVal] <= _tupleIgnoreThreshold)
		{
			index = _lookupTable[hashVal] + numRepeatsPerTuple[hashVal] - 1;
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
}


/**
 * Releases memory occupied by data structures in this file.
 */
void lookupTable2Delete()
{
	free(_lookupTable);
	_lookupTable = NULL;
	free(_refPos);
	_refPos = NULL;
	_seedLen = 0;
	_maxHitsPerQry = 0;
	_numDistinctTuples = 0;
	free(_seed);
	_seed = NULL;
	free(_numRepeatsPerTuple);
	_numRepeatsPerTuple = NULL;
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
 * @return 	List of best-matching reference coordinates.
 */
HitList *lookupTable2MapQry(const Query *query, uint qryLen, uint isRevComp)
{
	static int numQryTuples, i, j, hash, startIdx, endIdx, refIdx, refPos;
	static HitList *hitList, *tmpHitList, *tail;
	static char *seq;

	if (isRevComp == 1)
		seq = query->revCompSeq;
	else
		seq = query->seq;


	/* Step 1: Break the query into tuples. */
	hitList = NULL;
	tail = NULL;
	numQryTuples = qryLen - _seedLen + 1;
	for (i = 0; i < numQryTuples; ++i)
	{
		/* Step 2: Calculate the hash for each tuple. */
		/* Step 3: Find the reference coordinates of the hash in the hash
		 * table. */
		/* Step 4: Aggregate the reference coordinate positions. */
		strncpy(_seed, seq + i, _seedLen);
		hash = getHash(_seed, _seedLen);
		if (_numRepeatsPerTuple[hash] > _tupleIgnoreThreshold)
			continue;

		/* Add the tuples to the list. */
		startIdx = _lookupTable[hash];
		endIdx = startIdx + _numRepeatsPerTuple[hash] - 1;
		for (j = startIdx; j <= endIdx; ++j)
		{
			refMapGetIndexAndPos(_refPos[j], &refIdx, &refPos);
			tmpHitList = hitListCreateNode((char) refIdx, refPos - i, refPos);
			hitList = hitListPushBack(hitList, tail, tmpHitList);
			tail = tmpHitList;
		}
	}

	/* Step 5: Sort the coordinate positions. */
	hitList = hitListSort(hitList, hitListGetSize(hitList));

	/* Step 6: Find the best-matching positions. */
	tmpHitList = getBestScoringPosition(hitList);

	/* Step 7: Delete the HitList. */
	hitListDelete(hitList);

	/* Step 8: Return the best-matching positions. */
	return tmpHitList;
}


/**
 * Returns the positions of the best-matching reference sequences.
 *
 * @param list The HitList list from which best-matching positions are to
 * be found.
 * @param maxNumHits Max number of hits to be returned.
 * @return List of best-matching reference coordinates.
 */
static HitList *getBestScoringPosition(HitList *list)
{
	static int index, shift, numClusters, i, j, l, maxClusterSize, tmpSize;
	static int numMaxSizeClusters, *clusterSizeArr;
	static uint numRandomNumbers;
	static HitList *bestMatches2, *tmpBestMatches, *tmpBestMatches2, *tail;
	static HitList **bestMatches, *tmpHitList;

	/* If the list is empty, return NULL. */
	if (list == NULL)
		return NULL;
	/* If the list has only one node, return a copy of that node. */
	else if (list->next == NULL)
		return hitListDuplicateNode(list);
	/* If the list has more than one node, follow these steps below. */
	else
	{
		/* Step 1: Find number of clusters that are to be created. */
		tmpHitList = list;
		index = -1, shift = -1, numClusters = 0;
		while (tmpHitList->next != NULL)
		{
			if (tmpHitList->index != index || tmpHitList->shift != shift)
			{
				++numClusters;
				index = tmpHitList->index;
				shift = tmpHitList->shift;
			}
			tmpHitList = tmpHitList->next;
		}
		if (tmpHitList->index != index || tmpHitList->shift != shift)
			++numClusters;


		/* Step 2: Create clusters. */
		bestMatches = (HitList **) malloc(numClusters * sizeof(HitList *));
		if (bestMatches == NULL)
		{
			fprintf(stderr, "Error allocating enough memory in %s at "
					"line %d. Exiting.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmpBestMatches = NULL;
		i = -1;
		index = -1;
		shift = -1;
		tmpHitList = list;
		for (j = 0; j < numClusters; ++j)
			bestMatches[j] = NULL;
		tail = NULL;
		while (tmpHitList->next != NULL)
		{
			if (tmpHitList->index != index || tmpHitList->shift != shift)
			{
				++i;
				index = tmpHitList->index;
				shift = tmpHitList->shift;
				tail = NULL;
			}
			tmpBestMatches = hitListDuplicateNode(tmpHitList);
			bestMatches[i] = hitListPushBack(bestMatches[i], tail,
					tmpBestMatches);
			tail = tmpBestMatches;
			tmpHitList = tmpHitList->next;
		}
		if (tmpHitList->index != index || tmpHitList->shift != shift)
		{
			++i;
			tail = NULL;
		}
		tmpBestMatches = hitListDuplicateNode(tmpHitList);
		bestMatches[i] = hitListPushBack(bestMatches[i], tail, tmpBestMatches);


		/* Step 3: (1) Find max cluster size. (2) Store the size of each
		 * cluster in an array. */
		maxClusterSize = 0;
		clusterSizeArr = (int *) malloc(numClusters * sizeof(int));
		if (clusterSizeArr == NULL)
		{
			fprintf(stderr, "Error allocating emough memory in %s at "
					"line %d. Exiting.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		for (j = 0; j < numClusters; ++j)
		{
			tmpSize = hitListGetSize(bestMatches[j]);
			clusterSizeArr[j] = tmpSize;
			maxClusterSize = max(maxClusterSize, tmpSize);
		}


		/* Step 4: Find number of max size clusters */
		numMaxSizeClusters = 0;
		for (j = 0; j < numClusters; ++j)
		{
			if (clusterSizeArr[j] == maxClusterSize)
				++numMaxSizeClusters;
		}


		/* Step 5: Put first node of each cluster that has a max size
		 * in a separate list if its index matches one of the random numbers. */
		numRandomNumbers = _maxHitsPerQry;
		if (numRandomNumbers > numMaxSizeClusters)
			numRandomNumbers = numMaxSizeClusters;
		uint randomNumArr[numRandomNumbers];
		arrGetRandomNums(numRandomNumbers, 1, numMaxSizeClusters, randomNumArr);
		qsort(randomNumArr, numRandomNumbers, sizeof(uint), compare);

		bestMatches2 = NULL, tmpBestMatches2 = NULL;
		tail = NULL;
		l = 0;
		for (j = 0; j < numClusters; ++j)
		{
			if (clusterSizeArr[j] == maxClusterSize)
			{
				++l;
				if (arrSearch((int *) randomNumArr, numRandomNumbers, l) > -1)
				{
					tmpBestMatches2 = hitListDuplicateNode(bestMatches[j]);
					bestMatches2 = hitListPushBack(bestMatches2, tail,
							tmpBestMatches2);
					tail = tmpBestMatches2;
				}
			}
		}


		/* Step 6: Delete the clusters. */
		for (i = 0; i < numClusters; ++i)
			hitListDelete(bestMatches[i]);
		free(bestMatches);
		free(clusterSizeArr);


		/* Step 7: Return list of best-matching reference coordinates. */
		return bestMatches2;
	}
}


/**
 * Returns the integer difference between the given two parameters.
 *
 * @param a First parameter.
 * @param b Second parameter.
 * @return The difference between the given two parameters.
 */
static int compare(const void *a, const void *b)
{
	return (*(int *)a - *(int *)b);
}

