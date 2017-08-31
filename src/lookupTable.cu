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
#include "lookupTable.h"
#include "preprocess.h"
#include <stdio.h>
#include "array.h"


static RefPosList **_lookupTable = NULL;
static ushort *_tupleCount = NULL;
static int _numTuples = 0;
static int _seedLen = 0;
static int _maxNumHits = 0;


/**
 * Creates a lookup table.
 *
 * @param 	refFile		Reference file.
 * @param	seedLength	Seed length.
 * @param	maxNumHits	Max number of hits.
 */
void lookupTableCreate(const char *refFile, int seedLength, int maxNumHits)
{
	if (refFile == NULL)
	{
		fprintf(stderr, "Error: The given reference file is NULL. Exiting.\n");
		exit(EXIT_FAILURE);
	}
	else if (seedLength == 0)
	{
		fprintf(stderr, "Error: Seed length is 0. Exiting.\n");
		exit(EXIT_FAILURE);
	}

	_seedLen = seedLength;
	_maxNumHits = maxNumHits;
	_numTuples = (int) pow((float) DNA_ALPHABET_SIZE, (int) _seedLen);
	_lookupTable = (RefPosList **) malloc(_numTuples * sizeof(RefPosList *));
	RefPosList **tail = (RefPosList **) malloc(_numTuples * sizeof(RefPosList *));
	_tupleCount = (ushort *) malloc(_numTuples * sizeof(ushort));
	RefPosList *tmpRefPosList = NULL;
	int numIterations, i, j, offset = 0, numBases = 0, hashVal, lineLength;
	char refIndex = -1, line[MAX_LINE_LENGTH];
	for (i = 0; i < _numTuples; ++i)
	{
		_lookupTable[i] = NULL;
		tail[i] = NULL;
		_tupleCount[i] = 0;
	}

	/* Parse reference sequence file. */
	fprintf(stderr, "   Processing ref: ");
	FILE *filePtr = fopen(refFile, "r");
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
		{
			offset = 0;
			continue;
		}
		/* This line contains reference ID. */
		else if (line[offset] == '>')
		{
			offset = 0;
			numBases = 0;
			++refIndex;
			fprintf(stderr, "%d; ", refIndex);
//			if (refIndex > 7)
//				break;
		}
		/* This is a line containing sequence. */
		else
		{
			numIterations = lineLength - _seedLen + 1;
			/* Consider non-overlapping tuples only. */
			for (i = 0; i < numIterations; i = i + _seedLen)
			{
				hashVal = getHash(line + i, _seedLen);
				if (offset > 0)
					offset = 0;
				if(_tupleCount[hashVal] > TUPLE_IGNORE_THRES)
					continue;
				++_tupleCount[hashVal];
				if(_tupleCount[hashVal] > TUPLE_IGNORE_THRES)
				{
					refPosListDelete(_lookupTable[hashVal]);
					_lookupTable[hashVal] = NULL;
				}
				else
				{
					tmpRefPosList = refPosListCreate(refIndex, numBases + i);
					_lookupTable[hashVal] = refPosListPushBack(
							_lookupTable[hashVal], tail[hashVal], tmpRefPosList);
					tail[hashVal] = tmpRefPosList;
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
	fprintf(stderr, "\n");
	fclose(filePtr);
	free(tail);
}


/**
 * Frees memory occupied by lookup table.
 */
void lookupTableDelete()
{
	free(_tupleCount);
	_tupleCount = NULL;
	if (_lookupTable != NULL)
	{
		int i;
		for (i = 0; i < _numTuples; ++i)
			refPosListDelete(_lookupTable[i]);
		free(_lookupTable);
		_lookupTable = NULL;
	}
	_numTuples = 0;
	_seedLen = 0;
	_maxNumHits = 0;
}


/**
 * Searches the query sequence in the reference and returns
 * the best matching reference coordinates.
 *
 * @param hashTable Hash table containing hash as key and a list of reference
 * positions as values.
 * @param query Query object.
 * @param seqLength Sequence length.
 * @param isReverseComplement A value of '1' indicates that the
 * reverse complement of the query should be used for search; '0' indicates
 * otherwise.
 * @param seedLength Seed length to be used.
 * @param maxNumHits Max number of hits to be returned.
 * @return List of best-matching reference coordinates.
 */
HitList *lookupTableMapQry(const Query *query, uint qryLen, uint isRevComp)
{
	static int numQryTuples, i, hash;
	static RefPosList *tmpRefPosList;
	static HitList *hitList, *tmpHitList, *tail;
	static char *seq;
	char seed[_seedLen + 1];

	if (isRevComp == 1)
		seq = query->revCompSeq;
	else
		seq = query->seq;


	/* Step 1: Break the query into tuples. */
	hitList = NULL;
	tail = NULL;
	tmpRefPosList = NULL;
	numQryTuples = qryLen - _seedLen + 1;
	for (i = 0; i < numQryTuples; ++i)
	{
		/* Step 2: Calculate the hash for each tuple. */
		/* Step 3: Find the reference coordinates of the hash in the hash
		 * table. */
		/* Step 4: Aggregate the reference coordinate positions. */
		strncpy(seed, seq + i, _seedLen);
		hash = getHash(seed, _seedLen);
		tmpRefPosList = _lookupTable[hash];

		/* Add the tuples to the list. */
		while (tmpRefPosList != NULL)
		{
			tmpHitList = hitListCreateNode(tmpRefPosList->index,
					tmpRefPosList->position - i, tmpRefPosList->position);
			hitList = hitListPushBack(hitList, tail, tmpHitList);
			tail = tmpHitList;
			tmpRefPosList = tmpRefPosList->next;
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
	/* If the list is empty, return NULL. */
	if (list == NULL)
		return NULL;
	/* If the list has only one node, return a copy of that node. */
	else if (list->next == NULL)
	{
		HitList *tmpHitList = hitListDuplicateNode(list);
		return tmpHitList;
	}
	/* If the list has more than one node, follow these steps below. */
	else
	{
		/* Step 1: Find number of clusters that are to be created. */
		HitList *tmpHitList = list;
		int index = -1, shift = -1, numClusters = 0;
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
		HitList **bestMatches = (HitList **)
				malloc(numClusters * sizeof(HitList *));
		if (bestMatches == NULL)
		{
			fprintf(stderr, "Error allocating enough memory in %s at "
					"line %d. Exiting.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		HitList *tmpBestMatches = NULL;
		int i = -1, j = 0;
		index = -1;
		shift = -1;
		tmpHitList = list;
		for (j = 0; j < numClusters; ++j)
			bestMatches[j] = NULL;
		HitList *tail = NULL;
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
		int maxClusterSize = 0, tmpSize;
		int *clusterSizeArr = (int *) malloc(numClusters * sizeof(int));
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
		int numMaxSizeClusters = 0;
		for (j = 0; j < numClusters; ++j)
		{
			if (clusterSizeArr[j] == maxClusterSize)
				++numMaxSizeClusters;
		}


		/* Step 5: Put first node of each cluster that has a max size
		 * in a separate list if its index matches one of the random numbers. */
		uint numRandomNumbers = _maxNumHits;
		if (numRandomNumbers > numMaxSizeClusters)
			numRandomNumbers = numMaxSizeClusters;
		uint randomNumArr[numRandomNumbers];
		arrGetRandomNums(numRandomNumbers, 1, numMaxSizeClusters, randomNumArr);
		qsort(randomNumArr, numRandomNumbers, sizeof(uint), compare);

		HitList *bestMatches2 = NULL, *tmpBestMatches2 = NULL;
		tail = NULL;
		int l = 0;
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
