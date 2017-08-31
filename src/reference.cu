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

#include <stdio.h>
#include "reference.h"
#include "common.h"
#include <math.h>
#include "preprocess.h"
#include "array.h"
#include <assert.h>
#include <limits.h>

#define MAX_HITS			1100	/* This value should be close to
(TUPLE_IGNORE_THRES * Number of tuples per query). */
#define	MAX_CLUSTER_SIZE	900		/* This value should be less than or equal
to MAX_HITS. */
#define	MAX_RAND_NUMS		20 		/* This value should be less than or equal
to MAX_CLUSTER_SIZE. */

__shared__ uint randNumSeed;

/**
 * Creates a hash table from the reference sequences and returns the hash
 * table.
 *
 * @param refFile Reference file path. This file must be in FASTA format.
 * @param seedLength Seed length to be used. This value must be greater than 0.
 * Seed length must be greater than length of each reference.
 * @return Hash table containing the hash of seed as keys and list of reference
 * positions as values.
 */
RefPosList **refPreprocess(const char *refFile, uint seedLength)
{
	if (refFile == NULL)
	{
		fprintf(stderr, "Error: The given reference file is NULL. Exiting.\n");
		exit(EXIT_FAILURE);
	}
	else if (seedLength == 0)
		return NULL;


	/*
	 * Create the hash table.
	 */
	int numRefTuples = (int) pow((float) DNA_ALPHABET_SIZE, (int) seedLength);
	RefPosList **hashTable =
			(RefPosList **) malloc(numRefTuples * sizeof(RefPosList *));
	RefPosList **tail = (RefPosList **)
			malloc(numRefTuples * sizeof(RefPosList *));
	ushort *tupleCount = (ushort *) malloc(numRefTuples * sizeof(ushort));
	RefPosList *tmpRefPosList = NULL;
	int numIterations, i, j, offset = 0, numNewLines = 0;
	int numBases = 0;
	int hashVal;
	char refIndex = -1;
	ulong nameOffset = 0, seqOffset = 0, fragOffset = 0;
	FILE *filePtr = fopen(refFile, "r");
	char line[MAX_LINE_LENGTH];
	int lineLength;
	for (i = 0; i < numRefTuples; ++i)
	{
		hashTable[i] = NULL;
		tail[i] = NULL;
		tupleCount[i] = 0;
	}

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
			nameOffset = ftell(filePtr) - lineLength + offset - 1;
			seqOffset = ftell(filePtr);
			offset = 0;
			numBases = 0;
			numNewLines = 0;
			++refIndex;
		}
		/* This is a line containing sequence. */
		else
		{
			numIterations = lineLength - seedLength + 1;
			/* Consider non-overlapping tuples only. */
			for (i = 0; i < numIterations; i = i + seedLength)
			{
				hashVal = getHash(line + i, seedLength);
				if(tupleCount[hashVal] > TUPLE_IGNORE_THRES)
					continue;
				if (offset > 0)
				{
					fragOffset = seqOffset + numBases + i + numNewLines - 1;
					offset = 0;
				}
				else
					fragOffset = seqOffset + numBases + i + numNewLines;
				tmpRefPosList = refPosListCreate(refIndex, numBases + i);
				hashTable[hashVal] = refPosListPushBack(hashTable[hashVal],
						tail[hashVal], tmpRefPosList);
				tail[hashVal] = tmpRefPosList;
				++tupleCount[hashVal];
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
	fclose(filePtr);
	free(tail);
	free(tupleCount);

	return hashTable;
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
HitList *refSearchQuery(RefPosList **hashTable,
		const Query *query, uint seqLength, uint isReverseComplement,
		uint seedLength, uint maxNumHits)
{
	static int numQryTuples;
	static int i, count, isCountHigh;
	static RefPosList *tmpRefPosList, *tmpRefPosList2;
	static HitList *hitList, *tmpHitList, *tail;
	static char *seq;
	char seed[seedLength + 1];

	if (isReverseComplement == 1)
		seq = query->revCompSeq;
	else
		seq = query->seq;


	/* Step 1: Break the query into tuples. */
	hitList = NULL;
	tail = NULL;
	tmpRefPosList = NULL;
	numQryTuples = seqLength - seedLength + 1;
	for (i = 0; i < numQryTuples; ++i)
	{
		/* Step 2: Calculate the hash for each tuple. */
		/* Step 3: Find the reference coordinates of the hash in the hash
		 * table. */
		/* Step 4: Aggregate the reference coordinate positions. */
		strncpy(seed, seq + i, seedLength);
		tmpRefPosList = hashTable[getHash(seed, seedLength)];
		if (tmpRefPosList == NULL)
		{
			continue;
		}
		else if (tmpRefPosList->next == NULL)
		{
			tmpHitList = hitListCreateNode(
					tmpRefPosList->index,
					tmpRefPosList->position - i,
					tmpRefPosList->position);
			hitList = hitListPushBack(hitList, tail, tmpHitList);
			tail = tmpHitList;
		}
		else
		{
			/* If this tuple has more repeats than the threshold value,
			 * ignore the tuple. */
			tmpRefPosList2 = tmpRefPosList;
			count = 0;
			isCountHigh = 0;
			while (tmpRefPosList2 != NULL)
			{
				++count;
				if (count > TUPLE_IGNORE_THRES)
				{
					isCountHigh = 1;
					break;
				}
				tmpRefPosList2 = tmpRefPosList2->next;
			}
			if (isCountHigh == 1)
				continue;

			/* Add the tuples to the list. */
			while (tmpRefPosList != NULL)
			{
				tmpHitList = hitListCreateNode(
						tmpRefPosList->index,
						tmpRefPosList->position - i,
						tmpRefPosList->position);
				hitList = hitListPushBack(hitList, tail, tmpHitList);
				tail = tmpHitList;
				tmpRefPosList = tmpRefPosList->next;
			}
		}
	}


	/* Step 5: Sort the coordinate positions. */
	hitList = hitListSort(hitList, hitListGetSize(hitList));


	/* Step 6: Find the best-matching positions. */
	tmpHitList = getBestScoringPosition(hitList, maxNumHits);


	/* Step 7: Delete the HitList. */
	hitListDelete(hitList);


	/* Step 8: Return the best-matching positions. */
	return tmpHitList;
}


/**
 * Searches the paired query sequences in the reference.
 *
 * @param hashTable Hash table containing hash as key and a list of reference
 * positions as values.
 * @param qry1 First query object.
 * @param qry2 Second query object.
 * @param seqLen Sequence length.
 * @param seedLen Seed length to be used.
 * @param maxNumHits Max number of hits to be fetched.
 * @param minFragSize Minimum DNA fragment size used in sequencing.
 * @param maxFragSize Maximum DNA fragment size used in sequencing.
 * @param[out] result Array of @a HitList containing 8 elements. Element 1
 * will contain a list of hits for query 1; Element 2 will contain the
 * corresponding paired list of hits for query 2. Element 3 will contain a
 * list of hits for reverse complement of query 1; Element 4 will contain
 * the corresponding paired list of hits for query 2. Element 5 will contain a
 * list of hits for query 1; Element 6 will contain the corresponding paired
 * list of hits for the reverse complement of query 2. Element 6 will contain
 * a list of hits for reverse complement of query 1; Element 7 will contain
 * the corresponding paired list of hits for query 2. .
 */
void refSearchQuery_paired(RefPosList **hashTable, const Query *qry1,
		const Query *qry2, uint seqLen, uint seedLen, uint maxNumHits,
		uint minFragSize, uint maxFragSize, HitList **result)
{
	if (result == NULL || hashTable == NULL || qry1 == NULL || qry2 == NULL
			|| seqLen == 0 || seedLen == 0 || maxNumHits == 0 || minFragSize
			== 0 || maxFragSize == 0)
		return;

	static int isQry1RevComp, isQry2RevComp;

	/* read1_score + read2_score */
	isQry1RevComp = 0, isQry2RevComp = 0;
	getHits(hashTable, qry1, isQry1RevComp, qry2, isQry2RevComp, seqLen,
			seedLen, maxNumHits, minFragSize, maxFragSize, &(result[0]),
			&(result[1]));

	/* read1revcomp_score + read2_score */
	isQry1RevComp = 1, isQry2RevComp = 0;
	getHits(hashTable, qry1, isQry1RevComp, qry2, isQry2RevComp, seqLen,
			seedLen, maxNumHits, minFragSize, maxFragSize, &(result[2]),
			&(result[3]));

	/* read1_score + read2revcomp_score */
	isQry1RevComp = 0, isQry2RevComp = 1;
	getHits(hashTable, qry1, isQry1RevComp, qry2, isQry2RevComp, seqLen,
			seedLen, maxNumHits, minFragSize, maxFragSize, &(result[4]),
			&(result[5]));

	/* read1revcomp_score + read2revcomp_score */
	isQry1RevComp = 1, isQry2RevComp = 1;
	getHits(hashTable, qry1, isQry1RevComp, qry2, isQry2RevComp, seqLen,
			seedLen, maxNumHits, minFragSize, maxFragSize, &(result[6]),
			&(result[7]));
}


/**
 * Returns all hits for the given pair of queries.
 *
 * @param hashTable Hash table.
 * @param qry1 First @a Query object.
 * @param isQry1RevComp Indicates whether @a qry1 is a reverse complement.
 * @param qry2 Second @a Query object.
 * @param isQry2RevComp Indicates whether @a qry2 is a reverse complement.
 * @param seqLen Length of each query sequence.
 * @param seedLen Seed length.
 * @param maxNumHits Maximum number of hits to be fetched.
 * @param minFragSize Minimum DNA fragment size used in sequencing.
 * @param maxFragSize Maximum DNA fragment size used in sequencing.
 * @param[out] result1 @a HitList in which hits for @a qry1 should be stored.
 * @param[out] result2 @a HitList in which hits for @a qry2 should be stored.
 *
 */
static void getHits(RefPosList **hashTable, const Query *qry1,
		int isQry1RevComp, const Query *qry2, int isQry2RevComp, uint seqLen,
		uint seedLen, int maxNumHits, uint minFragSize, uint maxFragSize,
		HitList **result1, HitList **result2)
{
	static int numQryTuples;
	static int i, count, isCountHigh;
	static RefPosList *tmpRefPosList1, *tmpRefPosList2, *tmpRefPosList3;
	static HitList *hitList1, *hitList2, *tmpHitList1, *tmpHitList2, *tail;
	static char *qry1Seq, *qry2Seq;
	char qry1Seed[seedLen + 1], qry2Seed[seedLen + 1];

	/* Query 1 */
	if (isQry1RevComp == 0)
		qry1Seq = qry1->seq;
	else
		qry1Seq = qry1->revCompSeq;

	/* Break the query into tuples. */
	hitList1 = NULL, tail = NULL;
	tmpRefPosList1 = NULL, tmpRefPosList2 = NULL;
	numQryTuples = seqLen - seedLen + 1;
	for (i = 0; i < numQryTuples; ++i)
	{
		/* (1) Calculate the hash for each tuple.
		 * (2) Find the reference coordinates of the hash in the hash table.
		 * (3) Aggregate the reference coordinate positions. */
		strncpy(qry1Seed, qry1Seq + i, seedLen);
		tmpRefPosList1 = hashTable[getHash(qry1Seed, seedLen)];
		if (tmpRefPosList1 == NULL)
			continue;
		else if (tmpRefPosList1->next == NULL)
		{
			tmpHitList1 = hitListCreateNode(tmpRefPosList1->index,
					tmpRefPosList1->position - i, tmpRefPosList1->position);
			hitList1 = hitListPushBack(hitList1, tail, tmpHitList1);
			tail = tmpHitList1;
		}
		else
		{
			/* If this tuple has more repeats than the threshold value,
			 * ignore the tuple. */
			tmpRefPosList3 = tmpRefPosList1;
			count = 0;
			isCountHigh = 0;
			while (tmpRefPosList3 != NULL)
			{
				++count;
				if (count > TUPLE_IGNORE_THRES)
				{
					isCountHigh = 1;
					break;
				}
				tmpRefPosList3 = tmpRefPosList3->next;
			}
			if (isCountHigh == 1)
				continue;

			while (tmpRefPosList1->next != NULL)
			{
				tmpHitList1 = hitListCreateNode(tmpRefPosList1->index,
						tmpRefPosList1->position - i, tmpRefPosList1->position);
				hitList1 = hitListPushBack(hitList1, tail, tmpHitList1);
				tail = tmpHitList1;
				tmpRefPosList1 = tmpRefPosList1->next;
			}
			tmpHitList1 = hitListCreateNode(tmpRefPosList1->index,
					tmpRefPosList1->position - i, tmpRefPosList1->position);
			hitList1 = hitListPushBack(hitList1, tail, tmpHitList1);
			tail = tmpHitList1;
		}
	}
	hitList1 = hitListSort(hitList1, hitListGetSize(hitList1));
//	tmpHitList1 = getBestScoringPosition2(hitList1, -1);
	tmpHitList1 = getBestScoringPosition2(hitList1, maxNumHits);
	hitListDelete(hitList1);


	/* Query 2 */
	if (isQry2RevComp == 0)
		qry2Seq = qry2->seq;
	else
		qry2Seq = qry2->revCompSeq;

	/* Break the query into tuples. */
	hitList2 = NULL, tail = NULL;
	tmpRefPosList2 = NULL;
	for (i = 0; i < numQryTuples; ++i)
	{
		/* (1) Calculate the hash for each tuple.
		 * (2) Find the reference coordinates of the hash in the hash table.
		 * (3) Aggregate the reference coordinate positions. */
		strncpy(qry2Seed, qry2Seq + i, seedLen);
		tmpRefPosList2 = hashTable[getHash(qry2Seed, seedLen)];
		if (tmpRefPosList2 == NULL)
			continue;
		else if (tmpRefPosList2->next == NULL)
		{
			tmpHitList2 = hitListCreateNode(tmpRefPosList2->index,
					tmpRefPosList2->position - i, tmpRefPosList2->position);
			hitList2 = hitListPushBack(hitList2, tail, tmpHitList2);
			tail = tmpHitList2;
		}
		else
		{
			/* If this tuple has more repeats than the threshold value,
			 * ignore the tuple. */
			tmpRefPosList3 = tmpRefPosList2;
			count = 0;
			isCountHigh = 0;
			while (tmpRefPosList3 != NULL)
			{
				++count;
				if (count > TUPLE_IGNORE_THRES)
				{
					isCountHigh = 1;
					break;
				}
				tmpRefPosList3 = tmpRefPosList3->next;
			}
			if (isCountHigh == 1)
				continue;

			while (tmpRefPosList2->next != NULL)
			{
				tmpHitList2 = hitListCreateNode(tmpRefPosList2->index,
						tmpRefPosList2->position - i, tmpRefPosList2->position);
				hitList2 = hitListPushBack(hitList2, tail, tmpHitList2);
				tail = tmpHitList2;
				tmpRefPosList2 = tmpRefPosList2->next;
			}
			tmpHitList2 = hitListCreateNode(tmpRefPosList2->index,
					tmpRefPosList2->position - i, tmpRefPosList2->position);
			hitList2 = hitListPushBack(hitList2, tail, tmpHitList2);
			tail = tmpHitList2;
		}
	}
	hitList2 = hitListSort(hitList2, hitListGetSize(hitList2));
//	tmpHitList2 = getBestScoringPosition2(hitList2, -1);
	tmpHitList2 = getBestScoringPosition2(hitList2, maxNumHits);
	hitListDelete(hitList2);

	/* Include those hits from both queries that satisfy the contraints. */
	*result1 = NULL;
	*result2 = NULL;
	HitList *tmpHitListNode1, *tmpHitListNode2;
	hitList1 = tmpHitList1;
	HitList *tail1 = *result1;
	HitList *tail2 = *result2;
	while (hitList1 != NULL)
	{
		hitList2 = tmpHitList2;
		while (hitList2 != NULL)
		{
			if (((hitList2->offset + seqLen) - hitList1->offset) <= maxFragSize
					&& ((hitList2->offset + seqLen) - hitList1->offset)
							>= minFragSize)
			{
				tmpHitListNode1 = hitListDuplicateNode(hitList1);
				*result1 = hitListPushBack(*result1, tail1, tmpHitListNode1);
				tail1 = tmpHitListNode1;

				tmpHitListNode2 = hitListDuplicateNode(hitList2);
				*result2 = hitListPushBack(*result2, tail2, tmpHitListNode2);
				tail2 = tmpHitListNode2;
			}
			hitList2 = hitList2->next;
		}
		hitList1 = hitList1->next;
	}

	hitListDelete(tmpHitList1);
	hitListDelete(tmpHitList2);
}


/**
 * Returns the positions of the best-matching reference sequences.
 *
 * @param list The HitList list from which best-matching positions are to
 * be found.
 * @param maxNumHits Max number of hits to be returned.
 * @return List of best-matching reference coordinates.
 */
static HitList *getBestScoringPosition(HitList *list, uint maxNumHits)
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
		uint numRandomNumbers = maxNumHits;
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
 * Returns the positions of the best-matching reference sequences.
 *
 * @param list The HitList list from which best-matching positions are to
 * be found.
 * @param maxNumHits Max number of hits to be returned.
 * @return List of best-matching reference coordinates.
 */
static HitList *getBestScoringPosition2(HitList *list, int maxNumHits)
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
			fprintf(stderr, "Error allocating emough memory in %s at "
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
		int numRandomNumbers = maxNumHits;
		if (numRandomNumbers > numMaxSizeClusters)
			numRandomNumbers = numMaxSizeClusters;
		HitList *bestMatches2 = NULL, *tmpBestMatches2 = NULL;
		if (numRandomNumbers >= 0)
		{
			uint randomNumArr[numRandomNumbers];
			arrGetRandomNums(numRandomNumbers, 1, numMaxSizeClusters,
					randomNumArr);
			qsort(randomNumArr, numRandomNumbers, sizeof(uint), compare);

			tail = NULL;
			int l = 0;
			for (j = 0; j < numClusters; ++j)
			{
				if (clusterSizeArr[j] == maxClusterSize)
				{
					++l;
					if (arrSearch((int *) randomNumArr, numRandomNumbers, l)
							> -1)
					{
						tmpBestMatches2 = hitListDuplicateNode(bestMatches[j]);
						bestMatches2 = hitListPushBack(bestMatches2, tail,
								tmpBestMatches2);
						tail = tmpBestMatches2;
					}
				}
			}
		}
		else
		{
			tail = NULL;
			for (j = 0; j < numClusters; ++j)
			{
				if (clusterSizeArr[j] == maxClusterSize)
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


/**
 * This is a GPU kernel function. It searches the query sequence in the
 * reference and returns the best matching reference coordinates.
 *
 * @param keys Array containing index to @a vals array.
 * @param vals Array containing reference sequence offset.
 * @param numKeys Number of elements in @a keys. This number must be equal
 * to the number of elements in @a keys, otherwise it will lead to
 * invalid results.
 * @param numVals Number of elements in @a vals. This number must be equal
 * to the number of elements in @a vals, otherwise it will lead to
 * invalid results.
 * @param qrySeq Array containing query sequences.
 * @param[out] results Array containing matching reference sequence offsets.
 * @param qryLen Query sequence length.
 * @param seedLength Seed length to be used.
 * @param maxNumHits Max number of hits to be returned.
 * @param randNum A random number.
 */
__global__ void refSearchQuery2_kernel(uint *keys, uint *vals, int numKeys,
		int numVals, char *qrySeq, char *refIdxArr2, int *shiftArr2,
		int *posArr2, int qryLen, int seedLength, int maxNumHits, int randNum)
{
	__shared__ char seed[SEED_LENGTH];
	__shared__ char refIdxArr[MAX_HITS];
	__shared__ int shiftArr[MAX_HITS];
	__shared__ int posArr[MAX_HITS];

	long blkId = (blockIdx.y * gridDim.x) + blockIdx.x;
	randNumSeed = blkId + randNum;
	int numQryTuples = qryLen - seedLength + 1;
	int i, j, l = 0, hash;
	uint k, key1, key2;
	qrySeq += (blkId * qryLen);
	uint refIdxShiftBits = (sizeof(uint) * NUM_BITS_IN_A_BYTE) - REF_IDX_BITS;
	uint tupleIdxMask = pow_gpu2(2, POS_BITS) - 1;

	for (i = 0; i < MAX_HITS; ++i)
	{
		refIdxArr[i] = -1;
		shiftArr[i] = -1;
		posArr[i] = -1;
	}

	/* Step 1: Find matches and put them in array. */
	for (i = 0; i < numQryTuples; ++i)
	{
		strncpy_gpu((char *) seed, qrySeq + i, seedLength);
		hash = getHash_gpu2((char *) seed, seedLength);
		key1 = keys[hash];

		if (key1 != UINT_MAX)
		{
			j = 0;
			key2 = UINT_MAX;
			while (key2 == UINT_MAX)
			{
				++j;
				key2 = keys[hash + j];
			}

			for (k = key1; k < key2; ++k)
			{
				refIdxArr[l] = (char) (vals[k] >> refIdxShiftBits);
				posArr[l] = (vals[k] & tupleIdxMask) * seedLength;
				shiftArr[l] = posArr[l] - i;
				++l;

				insertionSort3_gpu(refIdxArr, shiftArr, posArr, l);
			}
		}
	}

	/* Step 2: Sort arrays. */
//	bubbleSort2(refIdxArr, shiftArr, posArr, l);
//	quickSort2_gpu(refIdxArr, shiftArr, posArr, l);

	/* Step 3: Create clusters. */
	/* Step 4: Pick the biggest clusters. */
	/* Step 5: Return the first positions of the biggest clusters. */
	refIdxArr2 += (blkId * maxNumHits);
	shiftArr2 += (blkId * maxNumHits);
	posArr2 += (blkId * maxNumHits);
	getBestMatches2(refIdxArr, shiftArr, posArr, l, maxNumHits, refIdxArr2,
			shiftArr2, posArr2);
}


/**
 * This is a wrapper function that wraps @a getBestMatches2 function. This
 * function has been created so that @a getBestMatches2 function can be
 * unit-tested.
 *
 * This function takes 3 arrays as input, creates clusters,
 * finds the biggest cluster, and returns the first hit in the biggest
 * cluster. If there is more than one biggest cluster, it selects one biggest
 * cluster randomly and returns the first hit from that cluster.
 *
 * @param	refIdxArr	Array containing the reference indexes. This array
 * is assumed to be sorted.
 * @param	shiftArr	Array containing the shifts. This array is assumed to
 * be sorted according to @a refIdxArr.
 * @param	posArr		Array containing reference tuple positions. This array
 * is assumed to be sorted according to @a refIdxArr and @a shiftArr.
 * @param	size		Number of elements in each of the array. This number
 * must be equal to the number of elements in each of @a refIdxArr,
 * @a shiftArr, and @a posArr, otherwise it will lead to invalid results.
 * @param	maxNumMatches	Maximum number of matches to be returned.
 * The maximum value this parameter can currently be is 20.
 * @param[out]	refIdxArr2	This array will store the reference index of the
 * best matching hits.
 * @param[out]	shiftArr2	This array will store the shift of the best
 * matching hits.
 * @param[out]	posArr2		This array will store the reference tuple positions
 * of the best matching hits.
 */
__global__ void getBestMatches2_wrap(char *refIdxArr, int *shiftArr,
		int *posArr, int size, int maxNumMatches,
		char *refIdxArr2, int *shiftArr2, int *posArr2)
{
	randNumSeed = (uint) clock();
	getBestMatches2(refIdxArr, shiftArr, posArr, size, maxNumMatches,
			refIdxArr2, shiftArr2, posArr2);
}


/**
 * This function takes 3 arrays as input, creates clusters,
 * finds the biggest cluster, and returns the first hit in the biggest
 * cluster. If there is more than one biggest cluster, it selects one biggest
 * cluster randomly and returns the first hit from that cluster.
 *
 * @param	refIdxArr	Array containing the reference indexes. This array
 * is assumed to be sorted.
 * @param	shiftArr	Array containing the shifts. This array is assumed to
 * be sorted according to @a refIdxArr.
 * @param	posArr		Array containing reference tuple positions. This array
 * is assumed to be sorted according to @a refIdxArr and @a shiftArr.
 * @param	size		Number of elements in each of the array. This number
 * must be equal to the number of elements in each of @a refIdxArr,
 * @a shiftArr, and @a posArr, otherwise it will lead to invalid results.
 * @param	maxNumMatches	Maximum number of matches to be returned.
 * The maximum value this parameter can currently be is 20.
 * @param[out]	refIdxArr2	This array will store the reference index of the
 * best matching hits.
 * @param[out]	shiftArr2	This array will store the shift of the best
 * matching hits.
 * @param[out]	posArr2		This array will store the reference tuple positions
 * of the best matching hits.
 */
__device__ void getBestMatches2(char *refIdxArr, int *shiftArr, int *posArr,
		int size, int maxNumMatches, char *refIdxArr2, int *shiftArr2,
		int *posArr2)
{
	if (size == 0)
	{
		int i;
		for (i = 0; i < maxNumMatches; ++i)
		{
			refIdxArr2[i] = -1;
			shiftArr2[i] = -1;
			posArr2[i] = -1;
		}
		return;
	}
	else if (size == 1)
	{
		refIdxArr2[0] = refIdxArr[0];
		shiftArr2[0] = shiftArr[0];
		posArr2[0] = posArr[0];

		int i;
		for (i = 1; i < maxNumMatches; ++i)
		{
			refIdxArr2[i] = -1;
			shiftArr2[i] = -1;
			posArr2[i] = -1;
		}
		return;
	}

	__shared__ ushort clusterStart[MAX_CLUSTER_SIZE];
	__shared__ ushort numHits[MAX_CLUSTER_SIZE];
	__shared__ uint randNums[MAX_RAND_NUMS];
	__shared__ ushort idxArr[MAX_CLUSTER_SIZE];

	/* Create clusters. */
	clusterStart[0] = 0;
	numHits[0] = 1;
	int i, j = 0;
	for (i = 1; i < size; ++i)
	{
		if ((refIdxArr[i] > refIdxArr[i - 1])
				|| (refIdxArr[i] == refIdxArr[i - 1]
				    && shiftArr[i] > shiftArr[i - 1]))
		{
			++j;
			clusterStart[j] = i;
			numHits[j] = 1;
		}
		else
			++numHits[j];
	}
	int numClusters = j + 1;

	/* Find the biggest cluster size. */
	int biggestClusterSize = 0;
	for (i = 0; i < numClusters; ++i)
		biggestClusterSize = max(biggestClusterSize, numHits[i]);

	/* Copy first hit of all biggest clusters to an array. */
	j = -1;
	for (i = 0; i < numClusters; ++i)
	{
		if (numHits[i] == biggestClusterSize)
		{
			++j;
			idxArr[j] = clusterStart[i];
		}
	}
	int numMaxClusters = j + 1;
	int numMatches = min(maxNumMatches, numMaxClusters);

	/* Randomly choose the required number of hits from among all the hits. */
	arrGetRandomNums_gpu(numMatches, 0, numMaxClusters - 1, randNums);
	for (i = 0; i < numMatches; ++i)
	{
		refIdxArr2[i] = refIdxArr[idxArr[randNums[i]]];
		shiftArr2[i] = shiftArr[idxArr[randNums[i]]];
		posArr2[i] = posArr[idxArr[randNums[i]]];
	}

	/* Assign '-1' to remaining elements. */
	for (i = numMatches; i < maxNumMatches; ++i)
	{
		refIdxArr2[i] = -1;
		shiftArr2[i] = -1;
		posArr2[i] = -1;
	}
}


/**
 * This is a wrapper function that wraps @a arrGetRandomNums_gpu
 * function. It's used to unit-test @a arrGetRandomNums_gpu function.
 *
 * @param n Number of random numbers to be created. This number should not
 * be greater than the range specified by lower and upper limit.
 * @param lowerLimit The lower limit of the random numbers.
 * @param upperLimit The upper limit of the random numbers.
 * @param[out] arr Array in which the random numbers will be stored. The size
 * of the array should be atleast @a n.
 */
__global__ void arrGetRandomNums_gpu_wrap(uint n, uint lowerLimit,
		uint upperLimit, uint *arr)
{
	randNumSeed = (uint) clock();
	arrGetRandomNums_gpu(n, lowerLimit, upperLimit, arr);
}


/**
 * Fetches the given number of random numbers between the given limits.
 *
 * @param n Number of random numbers to be created. This number should not
 * be greater than the range specified by lower and upper limit.
 * @param lowerLimit The lower limit of the random numbers.
 * @param upperLimit The upper limit of the random numbers.
 * @param[out] arr Array in which the random numbers will be stored. The size
 * of the array should be atleast @a n.
 */
__device__ void arrGetRandomNums_gpu(uint n, uint lowerLimit, uint upperLimit,
		uint *arr)
{
	int range = upperLimit - lowerLimit + 1;
	if (n <= 0 || lowerLimit > upperLimit || arr == NULL || n > range)
		return;

	int i = 0;
	uint randNum;
	while (i < n)
	{
		randNum = (uint) ((getRandNum() % range) + lowerLimit);
		if (arrSearch_gpu(arr, i, randNum) == -1)
		{
			arr[i] = randNum;
			++i;
		}
	}
}


/**
 * This is a wrapper function that wraps @a getRandNum function.
 * This function is used so that @a getRandNum can be unit-tested.
 *
 * This function generates random numbers and stores them in @a arr.
 *
 * @param[out] arr Array where random numbers will be stored.
 * The number of elements in this array must be equal to @a size.
 * @param size Number of elements in @a arr.
 */
__global__ void getRandNum_wrap(uint *arr, int size)
{
	randNumSeed = (uint) clock();

	int i;
	for (i = 0; i < size; ++i)
		arr[i] = getRandNum();
}


/**
 * Returns a pseudo-random number.
 *
 * @return Random number.
 * @note This algorithm has been taken from "The C Programming Language"
 * by Kernighan and Ritchie.
 */
__device__ uint getRandNum()
{
	randNumSeed = (randNumSeed * 1103515245) + 12345;
	return ((randNumSeed / 65536) % 32768);
}


/**
 * This is a wrapper function that wraps @a arrSearch_gpu
 * function. It's added so that @a arrSearch_gpu can be
 * unit-tested.
 *
 * @param arr Array in which the number will be searched.
 * @param arrSize Number of elements in the array. This number must be equal
 * to the number of elements in @a arr, otherwise it may result in invalid
 * behavior.
 * @param num The number to be searched in the array.
 * @param index Index of the searched number in the array; otherwise, -1.
 */
__global__ void arrSearch_gpu_wrap(uint *arr, int arrSize, uint num, int *index)
{
	*index = arrSearch_gpu(arr, arrSize, num);
}


/**
 * Returns the index of the given number if it is already present in the
 * given array; otherwise, returns -1.
 *
 * @param arr Array in which the number will be searched.
 * @param arrSize Number of elements in the array. This number must be
 * equal to the number of elements in @a arr, otherwise it may result in invalid
 * behavior.
 * @param num The number to be searched in the array.
 * @return Index of the searched number in the array; otherwise, -1.
 */
__device__ int arrSearch_gpu(uint *arr, int arrSize, uint num)
{
	if (arr == NULL || arrSize == 0)
		return -1;

	int i = 0;
	for (i = 0; i < arrSize; ++i)
	{
		if (arr[i] == num)
			return i;
	}

	return -1;
}


/**
 * This is a wrapper function that wraps @a bubbleSort2. This function has
 * been added so that @a bubbleSort2 can be unit-tested.
 *
 * Sorts the given arrays using Bubblesort. @a refIdxArr, @a shiftArr, and
 * @a posArr are parallel arrays. The arrays are first sorted by @a refIdxArr,
 * then by @a shiftArr, and finally by @a posArr.
 *
 * @param[in,out] refIdxArr Array containing reference index.
 * @param[in,out] shiftArr Array containing shifts.
 * @param[in,out] posArr Array containing reference tuple positions.
 * @param size Number of elements in each of @a refIdxArr, @a shiftArr, and
 * @a posArr. This number must be equal to number of elements in each of the
 * 3 arrays, otherwise it may lead to invalid results.
 */
__global__ void bubbleSort2_wrap(char *refIdxArr, int *shiftArr, int *posArr,
		uint size)
{
	bubbleSort2(refIdxArr, shiftArr, posArr, size);
}


/**
 * Sorts the given arrays using Bubblesort. @a refIdxArr, @a shiftArr, and
 * @a posArr are parallel arrays. The arrays are first sorted by @a refIdxArr,
 * then by @a shiftArr, and finally by @a posArr.
 *
 * @param[in,out] refIdxArr Array containing reference index.
 * @param[in,out] shiftArr Array containing shifts.
 * @param[in,out] posArr Array containing reference tuple positions.
 * @param size Number of elements in each of @a refIdxArr, @a shiftArr, and
 * @a posArr. This number must be equal to number of elements in each of the
 * 3 arrays, otherwise it may lead to invalid results.
 */
__device__ void bubbleSort2(char *refIdxArr, int *shiftArr, int *posArr,
		uint size)
{
	char tmpRefIdx;
	int tmpShift, tmpPos, i, j;

	for (i = 0; i < size; ++i)
	{
		for (j = size - 1; j > i; --j)
		{
			if (refIdxArr[j - 1] > refIdxArr[j])
			{
				tmpRefIdx = refIdxArr[j];
				tmpShift = shiftArr[j];
				tmpPos = posArr[j];

				refIdxArr[j] = refIdxArr[j - 1];
				shiftArr[j] = shiftArr[j - 1];
				posArr[j] = posArr[j - 1];

				refIdxArr[j - 1] = tmpRefIdx;
				shiftArr[j - 1] = tmpShift;
				posArr[j - 1] = tmpPos;
			}
			else if (refIdxArr[j - 1] == refIdxArr[j]
			             && shiftArr[j - 1] > shiftArr[j])
			{
				tmpShift = shiftArr[j];
				tmpPos = posArr[j];

				shiftArr[j] = shiftArr[j - 1];
				posArr[j] = posArr[j - 1];

				shiftArr[j - 1] = tmpShift;
				posArr[j - 1] = tmpPos;
			}
			else if (refIdxArr[j - 1] == refIdxArr[j]
			             && shiftArr[j - 1] == shiftArr[j]
			             && posArr[j - 1] > posArr[j])
			{
				tmpPos = posArr[j];
				posArr[j] = posArr[j - 1];
				posArr[j - 1] = tmpPos;
			}
		}
	}
}


/**
 * This is a wrapper function that wraps @a quickSort2_gpu function. This
 * function has been added so that @a quickSort2_gpu can be unit-tested.
 *
 * @param[in,out]	refIdxArr	Reference index array.
 * @param[in,out]	shiftArr	Shift array.
 * @param[in,out]	posArr		Position array.
 * @param			arrSize		Number of elements in each of @a refIdxArr,
 * @a shiftArr, and @a posArr.
 */
__global__ void quickSort2_gpu_wrap(char *refIdxArr, int *shiftArr, int *posArr,
		int arrSize)
{
	quickSort2_gpu(refIdxArr, shiftArr, posArr, arrSize);
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
 *
 * @note This code has been adapted from the code written by Darel Rex Finley
 * posted on the site http://alienryderflex.com/quicksort/.
 */
__device__ int quickSort2_gpu(char *refIdxArr, int *shiftArr, int *posArr,
		int arrSize)
{
#define  MAX_LEVELS  MAX_HITS
	__shared__ int begin[MAX_LEVELS];
	__shared__ int end[MAX_LEVELS];
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

			if (i == MAX_LEVELS - 1)
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
 * This is a wrapper function that wraps @a insertionSort3_gpu function. This
 * function has been added so that @a insertionSort3_gpu can be unit-tested.
 *
 * @param[in,out]	refIdxArr	Reference index array.
 * @param[in,out]	shiftArr	Shift array.
 * @param[in,out]	posArr		Position array.
 * @param			arrSize		Number of elements in each of @a refIdxArr,
 * @a shiftArr, and @a posArr.
 */
__global__ void insertionSort3_gpu_wrap(char *refIdxArr, int *shiftArr,
		int *posArr, int arrSize)
{
	insertionSort3_gpu(refIdxArr, shiftArr, posArr, arrSize);
}


/**
 * Performs insertion sort on the given arrays, first by @a refIdxArr, then
 * by @a shiftArr, and finally by @a posArr. The last element of the array is
 * the new element to be inserted.
 *
 * @note The array size is assumed to be at least 2, otherwise this function
 * will lead to incorrect results.
 *
 * @param[in,out]	refIdxArr	Reference index array.
 * @param[in,out]	shiftArr	Shift array.
 * @param[in,out]	posArr		Position array.
 * @param			arrSize		Number of elements in each of @a refIdxArr,
 * @a shiftArr, and @a posArr.
 */
__device__ void insertionSort3_gpu(char *refIdxArr, int *shiftArr, int *posArr,
		int arrSize)
{
	int i = arrSize - 2;
	int newRefIdx = refIdxArr[arrSize - 1];
	int newShift = shiftArr[arrSize - 1];
	int newPos = posArr[arrSize - 1];
	while (i >= 0 && ((newRefIdx < refIdxArr[i])
			|| (newRefIdx == refIdxArr[i] && newShift < shiftArr[i])
			|| (newRefIdx == refIdxArr[i] && newShift == shiftArr[i]
			     && newPos < posArr[i])))
	{
		--i;
	}

	int j;
	for (j = arrSize - 1; j > (i + 1); --j)
	{
		refIdxArr[j] = refIdxArr[j - 1];
		shiftArr[j] = shiftArr[j - 1];
		posArr[j] = posArr[j - 1];
	}
	refIdxArr[i + 1] = newRefIdx;
	shiftArr[i + 1] = newShift;
	posArr[i + 1] = newPos;
}


/**
 * This is a wrapper function that wraps @a getHash_gpu2 function. This function
 * has been added to unit-test @a getHash_gpu2 function.
 *
 * @param s Character string.
 * @param length Length of characters in @a s for which hash is to be created.
 * @param[out] hash Hash created from @a s.
 */
__global__ void getHash_gpu2_wrap(char *s, int length, int *hash)
{
	*hash = getHash_gpu2(s, length);
}


/**
 * Returns the hash value of the given string
 *
 * @param s Given string
 * @param length Length of the given string
 * @return Hash value
 */
__device__ static int getHash_gpu2(char *s, int length)
{
	int i, sum = 0;

	/* This section of code encodes the string into an integer.
	 * We use '4' as the base in the power function below because we
	 * have only 4 letters in our alphabet, i.e. A, C, G, and T.
	 *
	 * Example: Encoded value 'ACG' will be:
	 * = (4^0 * CODE_A) + (4^1 * CODE_C) + (4^2 * CODE_G)
	 * = (1 * CODE_A) + (4 * CODE_C) + (16 * CODE_G)
	 * (Assume CODE_A = 0, CODE_C = 1, and CODE_G = 2)
	 * = (1 * 0) + (4 * 1) + (16 * 2)
	 * = 0 + 4 + 32
	 * = 36
	 * That is, hash(ACG) = 36 */
	for (i = 0; i < length; ++i)
	{
		if (s[i] == 'A' || s[i] == 'a')
			sum += (int) pow((float) 4, i) * CODE_A;
		else if (s[i] == 'C' || s[i] == 'c')
			sum += (int) pow((float) 4, i) * CODE_C;
		else if (s[i] == 'G' || s[i] == 'g')
			sum += (int) pow((float) 4, i) * CODE_G;
		else if (s[i] == 'T' || s[i] == 't')
			sum += (int) pow((float) 4, i) * CODE_T;
		else if (s[i] == 'N' || s[i] == 'n')
			sum += (int) pow((float) 4, i) * CODE_N;
		else if (s[i] == '*')
			sum += (int) pow((float) 4, i) * CODE_STAR;
		else
			sum += 0;
	}

	return sum;
}


/**
 * This is a wrapper function that wraps @a serializeHashTable2 function. It
 * has been added so that @a serializeHashTable2 function can be unit-tested.
 *
 * This function serializes the given hash table. It stores the keys
 * in one array and all the values in another array.
 *
 * @param hashTable Hash table to be serialized.
 * @param numKeys Number of keys in the hash table. This number must be equal
 * to the number of keys in @a hashTable, otherwise this will lead to
 * incorrect results.
 * @param numPos Number of values in the hash table.
 * @param[out] keys Array contains the keys from the hash table. Each key is
 * the index into the @a vals array. @a keys should point to enough memory to
 * store @a numKeys elements, otherwise it will lead to invalid behavior.
 * @param[out] vals Array containing the values of the hash table. @a vals
 * should point to enough memory to store @a numPos elements, otherwise
 * it will lead to invalid behavior. Left-most REF_IDX_BITS bits represent
 * reference index; Middle SHIFT_BITS bits represent shift; Right-most POS_BITS
 * bits represent reference position of the tuple.
 * @param seedLen Seed length.
 * @param tupleIgnoreThres Threshold used to ignore tuples. If the number of
 * repeats of a tuple is more than this number, then that tuple will be
 * ignored and not serialized.
 */
void serializeHashTable2_wrap(RefPosList **hashTable, uint numKeys,
		uint numPos, uint *keys, uint *vals, uint seedLen,
		uint tupleIgnoreThres)
{
	serializeHashTable2(hashTable, numKeys, numPos, keys, vals, seedLen,
			tupleIgnoreThres);
}


/**
 * This function serializes the given hash table. It stores the keys
 * in one array and all the values in another array.
 *
 * @param		hashTable	Hash table to be serialized.
 * @param		numKeys		Number of keys in the hash table. This number
 * must be equal to the number of keys in @a hashTable, otherwise this will
 * lead to incorrect results.
 * @param		numPos		Number of values in the hash table.
 * @param[out]	keys 		Array contains the keys from the hash table.
 * Each key is the index into the @a vals array. @a keys should point to
 * enough memory to store @a numKeys elements, otherwise it will lead to
 * invalid behavior.
 * @param[out]	vals		Array containing the values of the hash table.
 * @a vals should point to enough memory to store @a numPos elements, otherwise
 * it will lead to invalid behavior. Left-most REF_IDX_BITS bits represent
 * reference index; Middle SHIFT_BITS bits represent shift; Right-most POS_BITS
 * bits represent reference position of the tuple.
 * @param		seedLen		Seed length.
 * @param		tupleIgnoreThres	Threshold used to ignore tuples. If the
 * number of repeats of a tuple is more than this number, then that tuple
 * will be ignored and not serialized.
 */
static void serializeHashTable2(RefPosList **hashTable, uint numKeys,
		uint numPos, uint *keys, uint *vals, uint seedLen,
		uint tupleIgnoreThres)
{
	RefPosList *tmpRefPosList;
	uint index = 0;
	int i;
	uint refIdx, tupleIdx;

	/* Find how many times a tuple appears in the hash table. */
	uint *numRepeats = (uint *) malloc(numKeys * sizeof(uint));
	int numKeys2 = numKeys - 1;
	FILE *repeatFile = fopen("/data/gpusw/numRepeats.txt", "w");
	for (i = 0; i < numKeys2; ++i)
	{
		tmpRefPosList = hashTable[i];
		numRepeats[i] = 0;
		while (tmpRefPosList != NULL)
		{
			++numRepeats[i];
			tmpRefPosList = tmpRefPosList->next;
		}
//		if (numRepeats[i] > 0)
//			fprintf(repeatFile, "%d\n", numRepeats[i]);
	}

	/* Serialize the hash table. */
	uint shiftBits = (sizeof(uint) * NUM_BITS_IN_A_BYTE) - REF_IDX_BITS;
	for (i = 0; i < numKeys2; ++i)
	{
		tmpRefPosList = hashTable[i];
		if (tmpRefPosList == NULL || numRepeats[i] > tupleIgnoreThres)
		{
			keys[i] = UINT_MAX;
			continue;
		}
		else
			keys[i] = index;

		while (tmpRefPosList != NULL)
		{
			refIdx = (uint) tmpRefPosList->index;
			refIdx = refIdx << shiftBits;
			tupleIdx = (uint) tmpRefPosList->position / seedLen;
			vals[index] = (refIdx | tupleIdx);

			++index;
			tmpRefPosList = tmpRefPosList->next;
		}
	}
	keys[numKeys - 1] = index;

	free(numRepeats);
}


/**
 * Prints serialized hash table.
 *
 * @param	keys	Array containing keys of hash table.
 * @param	numKeys	Number of elements in @a keys.
 * @param	vals	Array containing values of hash table.
 * @param	numVals	Number of elements in @a vals.
 * @param	seedLen	Seed length.
 */
void refPrintSerialHashTable(uint *keys, uint numKeys, uint *vals, uint numVals,
		int seedLen)
{
	int i, j, k, hash, tmpHash;
	uint key1, key2;
	char s[seedLen + 1];
	for (i = 0; i < numKeys; ++i)
	{
		key1 = keys[i];
//		if (i == 0)
//			fprintf(stderr, "i = 0, key1 = %u\n", keys[i]);
		if (key1 != UINT_MAX)
		{
			j = 0;
			hash = i;
			key2 = UINT_MAX;
			while (key2 == UINT_MAX)
			{
				++j;
				tmpHash = hash + j;
				if (tmpHash >= numKeys)
					break;
				key2 = keys[hash + j];
			}

			if (key2 == UINT_MAX)
				key2 = numVals;

//			if (i == 0)
//				fprintf(stderr, "i is really = 0\n");

			getUnhash(hash, s, seedLen);
			s[seedLen] = '\0';
			fprintf(stderr, "[seed=%s; key1=%d; key2=%d]:\t", s, key1, key2);
			for (k = key1; k < key2; ++k)
				fprintf(stderr, "%d\t", vals[k]);
			fprintf(stderr, ";\n");
		}
	}
}


/**
 * Prints the given hash table.
 *
 * @param	hashTable	Hash table to be printed.
 * @param	numKeys		Number of keys in @a hashTable.
 * @param	seedLen		Seed length.
 */
void refPrintHashTable(RefPosList **hashTable, int numKeys, int seedLen)
{
	int i;
	RefPosList *tmpRefPosList = NULL;
	char s[seedLen + 1];
	for (i = 0; i < numKeys; ++i)
	{
		tmpRefPosList = hashTable[i];
		if (tmpRefPosList != NULL)
		{
			getUnhash(i, s, seedLen);
			s[seedLen] = '\0';
			fprintf(stderr, "%s\t", s);
			while (tmpRefPosList != NULL)
			{
				fprintf(stderr, "(%d, %d); ", tmpRefPosList->index,
						tmpRefPosList->position);
				tmpRefPosList = tmpRefPosList->next;
			}
			fprintf(stderr, "\n");
		}
	}
}


/**
 * Frees memory occupied by hashTable.
 *
 * @param hashTable Hash table that is to be deleted.
 * @param numKeys Number of keys in the hash table. This number must be equal
 * to the number of keys in the hash table.
 */
void refDeleteHashTable(RefPosList **hashTable, uint numKeys)
{
	if (hashTable != NULL)
	{
		int i;
		for (i = 0; i < numKeys; ++i)
			refPosListDelete(hashTable[i]);
		free(hashTable);
	}
}


/**
 * This function wraps @a strncpy_gpu function. This function has been added
 * to unit-test @a strncpy_gpu function.
 *
 * @param t String where characters from @a s will be copied to.
 * @param s String where characters will be copied from.
 * @param n Number of characters to be copied from @a s to @a t.
 */
__global__ void strncpy_gpu_wrap(char *t, char *s, int n)
{
	strncpy_gpu(t, s, n);
}


/**
 * Copies @a length characters of @a s to @a t. A NULL character will be
 * appended to the end of @a t.
 *
 * @param t String where characters from @a s will be copied to. A NULL
 * character will be appended to the end of this string. Length
 * of this string must be at least @a n + 1 so that it can store @n characters
 * plus a NULL character.
 * @param s String where characters will be copied from.
 * @param n Number of characters to be copied from @a s to @a t.
 *
 * @note This piece of code has been adapted from "The C Programming
 * Language" book written by Kernighan and Ritchie.
 */
__device__ void strncpy_gpu(char *t, char *s, int n)
{
	while (n > 0 && *s != '\0')
	{
		*t = *s;
		++s;
		++t;
		--n;
	}
	*t = '\0';
}


/**
 * Calculates the value of base raised to the n-th power.
 *
 * This is a wrapper function that wraps @a pow_gpu2 function. This function
 * has been added so that @a pow_gpu2 can be unit-tested.
 *
 * @param 		base	Base value.
 * @param 		n		Exponent value.
 * @param[out]	p		The calculated value.
 */
__global__ void pow_gpu_wrap(int base, int n, int *p)
{
	*p = pow_gpu2(base, n);
}


/**
 * Calculates the value of base raised to the n-th power.
 *
 * @param 	base	Base value.
 * @param 	n		Exponent value.
 * @return			The calculated value.
 */
__device__ int pow_gpu2(int base, int n)
{
	int p = 1;
	while (n > 0)
	{
		p = p * base;
		--n;
	}

	return p;
}

