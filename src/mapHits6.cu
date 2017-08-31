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
#include "mapHits6.h"
#include "array.h"
#include <stdio.h>


/**
 * Returns the best mapped hits.
 *
 * @param 		numBestHits			Number of best hits to be returned.
 * @param		refIdx				Reference index.
 * @param		shift				Shifts.
 * @param		refPos				Reference positions.
 * @param		size				Number of elements in @a refIdx.
 * @param[out]	refIdx_bestMatch	Best-matching reference index. The number of
 * elements in this array must be equal to @a numBestHits.
 * @param[out]	shift_bestMatch		Best-matching shift. The number of elements
 * in this array must be equal to @a numBestHits.
 * @param[out]	refPos_bestMatch	Best-matching reference position. The
 * number of elements in this array must be equal to @a numBestHits.
 * @return		Number of best-hits.
 */
int mapHits6GetBestHits(int numBestHits, char *refIdx, int *shift, int *refPos,
		int size, char *refIdx_bestMatch, int *shift_bestMatch,
		int *refPos_bestMatch)
{
	/* Return 0 if there are no hits. */
	if (size == 0)
		return 0;

	/* Sort the list of mapped hits. */
	quickSort(0, size - 1, refIdx, shift, refPos);

	/* Create clusters. */
	int *cluster = (int *) calloc(size, sizeof(int));
	int numClusters = 0;
	createClusters(refIdx, shift, refPos, size, cluster, &numClusters);

	/* Find the biggest cluster size. */
	int *clusterSize = (int *) calloc(size, sizeof(int));
	int biggestClusterSize = 0;
	findBiggestClusterSize(cluster, numClusters, clusterSize, &biggestClusterSize,
			size);

	/* Find the first hit of the biggest clusters. */
	numBestHits = findBiggestClusters(numBestHits, refIdx_bestMatch,
			shift_bestMatch, refPos_bestMatch, biggestClusterSize, cluster,
			numClusters, clusterSize, refIdx, shift, refPos, size);

	free(cluster);
	free(clusterSize);
	return numBestHits;
}


/**
 * This is a wrapper function that wraps @a findBiggestClusters function. It
 * has been added so that @a findBiggestClusters can be unit-tested.
 *
 * @param 		numBestHits			Number of biggest clusters that are to
 * searched for.
 * @param[out]	refIdx_bestMatch	Best-matching reference index. The number
 * of elements in this array must be equal to @a numBestHits.
 * @param[out]	shift_bestMatch		Best-matching shift. The number of elements
 * in this array must be equal to @a numBestHits.
 * @param[out]	refPos_bestMatch	Best-matching reference position. The number
 * of elements in this array must be equal to @a numBestHits.
 * @param		refIdx				Reference Indices.
 * @param		shift				Shifts (reference position - query position).
 * @param		refPos				Reference position.
 * @param		size				Number of elements in @a refIdx.
 * @return 		Number of best-hits in @a refIdx_bestMatch.
 */
int findBiggestClusters_wrap6(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch, char *refIdx, int *shift,
		int *refPos, int size)
{
	int numClusters = 0;
	int *clusters = (int *) calloc(size, sizeof(int));
	createClusters(refIdx, shift, refPos, size, clusters, &numClusters);
	int biggestClusterSize = 0;
	int *clusterSize = (int *) calloc(size, sizeof(int));
	findBiggestClusterSize(clusters, numClusters, clusterSize,
			&biggestClusterSize, size);
	int numHits = findBiggestClusters(numBestHits, refIdx_bestMatch,
			shift_bestMatch, refPos_bestMatch, biggestClusterSize, clusters,
			numClusters, clusterSize, refIdx, shift, refPos, size);
	free(clusters);
	free(clusterSize);
	return numHits;
}


/**
 * Finds the biggest clusters.
 *
 * @param 		numBestHits			Number of biggest clusters that are to
 * searched for.
 * @param[out]	refIdx_bestMatch	Best-matching reference index. The number
 * of elements in this array must be equal to @a numBestHits.
 * @param[out]	shift_bestMatch		Best-matching shift. The number of elements
 * in this array must be equal to @a numBestHits.
 * @param[out]	refPos_bestMatch	Best-matching reference position. The number
 * of elements in this array must be equal to @a numBestHits.
 * @param		biggestClusterSize	Biggest cluster size.
 * @param[out]	cluster				Stores cluster indices.
 * @param		numClusters			Number of clusters.
 * @param		clusterSize			Contains cluster sizes.
 * @param		refIdx				Reference Indices.
 * @param		shift				Shifts (reference position - query position).
 * @param		refPos				Reference position.
 * @param		size				Number of elements in @a refIdx.
 * @return 		Number of best-hits in @a refIdx_bestMatch.
 */
static int findBiggestClusters(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch, int biggestClusterSize,
		int *cluster, int numClusters, int *clusterSize, char *refIdx,
		int *shift, int *refPos, int size)
{
	int i, j = 0;
	int *biggestCluster = (int *) calloc(size, sizeof(int));
	for (i = 0; i < numClusters; ++i)
	{
		if (clusterSize[i] == biggestClusterSize)
		{
			biggestCluster[j] = cluster[i];
			++j;
		}
	}
	int numBiggestClusters = j;

	if (numBestHits > numBiggestClusters)
		numBestHits = numBiggestClusters;

	int *randNumArr = (int *) calloc(size, sizeof(int));
	arrGetRandomNums((uint) numBestHits, 0, (uint) (numBiggestClusters - 1),
			(uint *) randNumArr);
	qsort(randNumArr, numBestHits, sizeof(int), compare);
	int tmp;
	for (i = 0; i < numBestHits; ++i)
	{
		tmp = biggestCluster[randNumArr[i]];
		refIdx_bestMatch[i] = refIdx[tmp];
		shift_bestMatch[i] = shift[tmp];
		refPos_bestMatch[i] = refPos[tmp];
	}
	free(randNumArr);
	free(biggestCluster);
	return numBestHits;
}


/**
 * This is a wrapper function that wraps @a findBiggestClusterSize function.
 * It has been added so that @a findBiggestClusterSize can be unit-tested.
 *
 * @param		refIdx				Reference indices.
 * @param		shift				Shift (reference position - query position).
 * @param		refPos				Reference position.
 * @param		size				Number of elements in @a refIdx.
 * @param[out]	cluster				Stores cluster indices.
 * @param[out]	numClusters			Number of clusters.
 * @param[out]	clusterSize			Stores cluster sizes.
 * @param[out]	biggestClusterSize	Biggest cluster size.
 */
void findBiggestClusterSize_wrap6(char *refIdx, int *shift, int *refPos,
		int size, int *cluster, int *numClusters, int *clusterSize,
		int *biggestClusterSize)
{
	createClusters(refIdx, shift, refPos, size, cluster, numClusters);
	findBiggestClusterSize(cluster, *numClusters, clusterSize, biggestClusterSize,
			size);
}


/**
 * Find the biggest cluster size.
 *
 * @param	cluster				Stores indices of clusters.
 * @param	numClusters			Number of clusters.
 * @param	clusterSize			Cluster sizes.
 * @param	biggestClusterSize	Biggest cluster size.
 * @param	size				Number of hits.
 */
static void findBiggestClusterSize(int *cluster, int numClusters,
		int *clusterSize, int *biggestClusterSize, int size)
{
	if (numClusters == 1)
	{
		*biggestClusterSize = size;
		clusterSize[0] = size;
	}
	else
	{
		int i;
		for (i = 1; i < numClusters; ++i)
		{
			clusterSize[i - 1] = cluster[i] - cluster[i - 1];
			if (*biggestClusterSize < clusterSize[i - 1])
				*biggestClusterSize = clusterSize[i - 1];
		}
		clusterSize[i - 1] = size - cluster[i - 1];
		if (*biggestClusterSize < clusterSize[i - 1])
			*biggestClusterSize = clusterSize[i - 1];
	}
}


/**
 * This is a wrapper function that wraps @a createClusters function. It has
 * been added so that @a createClusters can be unit-tested.
 *
 * @param	refIdx	Reference indices.
 * @param	shift	Shifts (reference position - query position).
 * @param	refPos	Reference position.
 * @param	size	Number of elements in @a refIdx.
 * @param[out]	cluster		Stores indices of clusters.
 * @param[out]	numClusters	Number of clusters found.
 */
void createClusters_wrap6(char *refIdx, int *shift, int *refPos, int size,
		int *cluster, int *numClusters)
{
	createClusters(refIdx, shift, refPos, size, cluster, numClusters);
}


/**
 * Creates clusters. It assumes that @a refIdx, @a shift, @a refPos are
 * sorted.
 *
 * @param	refIdx	Reference indices.
 * @param	shift	Shifts (reference position - query position).
 * @param	refPos	Reference position.
 * @param	size	Number of elements in @a refIdx.
 * @param[out]	cluster		Stores indices of clusters.
 * @param[out]	numClusters	Number of clusters found.
 */
static void createClusters(char *refIdx, int *shift, int *refPos, int size,
		int *cluster, int *numClusters)
{
	cluster[0] = 0;
	++(*numClusters);
	int i, j = 1;
	for (i = 1; i < size; ++i)
	{
		if (refIdx[i] == refIdx[i - 1] && shift[i] == shift[i - 1])
		{ }
		else
		{
			cluster[j] = i;
			++j;
		}
	}
	*numClusters = j;
}


/**
 * This is a wrapper function that wraps @a quickSort function. It has
 * been added so that @a quickSort function can be unit-tested.
 *
 * @param[out]	refIdx	Reference index.
 * @param[out]	shift	Shift.
 * @param[out]	refPos	Reference position.
 * @param		size	Number of elements in @a refIdx.
 */
void quickSort_wrap6(char *refIdx, int *shift, int *refPos, int size)
{
	quickSort(0, size - 1, refIdx, shift, refPos);
}


/**
 * Sorts the list of mapped hits using the QuickSort algorithm.
 *
 * @param		left		Beginning index
 * @param		right		Ending index
 * @param[out]	refIdx		Reference indices.
 * @param[out]	shift		Shifts (reference position - query position).
 * @param[out]	refPos		Reference position.
 */
static void quickSort(int left, int right, char *refIdx, int *shift, int *refPos)
{
	int i, last;

	if (left >= right)
		return;
	swap(left, (left + right) / 2, refIdx, shift, refPos);
	last = left;
	for (i = left + 1; i <= right; ++i)
	{
		if ((refIdx[i] < refIdx[left])
				|| (refIdx[i] == refIdx[left] && shift[i] < shift[left])
				|| (refIdx[i] == refIdx[left] && shift[i] == shift[left]
				    && refPos[i] < refPos[left]))
			swap(++last, i, refIdx, shift, refPos);
	}
	swap(left, last, refIdx, shift, refPos);
	quickSort(left, last - 1, refIdx, shift, refPos);
	quickSort(last + 1, right, refIdx, shift, refPos);
}


/**
 * This is a wrapper function that wraps @a swap function. It has been added
 * so that @a swap function can be unit-tested.
 *
 * @param 		i		Index of one of the elements to be swapped.
 * @param 		j		Index of the other element to be swapped.
 * @param[out]	refIdx	Reference index.
 * @param[out]	shift	Shift.
 * @param[out]	refPos	Reference position.
 */
void swap_wrap6(int i, int j, char *refIdx, int *shift, int *refPos)
{
	swap(i, j, refIdx, shift, refPos);
}


/**
 * Swaps element with index @a i with element with index @a j.
 *
 * @param 		i		Index of the first element to be swapped.
 * @param 		j		Index of the second element to be swapped.
 * @param[out]	refIdx	Reference indices.
 * @param[out]	shift	Shifts (reference position - query position).
 * @param[out]	refPos	Reference position.
 */
static void swap(int i, int j, char *refIdx, int *shift, int *refPos)
{
	char tmp = refIdx[i];
	refIdx[i] = refIdx[j];
	refIdx[j] = tmp;

	int tmp2 = shift[i];
	shift[i] = shift[j];
	shift[j] = tmp2;

	tmp2 = refPos[i];
	refPos[i] = refPos[j];
	refPos[j] = tmp2;
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

