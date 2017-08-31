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

#ifndef MAPHITS6_H_
#define MAPHITS6_H_

void createClusters_wrap6(char *refIdx, int *shift, int *refPos, int size,
		int *cluster, int *numClusters);

static void createClusters(char *refIdx, int *shift, int *refPos, int size,
		int *cluster, int *numClusters);

void findBiggestClusterSize_wrap6(char *refIdx, int *shift, int *refPos,
		int size, int *cluster, int *numClusters, int *clusterSize,
		int *biggestClusterSize);

static void findBiggestClusterSize(int *cluster, int numClusters,
		int *clusterSize, int *biggestClusterSize, int size);

int mapHits6GetBestHits(int numBestHits, char *refIdx, int *shift, int *refPos,
		int size, char *refIdx_bestMatch, int *shift_bestMatch,
		int *refPos_bestMatch);

void quickSort_wrap6(char *refIdx, int *shift, int *refPos, int size);

static void quickSort(int left, int right, char *refIdx, int *shift, int *refPos);

void swap_wrap6(int i, int j, char *refIdx, int *shift, int *refPos);

static void swap(int i, int j, char *refIdx, int *shift, int *refPos);

static int compare(const void *a, const void *b);

int findBiggestClusters_wrap6(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch, char *refIdx, int *shift,
		int *refPos, int size);

static int findBiggestClusters(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch, int biggestClusterSize,
		int *cluster, int numClusters, int *clusterSize, char *refIdx,
		int *shift, int *refPos, int size);

#endif /* MAPHITS6_H_ */
