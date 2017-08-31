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

#ifndef MAPHITS2_H_
#define MAPHITS2_H_

#include "common.h"

static int compare(const void *a, const void *b);

int findBiggestClusters_wrap2(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch);

static int findBiggestClusters(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch);

int findBiggestClusters2_wrap2(int numBestHits, char *refIdx_bestMatch,
		int *refPos_bestMatch);

static int findBiggestClusters2(int numBestHits, char *refIdx_bestMatch,
		int *refPos_bestMatch);

int findBiggestClusterSize_wrap2();

static void findBiggestClusterSize();

int mapHits2GetBestHits(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch);

int mapHits2GetBestHits2(int numBestHits, char *refIdx_bestMatch,
		int *refPos_bestMatch);

void mapHits2AddHit_wrap(char **refIdxArr, int **shiftArr, int **refPosArr,
		int **numClusterHits, int *size, char refIdx, int shift, int refPos);

void mapHits2AddHit(char refIdx, int shift, int refPos);

void mapHits2AddHit2_wrap(char **ref, uchar **hits, int **pos,
		int *biggestClusterSize, char refIdx, int qryPos, int refPos);

void mapHits2AddHit2(char refIdx, int qryPos, int refPos);

void mapHits2Create(int seedLen);

void mapHits2Create2(const char *refFile, int seedLen);

void mapHits2Reset();

void mapHits2Delete();

void mapHits2Delete2();

int search_wrap(char refIdx, int shift, int *nextHighestIdx);

static int search(char refIdx, int shift, int *nextHighestIdx);

void createRefMap_wrap(int **lastPos);

static void createRefMap(const char *refFile);

#endif /* MAPHITS2_H_ */
