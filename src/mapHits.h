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

#ifndef MAPHITS_H_
#define MAPHITS_H_

void createClusters_wrap(int **cluster);
static void createClusters();
int findBiggestClusterSize_wrap();
static void findBiggestClusterSize();
int mapHitsGetBestHits(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch);
void mapHitsAddHit(char refIdx, int shift, int refPos);
void mapHitsDelete();
void mapHitsCreate(int seedLen);
void quickSort_wrap(char **refIdx, int **shift, int **refPos);
static void quickSort(int left, int right);
void swap_wrap(int i, int j, char **refIdx, int **shift, int **refPos);
static void swap(int i, int j);
static int compare(const void *a, const void *b);
int findBiggestClusters_wrap(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch);
static int findBiggestClusters(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch);
void mapHitsAddHit_wrap(char **refIdxArr, int **shiftArr, int **refPosArr,
		int *size, char refIdx, int shift, int refPos);
void mapHitsPrint();
void mapHitsReset();

#endif /* MAPHITS_H_ */
