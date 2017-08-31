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

#ifndef MAPHITS5_H_
#define MAPHITS5_H_

void createClusters_wrap5(int **cluster);
static void createClusters();
int findBiggestClusterSize_wrap5();
static void findBiggestClusterSize();
int mapHits5GetBestHits(int numBestHits, char *refIdx, int *refPos);
void mapHits5AddHit(char refIdx, int refPos);
void mapHits5Delete();
void mapHits5Create(int seedLen);
void quickSort_wrap5(char **refIdx, int **refPos);
static void quickSort(int left, int right);
void swap_wrap5(int i, int j, char **refIdx, int **refPos);
static void swap(int i, int j);
static int compare(const void *a, const void *b);
int findBiggestClusters_wrap5(int numBestHits, char *refIdx, int *refPos);
static int findBiggestClusters(int numBestHits, char *refIdx, int *refPos);
void mapHits5AddHit_wrap(char **refIdxArr, int **refPosArr, int *size,
		char refIdx, int refPos);
void mapHits5Print();
void mapHits5Reset();

#endif /* MAPHITS5_H_ */
