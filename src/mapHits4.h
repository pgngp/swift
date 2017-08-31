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

#ifndef MAPHITS4_H_
#define MAPHITS4_H_

void mapHits4Create(int seedLen);

void mapHits4Delete();

void mapHits4Reset();

void mapHits4AddHit_wrap(ulong **hits, int *size, char refIdx, int shift,
		int refPos);

void mapHits4AddHit(char refIdx, int shift, int refPos);

int mapHits4GetBestHits(int numBestHits, char *refIdx, int *shift, int *refPos);

void createClusters_wrap4(int **cluster);

static void createClusters();

int findBiggestClusterSize_wrap4();

static void findBiggestClusterSize();

int findBiggestClusters_wrap4(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch);

static int findBiggestClusters(int numBestHits, char *refIdx_bestMatch,
		int *shift_bestMatch, int *refPos_bestMatch);

static int compare_ulong(const void *a, const void *b);

static int compare_int(const void *a, const void *b);

#endif /* MAPHITS4_H_ */
