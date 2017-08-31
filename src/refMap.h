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

#ifndef REFMAP_H_
#define REFMAP_H_

void refMapCreate_wrap(const char *file, uint **nameOffset, uint **seqOffset,
		int **numBasesPerLine, int *size);

void refMapCreate(const char *file);

uint refMapGetFileOffset(uint refNameOffset, uint pos);

void refMapGetOffsets(int index, uint pos, uint *nameOffset, uint *seqOffset);

void refMapGetIndexAndPos(uint seqOffset, int *index, int *pos);

void refMapFree();

int binarySearch_wrapper(uint x, uint *arr, uint arrSize);

static int binarySearch(uint x, uint *arr, uint arrSize);

void refMapPrint2();

#endif /* REFMAP_H_ */
