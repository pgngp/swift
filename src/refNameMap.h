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

#ifndef REFNAMEMAP_H_
#define REFNAMEMAP_H_

void refNameMapCreate(const char *refFile);

void refNameMapCreate_wrap(const char *refFile, int **namePos,
		int **seqStartPos, int **seqEndPos, int *size);

void refNameMapDelete();

int refNameMapGetNameOffset(uint seqOffset);

void refNameMapPrint();

#endif /* REFNAMEMAP_H_ */
