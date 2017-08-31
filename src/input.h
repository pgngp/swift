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

/**
 * @file
 *
 * @section DESCRIPTION
 *
 * This file contains functions that fetch sequences from input files.
 */

#ifndef INPUT_H_
#define INPUT_H_

#include <stdio.h>

int getSeqCount(char *fileName);

int getSeqSize(char *fileName);

void getSeq(FILE *filePtr, ulong fileOffset, char *seq, int maxLength);

void getSeq2(FILE *filePtr, ulong fileOffset, char *seq, int maxLength);

void getRefSeq(FILE *filePtr, ulong fileOffset, char *seq, int maxLength);

void getRefSeq2(FILE *filePtr, ulong fileOffset, char *seq, int maxLength);

void getOffsets(FILE *filePtr, uint *queryNameOffset, uint *querySeqOffset,
		uint *refNameOffset, uint *refSeqOffset, uint *refSeqFragOffset,
		int *refDistance, uint *isReverseComplement, int *refIdx);

void getOffsets_paired(FILE *filePtr, ulong *queryNameOffset,
        ulong *querySeqOffset, int *mateIdx, ulong *refNameOffset,
        ulong *refSeqOffset, ulong *refSeqFragOffset, int *refDistance,
        uint *isRevComp, int *refIdx);

void getOffsets_paired2(FILE *filePtr, uint *queryNameOffset,
        uint *querySeqOffset, int *mateIdx, uint *refNameOffset,
        uint *refSeqOffset, uint *refSeqFragOffset, int *refDistance,
        uint *isRevComp, int *refIdx, int *qryIdx);

void getSeqName(FILE *filePtr, ulong fileOffset, char *name);

#endif /* INPUT_H_ */
