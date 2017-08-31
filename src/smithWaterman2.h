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

#ifndef SMITHWATERMAN2_H_
#define SMITHWATERMAN2_H_

__global__ void smithWaterman2(char *refSeq, int refSeqLength, char *querySeq,
		int querySeqLength, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *maxScore, int *startPos, int *endPos,
		int *alignLength, char *refConsensus, char *queryConsensus,
		int numQueries, int numReferences);

__device__ void computeMaxScore2(int querySeqLength, int refSeqLength,
		char *topLeft, char *left, char *top, const char *querySeq,
		const char *refSeq, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *maxScore, int *cellMax);

__device__ void backtrack2(char *queryConsensus, char *refConsensus,
		int *startPos, int *endPos, int *alignLength, char *topLeft, char *left,
		char *top, int querySeqLength, int refSeqLength, int cellMax,
		char *querySeq, char *refSeq);

__device__ int toupper2(int ch);

__device__ int isMatch2(char a, char b);

#endif /* SMITHWATERMAN2_H_ */
