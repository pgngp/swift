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
 * This file contains device functions that implement the Smith-Waterman
 * algorithm.
 */

#ifndef SMITHWATERMAN_H_
#define SMITHWATERMAN_H_


__device__ int toupper(int ch);
int toupper_cpu(int ch);
__global__ void smithWaterman(char *refSeq, int refSeqLength, char *querySeq,
		int querySeqLength, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *maxScore, int *startPos, int *endPos,
		int *alignLength, char *refConsensus, char *queryConsensus,
		int numQueries, int numReferences);
__device__ void computeMaxScore(int querySeqLength, int refSeqLength,
		char *topLeft, char *left, char *top, const char *querySeq,
		const char *refSeq, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *maxScore, int *cellMax);
__device__ void computeMaxScore2(int maxQryLen, int refSeqLength,
		char *topLeft, char *left, char *top, const char *querySeq,
		const char *refSeq, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *maxScore, int *cellMax);
__device__ void backtrack(char *queryConsensus, char *refConsensus,
		int *startPos, int *endPos, int *alignLength, char *topLeft,
		char *left, char *top, int querySeqLength, int refSeqLength,
		int cellMax, char *querySeq, char *refSeq);
__device__ void backtrack2(char *queryConsensus, char *refConsensus,
		int *startPos, int *endPos, int *alignLength, char *topLeft, char *left,
		char *top, int maxQryLen, int refSeqLength, int cellMax, char *querySeq,
		char *refSeq);
__global__ void computeMaxScore_wrapFunc(int qryLen, int refLen, char *topLeft,
		char *left, char *top, const char *qrySeq, const char *refSeq,
		float match, float mismatch, float gapOpenPenalty, float gapExtPenalty,
		float *maxScore, int *cellMax);
__global__ void backtrack_wrapFunc(char *qryConsensus, char *refConsensus,
		int *startPos, int *endPos, int *alignLen, char *topLeft, char *left,
		char *top, int qryLen, int refLen, int cellMax, char *qrySeq,
		char *refSeq);
void smithWaterman_cpu(float *HMatrix, int *rowBacktracker,
		int *colBacktracker, char *refSeq, int refLen, char *qrySeq,
		int qryLen, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *alignScore, int *startPos, int *endPos,
		int *alignLen, char *refConsensus, char *qryConsensus);
void computeMaxScore_cpu_wrapFunc(float *H, int *rowBacktracker,
		int *colBacktracker, char *refSeq, int refLen, char *qrySeq,
		int qryLen, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *alignScore, int *maxRow, int *maxCol);
static void computeMaxScore_cpu(float *H, int *rowBacktracker,
		int *colBacktracker, char *refSeq, int refLen, char *qrySeq,
		int qryLen, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *alignScore, int *maxRow, int *maxCol);
void backtrack_cpu_wrapFunc(int *rowBacktracker, int *colBacktracker,
		int maxRow, int maxCol, char *refSeq, int refLen, char *qrySeq,
		int qryLen, int *startPos, int *endPos, int *alignLen,
		char *refConsensus, char *qryConsensus);
static void backtrack_cpu(int *rowBacktracker, int *colBacktracker,
		int maxRow, int maxCol, char *refSeq, int refLen, char *qrySeq,
		int qryLen, int *startPos, int *endPos, int *alignLen,
		char *refConsensus, char *qryConsensus);
static float findArrMax(float array[], int size, int *index);
__global__ void isMatch_wrap(char *a, char *b, int *retVal);
__device__ int isMatch(char a, char b);
int isMatch_cpu(char a, char b);

#endif /* SMITHWATERMAN_H_ */
