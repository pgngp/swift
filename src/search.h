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

#ifndef SEARCH_H_
#define SEARCH_H_

#include "query.h"
#include "hitList.h"
#include <stdio.h>

void searchQueries(char *queryFileName, int numQueries, int querySeqLength,
		char *refFileName, int seedLength, char *matchFile, uint maxNumHits,
		int *numMatches);

void searchQueries2(char *qryFileName, int numQrs, int qryLen,
		char *refFileName, int seedLen, char *matchFile, uint maxNumHits,
		int *numMatches);

void searchQueries3(char *qryFileName, int numQrs, int qryLen,
		char *refFileName, int seedLen, char *matchFile, uint maxNumHits,
		int *numMatches, int numThreads, int tupleIgnoreThres);

void searchQueries4(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches);

void searchQueries4_2(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches);

void searchQueries4_3(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches, int numThreads,
		int tupleIgnoreThres);

void searchQueries5(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches, int tupleIgnoreThres);

void searchQueries9(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches, int tupleIgnoreThres);

void searchQueries10(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches, int tupleIgnoreThres);

void searchQueries11(char *qryFileName, char *refFileName, char *matchFile,
		int seedLen, int maxNumHits, int *numMatches, int tupleIgnoreThres);

void searchQueries_gpu(char *qryFileName, int numQueries, int qryLength,
		char *refFileName, int numReferences, int refLength, int seedLength,
		char *matchFile, uint maxNumHits);

void searchQueries_paired(char *qryFile1, char *qryFile2, int numQrs,
		int qryLen, char *refFile, int numRefs, int refLen, int seedLen,
		char *matchFile, uint maxNumHits, uint minDistance, uint maxDistance);

void searchQueries_paired2(char *qryFile1, char *qryFile2, int numQrs,
		int qryLen, char *refFile, int seedLen, char *matchFile,
		uint maxNumHits, int *numMatches);

static void printHeuristicMatches(Query *qry, int isReverseComplement,
		HitList *matches, FILE *filePtr, int seedLen, int *numMatches);

static void printHeuristicMatches2(Query **qry, char *refIdxArr, int *shiftArr,
		int *posArr, long numQrsPerIter, int maxNumHits, FILE *filePtr,
		int seedLen);

static void printHeuristicMatches3(Query *qry, int isRevComp, int mateIdx,
		HitList *matches, FILE *filePtr, int seedLen, int *numMatches);

static void printHeuristicMatches4(Query *qry, int isRevComp, char *refIdx,
		int *shift, int *refPos, int numMatches, FILE *filePtr);

static void printHeuristicMatches5(Query **qry, char *refIdx, int *refPos,
		int numResults, int *numMatches, FILE *filePtr);

static void printHeuristicMatches5_2(Query **qry, char *refIdx, int *refPos,
		int numQrs, int maxHits, int *numMatches, FILE *filePtr);

static void printHeuristicMatches6(Query *qry, int strand, char *refIdx,
		int *refPos, int maxHitsPerQry, FILE *filePtr);

static void printHeuristicMatches10(Query *qry, int strand, char *refIdx,
		int *refPos, int *refSize, int numHits, int *numMatches, FILE *filePtr);

static void printHeuristicMatches11(Query *qry, int strand, char *refIdx,
		int *refPos, int numHits, int *numMatches, FILE *filePtr);

static void printHeuristicMatches_paired(const Query *qry1, const Query *qry2,
		HitList **matches, FILE *filePtr);

void getNumIters_wrap(int *numIter, long *numQrsPerIter, long *numQrsLastIter,
		int numQrs, int maxNumHits, long gpuTotalMem, int gpuMemMargin,
		int keysTotalMem, int valsTotalMem, int qrySeqMem, int refIdxMem,
		int shiftMem, int posMem);

static void getNumIters(int *numIter, long *numQrsPerIter, long *numQrsLastIter,
		int numQrs, int maxNumHits, long gpuTotalMem, int gpuMemMargin,
		int keysTotalMem, int valsTotalMem, int qrySeqMem, int refIdxMem,
		int shiftMem, int posMem);

void getNumBlocks_wrap(int *blocksX, int *blocksY, long availBlocksX,
		long availBlocksY, long numQrsPerIter);

static void getNumBlocks(int *blocksX, int *blocksY, long availBlocksX,
		long availBlocksY, long numQrsPerIter);

void *getHits(void *arg);

void *getHits2(void *arg);

#endif /* SEARCH_H_ */
