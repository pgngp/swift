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

#ifndef ALIGN_H_
#define ALIGN_H_

int getMemPerRefQryPair(int refLength, int queryLength);

void alignQueries(char *queryFileName, int querySeqLength, char *refFileName,
		char *matchFile, int outFormat, float match, float mismatch,
		float gapOpenPenalty, float gapExtPenalty, const char *outputFile,
		int numMatches);

void alignQueries2(char *queryFileName, int querySeqLength, char *refFileName,
		char *matchFile, int outFormat, float match, float mismatch,
		float gapOpenPenalty, float gapExtPenalty, const char *outputFile,
		int numMatches);

void alignQueries3(char *queryFileName, int maxQryLen, char *refFileName,
		char *matchFile, int outFormat, float match, float mismatch,
		float gapOpenPenalty, float gapExtPenalty, const char *outputFile,
		int numMatches);

void alignQueries_paired(char *qryFile1, char *qryFile2, int numQrs,
		int qryLen, char *refFile, int numRefs, int refLen, char *matchFile,
		int outFormat, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, const char *outputFile, uint minFragSize,
		uint maxFragSize);

void alignQueries_paired2(char *qry1File, char *qry2File, int qryLen,
		char *refFile, char *matchFile, int outFormat, float match,
		float mismatch, float gapOpenPenalty, float gapExtPenalty,
		const char *outputFile, uint minFragSize, uint maxFragSize,
		int numMatches);

void alignQueries_cpu(char *qryFile, int numQrs, int qryLen, char *refFile,
		char *matchFile, int outFormat, float match, float mismatch,
		float gapOpenPenalty, float gapExtPenalty, const char *outputFile);

void alignQueries_cpu_mt(char *qryFile, int numQrs, int qryLen, char *refFile,
		char *matchFile, int outFormat, float match, float mismatch,
		float gapOpenPenalty, float gapExtPenalty, const char *outputFile);

void alignQueries_paired_cpu(char *qry1File, char *qry2File, int numQrs,
		int qryLen, char *refFile, int numRefs, int refLen, char *matchFile,
		int outFormat, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, const char *outputFile);

int getNumQueriesPerIteration(cudaDeviceProp *deviceProp,
		int memPerRefQueryPair, int numTotalQueries);

int getNumQueriesPerIteration_paired(cudaDeviceProp *deviceProp,
		int memPerRefQueryPair, int numTotalQueries);

int getNumSeqsInCurrentIteration(int iterationIndex, int numTotalSeqs,
		int numIterations, int maxSeqsPerIteration);

int getNumSeqsInLastIteration(int numTotalSeqs, int numIterations,
		int maxSeqsPerIteration);

void alignSequences(char *querySeq_d, char *refSeq_d, int querySeqLength,
		int refSeqLength, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *score_d, int *alignStart_d,
		int *alignEnd_d, int *alignLength_d, char *queryConsensus_d,
		char *refConsensus_d, int queriesPerIteration, int refsPerIteration,
		int maxThreadsPerBlock);

static int getNumMatches(const char *filteredResultsFile);

void chooseBestAlignment_wrap(char *qryName, int qryNameLen, float *score,
		int size);

static void chooseBestAlignment(char *qryName, int qryNAmeLen, float *score,
		int size);

static void chooseBestAlignment2(char *qryName, int qryNameLen, float *score,
		int size);

void chooseBestAlignment_paired_wrap(char *qryName, int qryNameLen,
		float *score, int size);

static void chooseBestAlignment_paired(char *qryName, int qryNameLen,
		float *score, int size);

void chooseBestAlignment_paired2_wrap(char *qryName, int qryNameLen,
		float *score, int *alignStart, int *alignEnd, int size,
		uint minFragSize, uint maxFragSize);

static void chooseBestAlignment_paired2(char *qryName, int qryNameLen,
		float *score, int *alignStart, int *alignEnd, int size,
		uint minFragSize, uint maxFragSize);

void chooseBestAlignment_paired3_wrap(char *qryName, int qryNameLen,
		uint *revComp, char *refName, int refNameLen, float *score,
		int *refDistance, int *alignStart, int *alignEnd, int size,
		uint minFragSize, uint maxFragSize);

static void chooseBestAlignment_paired3(char *qryName, int qryNameLen,
		uint *revComp, char *refName, int refNameLen, float *score,
		int *refDistance, int *alignStart, int *alignEnd, int size,
		uint minFragSize, uint maxFragSize);

int qryNameCmp_wrap(char *s, char *t, int size);

static int qryNameCmp(char *s, char *t, int size);

void getMatePairName(const char *qryName, char *mateName);

void getNumQrsPerIter_paired_wrap(const char *file, int maxNumQrsPerIter,
		int **numQrsPerIter, int *arrSize);

static void getNumQrsPerIter_paired(const char *file, int maxNumQrsPerIter,
		int **numQrsPerIter, int *arrSize);

void *getAlignments(void *arg);

#endif /* ALIGN_H_ */
