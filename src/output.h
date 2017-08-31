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
 * This file contains functions that print the output in different formats.
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

void printOutput(int numQueries, char *queryName, char *refName, float *score,
		int *refDistance, int *alignStart, int *alignEnd, int *alignLength,
		int combinedLength, char *queryConsensus, char *refConsensus,
		FILE *outputFilePtr);

void printOutput_paired(int numQueries, char *queryName, char *refName,
		float *score, int *refDistance, int *alignStart, int *alignEnd,
		int *alignLength, int combinedLength, char *queryConsensus,
		char *refConsensus, FILE *outputFilePtr);

void printOutput_paired2(int numQueries, int *qryIdx, char *queryName,
        char *refName, float *score, int *refDistance, int *alignStart,
        int *alignEnd, int *alignLength, int combinedLength,
        char *queryConsensus, char *refConsensus, FILE *outputFilePtr);

void printOutput_paired_cpu(int numQueries, char *queryName, char *refName,
		float *score, int *refDistance, int *alignStart, int *alignEnd,
		int *alignLength, int combinedLength, char *queryConsensus,
		char *refConsensus, FILE *outputFilePtr, int isFirstQry);

void printOutputSam(int numQueries, char *queryName, char *refName,
		float *score, int *alignStart, int *alignEnd, int *alignLength,
		int combinedLength, char *queryConsensus, char *refConsensus,
		FILE *outputFilePtr);

void printCigar(char *queryConsensus, char *refConsensus, int alignLength,
		FILE *outputFilePtr);

#endif /* OUTPUT_H_ */
