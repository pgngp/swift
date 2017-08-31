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
 * This file contains functions allocate memory on the host and device. It
 * also contains functions that copy memory between host and device.
 */

#ifndef MEMORY_H_
#define MEMORY_H_

void copyMemHostToDevice(
		char *seq_d,
		char *seq,
		int numBytes);

void allocateHostMem(
		float **score,
		int **refDistance,
		int **alignStart,
		int **alignEnd,
		int **alignLength,
		char **queryConsensus,
		char **refConsensus,
		char **queryName,
		char **querySeq,
		char **refName,
		char **refSeq,
		int numQueries,
		int maxQryNameLength,
		int maxRefNameLength,
		int querySeqLength,
		int refSeqLength,
		int combinedLength);

void allocateDeviceMem(
		char **querySeq,
		char **refSeq,
		float **score,
		float **hMatrix,
		int **cellBacktracker,
		int **alignStart,
		int **alignEnd,
		int **alignLength,
		char **queryConsensus,
		char **refConsensus,
		int numQueries,
		int querySeqLength,
		int refSeqLength,
		int combinedLength);

void copyMemDeviceToHost(
		float *score,
		float *score_d,
		int *alignStart,
		int *alignStart_d,
		int *alignEnd,
		int *alignEnd_d,
		int *alignLength,
		int *alignLength_d,
		char *queryConsensus,
		char *queryConsensus_d,
		char *refConsensus,
		char *refConsensus_d,
		int numQueries,
		int numRefs,
		int combinedLength);

#endif /* MEMORY_H_ */
