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

#include "smithWaterman2.h"
#include "common.h"
#include <stdio.h>


/**
 * Perform Smith-Waterman alignment on the given query and reference
 * sequences.
 * @param refSeq Reference sequence.
 * @param refSeqLength Reference sequence length.
 * @param querySeq Query sequence.
 * @param querySeqLength Query sequence length.
 * @param match Match value.
 * @param mismatch Mismatch value.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param maxScore Array where alignment scores will be stored.
 * @param startPos Array where alignment start positions will be stored.
 * @param endPos Array where alignment end positions will be stored.
 * @param alignLength Array where alignment lengths will be stored.
 * @param refConsensus Character pointer where reference consensus sequence
 * will be stored.
 * @param queryConsensus Character pointer where query consensus sequence will
 * be stored.
 * @param numQueries Number of query sequences stored in querySeq.
 * @param numReferences Number of reference sequences stored in refSeq.
 */
__global__ void smithWaterman2(char *refSeq, int refSeqLength, char *querySeq,
		int querySeqLength, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *maxScore, int *startPos, int *endPos,
		int *alignLength, char *refConsensus, char *queryConsensus,
		int numQueries, int numReferences)
{
	__shared__ char topLeft[3801];
	__shared__ char left[3801];
	__shared__ char top[3801];
	if (threadIdx.x == 0)
	{
		int i = 0;
		for (i = 0; i < 3801; ++i)
		{
			topLeft[i] = 0;
			left[i] = 0;
			top[i] = 0;
		}
	}
	__syncthreads();

	/* Compute max score */
	int cellMax = 0;
	computeMaxScore2(querySeqLength, refSeqLength, topLeft, left, top, querySeq,
			refSeq, match, mismatch, gapOpenPenalty, gapExtPenalty,
			maxScore, &cellMax);
	__syncthreads();

	if (threadIdx.x > 0)
		return;

	/* Backtrack */
	backtrack2(queryConsensus, refConsensus, startPos, endPos, alignLength,
			topLeft, left, top, querySeqLength, refSeqLength, cellMax,
			querySeq, refSeq);
}


/**
 * Computes maximum alignment score and creates a matrix that can be used
 * for backtracking to create alignment sequences.
 *
 * @param querySeqLength Query sequence length.
 * @param refSeqLength Reference sequence length.
 * @param topLeft A bit value of '1' indicates that the cell was computed
 * using the top-left cell.
 * @param left A bit value of '1' indicates that the cell was computed
 * using the left cell.
 * @param top A bit value of '1' indicates that the cell was computed
 * using the top cell.
 * @param querySeq Query sequence.
 * @param refSeq Reference sequence.
 * @param match Match value.
 * @param mismatch Mismatch value.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param maxScore Array where alignment scores will be stored.
 * @param cellMax Cell containing the max score.
 */
__device__ void computeMaxScore2(int querySeqLength, int refSeqLength,
		char *topLeft, char *left, char *top, const char *querySeq,
		const char *refSeq, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *maxScore, int *cellMax)
{
	__shared__ char querySeqShared[100];
	__shared__ char refSeqShared[300];
	__shared__ float HShared[3 * (100 + 1)];
	__shared__ float scores[100 + 1];
	__shared__ int cell[100 + 1];

	int tx = threadIdx.x;

	/* Copy query sequence from global memory to shared memory */
	if (tx < querySeqLength)
	{
		int querySeqStartPos = blockIdx.x * (querySeqLength + 1);
		querySeqShared[tx] = toupper2(querySeq[querySeqStartPos + tx]);
	}

	/* Copy reference sequence from global memory to shared memory */
	int refSeqLength2 = refSeqLength + 1;
	if (tx < refSeqLength)
	{
		int refSeqStartPos = blockIdx.x * (refSeqLength + 1);
		refSeqShared[tx] = toupper2(refSeq[refSeqStartPos + tx]);
	}

	int querySeqLength2 = querySeqLength + 1;
	if (tx < querySeqLength2)
	{
		HShared[(0 * querySeqLength2) + tx] = 0.0f;
		HShared[(1 * querySeqLength2) + tx] = 0.0f;
		HShared[(2 * querySeqLength2) + tx] = 0.0f;
		scores[tx] = 0.0f;
		cell[tx] = 0;
	}
	__syncthreads();


	int numIterations = refSeqLength + querySeqLength + 1;
	int currentHOffset, prevHOffset, prevPrevHOffset;
	float topLeftCellScore, topCellScore, leftCellScore, currentCellScore;
	int hLength = 3 * querySeqLength2;
	prevPrevHOffset = 0 * querySeqLength2;
	prevHOffset = 1 * querySeqLength2;
	currentHOffset = 2 * querySeqLength2;
	int refIndex, queryIndex, i, currentCellIndex, topCellIndex, leftCellIndex;
	int byte = 0, bit = 0;
	int charSize = (int) (sizeof(char) * NUM_BITS_IN_A_BYTE);

	for (i = 2; i < numIterations; ++i)
	{
		refIndex = i - tx;
		queryIndex = tx;
		if (refIndex > 0 && refIndex <= refSeqLength && queryIndex > 0
				&& queryIndex <= querySeqLength)
		{
			if (isMatch2(refSeqShared[refIndex - 1],
					querySeqShared[queryIndex - 1]) == 0)
				topLeftCellScore = HShared[prevPrevHOffset + tx - 1] + match;
			else
				topLeftCellScore = HShared[prevPrevHOffset + tx - 1] + mismatch;

			currentCellIndex = ((refSeqLength + 1) * tx) + refIndex;
			topCellIndex = currentCellIndex - refSeqLength2;
			byte = topCellIndex / charSize;
			bit = charSize - 1 - (topCellIndex % charSize);
			if ((top[byte] & (1 << bit)) > 0)
				topCellScore = HShared[prevHOffset + tx - 1] + gapExtPenalty;
			else
				topCellScore = HShared[prevHOffset + tx - 1] + gapOpenPenalty;

			leftCellIndex = currentCellIndex - 1;
			byte = leftCellIndex / charSize;
			bit = charSize - 1 - (leftCellIndex % charSize);
			if ((left[byte] & (1 << bit)) > 0)
				leftCellScore = HShared[prevHOffset + tx] + gapExtPenalty;
			else
				leftCellScore = HShared[prevHOffset + tx] + gapOpenPenalty;

			currentCellScore = 0.0f;
			byte = currentCellIndex / charSize;
			bit = charSize - 1 - (currentCellIndex % charSize);
			if (currentCellScore < topLeftCellScore)
			{
				currentCellScore = topLeftCellScore;
				topLeft[byte] |= (1 << bit); /* Set the bit in 'topLeft'. */
			}
			if (currentCellScore < topCellScore)
			{
				currentCellScore = topCellScore;
				top[byte] |= (1 << bit); /* Set the bit in 'top'. */
				topLeft[byte] &= ~(1 << bit); /* Clear the bit in 'topLeft'. */
			}
			if (currentCellScore < leftCellScore)
			{
				currentCellScore = leftCellScore;
				left[byte] |= (1 << bit); /* Set the bit in 'left'. */
				topLeft[byte] &= ~(1 << bit); /* Clear the bit in 'topLeft'. */
				top[byte] &= ~(1 << bit); /* Clear the bit in 'top'. */
			}
			HShared[currentHOffset + tx] = currentCellScore;

			if (scores[tx] < currentCellScore)
			{
				scores[tx] = currentCellScore;
				cell[tx] = currentCellIndex;
			}
		}
		__syncthreads();

		prevPrevHOffset = prevHOffset;
		prevHOffset = currentHOffset;
		currentHOffset = currentHOffset + querySeqLength2;
		if (currentHOffset >= hLength)
			currentHOffset = 0;
	}


	if (tx == 0)
	{
		float scoreMax = 0.0f;
		*cellMax = 0;
		for (i = 0; i < querySeqLength2; ++i)
		{
			if (scoreMax < scores[i])
			{
				scoreMax = scores[i];
				*cellMax = cell[i];
			}
		}
		maxScore[blockIdx.x] = scoreMax;
	}
}


/**
 * Backtracks the alignment score computation matrix to create the
 * alignment string.
 *
 * @param queryConsensus Character pointer where the query consensus sequence
 * will be stored.
 * @param refConsensus Character pointer where reference consensus sequence
 * will be stored.
 * @param startPos Array where alignment start positions will be stored.
 * @param endPos Array where alignment end positions will be stored.
 * @param alignLength Array where alignment length will be stored.
 * @param topLeft A bit value of '1' indicates that the cell was computed
 * using the top-left cell.
 * @param left A bit value of '1' indicates that the cell was computed
 * using the left cell.
 * @param top A bit value of '1' indicates that the cell was computed
 * using the top cell.
 * @param querySeqLength Length of query sequence.
 * @param refSeqLength Length of reference sequence.
 * @param cellMax Cell containing the max score.
 * @param querySeq Query sequence.
 * @param refSeq Reference sequence.
 */
__device__ void backtrack2(char *queryConsensus, char *refConsensus,
		int *startPos, int *endPos, int *alignLength, char *topLeft, char *left,
		char *top, int querySeqLength, int refSeqLength, int cellMax,
		char *querySeq, char *refSeq)
{
	startPos[blockIdx.x] = 0;
	endPos[blockIdx.x] = 0;
	alignLength[blockIdx.x] = 0;

	int currentCellIndex = cellMax;
	int currentRow = currentCellIndex / (refSeqLength + 1);
	int currentCol = currentCellIndex % (refSeqLength + 1);
	int charSize = (int) (sizeof(char) * NUM_BITS_IN_A_BYTE);
	int byte = currentCellIndex / charSize;
	int bit = charSize - 1 - (currentCellIndex % charSize);
	int nextCell;
	if ((top[byte] & (1 << bit)) > 0)
		nextCell = currentCellIndex - refSeqLength - 1;
	else if ((left[byte] & (1 << bit)) > 0)
		nextCell = currentCellIndex - 1;
	else if ((topLeft[byte] & (1 << bit)) > 0)
		nextCell = currentCellIndex - refSeqLength - 2;
	int nextRow = nextCell / (refSeqLength + 1);
	int nextCol = nextCell % (refSeqLength + 1);
	int tick = -1;
	int colIncrementer = -1;
	int numRefBases = -1;

	refSeq += (blockIdx.x * (refSeqLength + 1));
	querySeq += (blockIdx.x * (querySeqLength + 1));
	refConsensus += (blockIdx.x * (refSeqLength + querySeqLength + 2));
	queryConsensus += (blockIdx.x * (refSeqLength + querySeqLength + 2));

	while ((currentRow != nextRow || currentCol != nextCol)
			&& nextCol != 0 && nextRow != 0)
	{
		++tick;
		++colIncrementer;

		if (nextCol == currentCol)
		{
			/* Deletion in reference */
			refConsensus[tick] = '-';
		}
		else
		{
			/* Match/mismatch in reference */
			refConsensus[tick] = refSeq[currentCol - 1];
			++numRefBases;
		}

		if (nextRow == currentRow)
		{
			/* Deletion in query */
			queryConsensus[tick] = '-';
		}
		else
		{
			/* Match/mismatch in query */
			queryConsensus[tick] = querySeq[currentRow - 1];
		}

		startPos[blockIdx.x] = currentCol;
		currentRow = nextRow;
		currentCol = nextCol;
		currentCellIndex = nextCell;
		byte = currentCellIndex / charSize;
		bit = charSize - 1 - (currentCellIndex % charSize);
		if ((top[byte] & (1 << bit)) > 0)
			nextCell = currentCellIndex - refSeqLength - 1;
		else if ((left[byte] & (1 << bit)) > 0)
			nextCell = currentCellIndex - 1;
		else if ((topLeft[byte] & (1 << bit)) > 0)
			nextCell = currentCellIndex - refSeqLength - 2;
		else
		{
			currentRow = nextRow;
			currentCol = nextCol;
			break;
		}
		nextRow = nextCell / (refSeqLength + 1);
		nextCol = nextCell % (refSeqLength + 1);
	}
	if ((currentRow != nextRow || currentCol != nextCol)
			&& (nextCol == 0 || nextRow == 0))
	{
		++tick;
		++colIncrementer;
		if (nextCol == currentCol)
		{
			/* Deletion in reference */
			refConsensus[tick] = '-';
		}
		else
		{
			/* Match/mismatch in reference */
			refConsensus[tick] = refSeq[currentCol - 1];
			++numRefBases;
		}

		if (nextRow == currentRow)
		{
			/* Deletion in query */
			queryConsensus[tick] = '-';
		}
		else
		{
			/* Match/mismatch in query */
			queryConsensus[tick] = querySeq[currentRow - 1];
		}

		startPos[blockIdx.x] = currentCol;
	}
	endPos[blockIdx.x] = startPos[blockIdx.x] + numRefBases;
	alignLength[blockIdx.x] = colIncrementer + 1;
}


/**
 * Returns the upper-case of the given character.
 *
 * @param ch Character which needs to be converted to upper-case.
 * @return Upper-case letter.
 *
 * Code from koders.com.
 */
__device__ int toupper2(int ch)
{
	if ((unsigned int)(ch - 'a') < 26u)
		ch += 'A' - 'a';
	return ch;
}


/**
 * Determines whether base @a is equal to base @a b. If they are, returns
 * 0; otherwise returns -1.
 *
 * @note Base 'N' matches any base, including 'N'.
 *
 * @param	a	First base.
 * @param	b	Second base.
 * @returns	0 if both bases are equal and -1 if they are not.
 */
__device__ int isMatch2(char a, char b)
{
	if (a == b || a == 'N' || b == 'N')
		return 0;
	else
		return -1;
}

