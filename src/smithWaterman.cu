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

#include "smithWaterman.h"
#include "common.h"
#include <stdio.h>


/**
 * Performs Smith-Waterman alignment on a CPU.
 *
 * @param HMatrix Computation matrix.
 * @param rowBacktracker Matrix used for backtracking rows.
 * @param colBacktracker Matrix used for backtracking columns.
 * @param refSeq Reference sequence.
 * @param refLen Reference sequence length.
 * @param qrySeq Query sequence.
 * @param qryLen Query sequence length.
 * @param match Match score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap open penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param[out] alignScore Alignment score.
 * @param[out] Alignment start position.
 * @param[out] Alignment end position.
 * @param[out] Alignment length.
 * @param[out] Reference consensus sequence.
 * @param[out] Query consensus sequence.
 */
void smithWaterman_cpu(float *HMatrix, int *rowBacktracker,
		int *colBacktracker, char *refSeq, int refLen, char *qrySeq,
		int qryLen, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *alignScore, int *startPos, int *endPos,
		int *alignLen, char *refConsensus, char *qryConsensus)
{
	int maxRow = 0, maxCol = 0;
	computeMaxScore_cpu(HMatrix, rowBacktracker, colBacktracker, refSeq,
			refLen, qrySeq, qryLen, match, mismatch, gapOpenPenalty,
			gapExtPenalty, alignScore, &maxRow, &maxCol);

	backtrack_cpu(rowBacktracker, colBacktracker, maxRow, maxCol, refSeq,
			refLen, qrySeq, qryLen, startPos, endPos, alignLen, refConsensus,
			qryConsensus);
}


/**
 * This is a wrapper function for @a computeMaxScore_cpu. This function
 * has been added so that @a computeMaxScore_cpu can be unit-tested.
 *
 * @param H Computation matrix.
 * @param rowBacktracker Matrix used for backtracking rows.
 * @param colBacktracker Matrix used for backtracking columns.
 * @param refSeq Reference sequence.
 * @param refLen Reference sequence length.
 * @param qrySeq Query sequence.
 * @param qryLen Query sequence length.
 * @param match Match Score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extention penalty.
 * @param[out] alignScore Alignment score.
 * @param[out] maxRow Row containing the max value.
 * @param[out] maxCol Column containing the max value.
 */
void computeMaxScore_cpu_wrapFunc(float *H, int *rowBacktracker,
		int *colBacktracker, char *refSeq, int refLen, char *qrySeq,
		int qryLen, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *alignScore, int *maxRow, int *maxCol)
{
	int i;
	int size = (refLen + 1) * (qryLen + 1);
	for (i = 0; i < size; ++i)
		H[i] = 0.0f;

	computeMaxScore_cpu(H, rowBacktracker, colBacktracker, refSeq, refLen,
			qrySeq, qryLen, match, mismatch, gapOpenPenalty, gapExtPenalty,
			alignScore, maxRow, maxCol);
}


/**
 * Computes alignment score.
 *
 * @param H Computation matrix.
 * @param rowBacktracker Matrix used for backtracking rows.
 * @param colBacktracker Matrix used for backtracking columns.
 * @param refSeq Reference sequence.
 * @param refLen Reference sequence length.
 * @param qrySeq Query sequence.
 * @param qryLen Query sequence length.
 * @param match Match Score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extention penalty.
 * @param[out] alignScore Alignment score.
 * @param[out] maxRow Row containing the max value.
 * @param[out] maxCol Column containing the max value.
 */
static void computeMaxScore_cpu(float *H, int *rowBacktracker,
		int *colBacktracker, char *refSeq, int refLen, char *qrySeq,
		int qryLen, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *alignScore, int *maxRow, int *maxCol)
{
	int row, col, index;
	float tmpScoreArr[4];
	tmpScoreArr[3] = 0.0f;
	int refLen2 = refLen + 1;
	int topLeftCellIdx, topCellIdx, leftCellIdx, currentCellIdx;
	*alignScore = 0;

	for (row = 1; row <= qryLen; ++row)
	{
		for (col = 1; col <= refLen; ++col)
		{
			topLeftCellIdx = ((row - 1) * refLen2) + (col - 1);
			topCellIdx = topLeftCellIdx + 1;
			leftCellIdx = (row * refLen2) + (col - 1);
			currentCellIdx = (row * refLen2) + col;

			if (refSeq[col - 1] == '\0')
				break;
			else if (isMatch_cpu(qrySeq[row - 1], toupper_cpu(refSeq[col - 1])) == 0)
				tmpScoreArr[0] = H[topLeftCellIdx] + match;
			else
				tmpScoreArr[0] = H[topLeftCellIdx] + mismatch;

			if (rowBacktracker[topCellIdx] == (row - 2)
					&& colBacktracker[topCellIdx] == col)
				tmpScoreArr[1] = H[topCellIdx] + gapExtPenalty;
			else
				tmpScoreArr[1] = H[topCellIdx] + gapOpenPenalty;

			if (rowBacktracker[leftCellIdx] == row
					&& colBacktracker[leftCellIdx] == (col - 2))
				tmpScoreArr[2] = H[leftCellIdx] + gapExtPenalty;
			else
				tmpScoreArr[2] = H[leftCellIdx] + gapOpenPenalty;

			H[(row * refLen2) + col] = findArrMax(tmpScoreArr, 4, &index);
			switch (index)
			{
			case 0: /* Score stems from a match/mismatch. */
				rowBacktracker[currentCellIdx] = row - 1;
				colBacktracker[currentCellIdx] = col - 1;
				break;
			case 1: /* Score stems from a deletion in reference. */
				rowBacktracker[currentCellIdx] = row - 1;
				colBacktracker[currentCellIdx] = col;
				break;
			case 2: /* Score stems from a deletion in query. */
				rowBacktracker[currentCellIdx] = row;
				colBacktracker[currentCellIdx] = col - 1;
				break;
			case 3: /* (row, col) is the beginning of a sequence. */
				rowBacktracker[currentCellIdx] = row;
				colBacktracker[currentCellIdx] = col;
				break;
			}

			if (*alignScore < H[(row * refLen2) + col])
			{
				*alignScore = H[(row * refLen2) + col];
				*maxRow = row;
				*maxCol = col;
			}
		}
	}
}


/**
 * Finds max value in the given array.
 *
 * @param array Array in which the max value is to be searched.
 * @param size Number of elements in @a array.
 * @param[out] index Index of the array element containing the max value.
 * @return Max value in the array.
 */
static float findArrMax(float array[], int size, int *index)
{
	float max;
	int i;
	max = -1.0f;

	for (i = 0; i < size; i++)
	{
		if (max < array[i])
		{
			max = array[i];
			*index = i;
		}
	}

	return max;
}


/**
 * This is a wrapper function for @a backtrack_cpu function. This function
 * has been added so that @a backtrack_cpu function can be unit-tested.
 *
 * @param rowBacktracker Row backtracking matrix used to create the alignment
 * strings.
 * @param colBacktracker Column backtracking matrix used to create the
 * alignment strings.
 * @param maxRow Row in @a backtracker matrix that contains the max value.
 * @param maxCol Col in @a backtracker matrix that contains the max value.
 * @param refSeq Reference sequence.
 * @param refLen Reference sequence length.
 * @param qrySeq Query sequence.
 * @param qryLen Query sequence length.
 * @param startPos Alignment start position.
 * @param endPos Alignment end position.
 * @param alignLen Alignment length.
 * @param refConsensus Reference consensus sequence.
 * @param qryConsensus Query consensus sequence.
 */
void backtrack_cpu_wrapFunc(int *rowBacktracker, int *colBacktracker,
		int maxRow, int maxCol, char *refSeq, int refLen, char *qrySeq,
		int qryLen, int *startPos, int *endPos, int *alignLen,
		char *refConsensus, char *qryConsensus)
{
	backtrack_cpu(rowBacktracker, colBacktracker, maxRow, maxCol, refSeq,
			refLen, qrySeq, qryLen, startPos, endPos, alignLen, refConsensus,
			qryConsensus);
}


/**
 * Backtracks the computation matrix to create an alignment string.
 *
 * @param rowBacktracker Row backtracking matrix used to create the alignment
 * strings.
 * @param colBacktracker Column backtracking matrix used to create the
 * alignment strings.
 * @param maxRow Row in @a backtracker matrix that contains the max value.
 * @param maxCol Col in @a backtracker matrix that contains the max value.
 * @param refSeq Reference sequence.
 * @param refLen Reference sequence length.
 * @param qrySeq Query sequence.
 * @param qryLen Query sequence length.
 * @param startPos Alignment start position.
 * @param endPos Alignment end position.
 * @param alignLen Alignment length.
 * @param refConsensus Reference consensus sequence.
 * @param qryConsensus Query consensus sequence.
 */
static void backtrack_cpu(int *rowBacktracker, int *colBacktracker,
		int maxRow, int maxCol, char *refSeq, int refLen, char *qrySeq,
		int qryLen, int *startPos, int *endPos, int *alignLen,
		char *refConsensus, char *qryConsensus)
{
	int refLen2 = refLen + 1;
	int currentRow = maxRow, currentCol = maxCol;
	int nextRow = rowBacktracker[(currentRow * refLen2) + currentCol];
	int nextCol = colBacktracker[(currentRow * refLen2) + currentCol];
	int tick = -1, colIncrementer = -1, numRefBases = -1;

	while ((currentRow != nextRow || currentCol != nextCol)
			&& nextCol != 0 && nextRow != 0)
	{
		++tick;
		++colIncrementer;

		if (nextCol == currentCol)
			/* Deletion in reference */
			refConsensus[tick] = '-';
		else
		{
			/* Match/mismatch in reference */
			refConsensus[tick] = refSeq[currentCol - 1];
			++numRefBases;
		}

		if (nextRow == currentRow)
			/* Deletion in query */
			qryConsensus[tick] = '-';
		else
			/* Match/mismatch in query */
			qryConsensus[tick] = qrySeq[currentRow - 1];

		*startPos = currentCol;
		currentRow = nextRow;
		currentCol = nextCol;
		nextRow = rowBacktracker[(currentRow * refLen2) + currentCol];
		nextCol = colBacktracker[(currentRow * refLen2) + currentCol];
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
			qryConsensus[tick] = '-';
		}
		else
		{
			/* Match/mismatch in query */
			qryConsensus[tick] = qrySeq[currentRow - 1];
		}
		*startPos = currentCol;
	}
	*endPos = *startPos + numRefBases;
	*alignLen = colIncrementer + 1;
}


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
__global__ void smithWaterman(char *refSeq, int refSeqLength, char *querySeq,
		int querySeqLength, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *maxScore, int *startPos, int *endPos,
		int *alignLength, char *refConsensus, char *queryConsensus,
		int numQueries, int numReferences)
{
	__shared__ char topLeft[3801]; /* (101 * 301) / 8 = 3801. */
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
	computeMaxScore(querySeqLength, refSeqLength, topLeft, left, top, querySeq,
			refSeq, match, mismatch, gapOpenPenalty, gapExtPenalty,
			maxScore, &cellMax);
	__syncthreads();

	if (threadIdx.x > 0)
		return;

	/* Backtrack */
	backtrack(queryConsensus, refConsensus, startPos, endPos, alignLength,
			topLeft, left, top, querySeqLength, refSeqLength, cellMax,
			querySeq, refSeq);
}


///**
// * Perform Smith-Waterman alignment on the given query and reference
// * sequences.
// * @param refSeq Reference sequence.
// * @param refSeqLength Reference sequence length.
// * @param querySeq Query sequence.
// * @param maxQryLen Max size of query sequence.
// * @param match Match value.
// * @param mismatch Mismatch value.
// * @param gapOpenPenalty Gap opening penalty.
// * @param gapExtPenalty Gap extension penalty.
// * @param maxScore Array where alignment scores will be stored.
// * @param startPos Array where alignment start positions will be stored.
// * @param endPos Array where alignment end positions will be stored.
// * @param alignLength Array where alignment lengths will be stored.
// * @param refConsensus Character pointer where reference consensus sequence
// * will be stored.
// * @param queryConsensus Character pointer where query consensus sequence will
// * be stored.
// * @param numQueries Number of query sequences stored in querySeq.
// * @param numReferences Number of reference sequences stored in refSeq.
// */
//__global__ void smithWaterman2(char *refSeq, int refSeqLength, char *querySeq,
//		int maxQryLen, float match, float mismatch, float gapOpenPenalty,
//		float gapExtPenalty, float *maxScore, int *startPos, int *endPos,
//		int *alignLength, char *refConsensus, char *queryConsensus,
//		int numQueries, int numReferences)
//{
//#define	ARR_SIZE	7563 /* (201 * 301) / 8 = 7563. */
//	__shared__ char topLeft[ARR_SIZE];
//	__shared__ char left[ARR_SIZE];
//	__shared__ char top[ARR_SIZE];
//	if (threadIdx.x == 0)
//	{
//		int i = 0;
//		for (i = 0; i < ARR_SIZE; ++i)
//		{
//			topLeft[i] = 0;
//			left[i] = 0;
//			top[i] = 0;
//		}
//	}
//	__syncthreads();
//
//	/* Compute max score */
//	int cellMax = 0;
//	computeMaxScore2(maxQryLen, refSeqLength, topLeft, left, top, querySeq,
//			refSeq, match, mismatch, gapOpenPenalty, gapExtPenalty,
//			maxScore, &cellMax);
//	__syncthreads();
//
//	if (threadIdx.x > 0)
//		return;
//
//	/* Backtrack */
//	backtrack2(queryConsensus, refConsensus, startPos, endPos, alignLength,
//			topLeft, left, top, maxQryLen, refSeqLength, cellMax,
//			querySeq, refSeq);
//}


/**
 * This is a wrapper function that wraps the computeMaxScore function.
 *
 * @param qryLen Query length.
 * @param refLen Reference length.
 * @param topLeft A bit value of '1' indicates that the cell was computed
 * using the top-left cell.
 * @param left Character array that contains a '1' bit if the corresponding
 * cell was computed using its left cell.
 * @param top Character array that contains a '1' bit if the corresponding
 * cell was computed using its top cell.
 * @param qrySeq Query sequence.
 * @param refSeq Reference sequence.
 * @param match Match score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param maxScore Alignment score.
 * @param cellMax Index of the cell containing the max alignment score.
 */
__global__ void computeMaxScore_wrapFunc(int qryLen, int refLen, char *topLeft,
		char *left, char *top, const char *qrySeq, const char *refSeq,
		float match, float mismatch, float gapOpenPenalty, float gapExtPenalty,
		float *maxScore, int *cellMax)
{
	computeMaxScore(qryLen, refLen, topLeft, left, top, qrySeq, refSeq, match,
			mismatch, gapOpenPenalty, gapExtPenalty, maxScore, cellMax);
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
__device__ void computeMaxScore(int querySeqLength, int refSeqLength,
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
		querySeqShared[tx] = toupper(querySeq[querySeqStartPos + tx]);
	}

	/* Copy reference sequence from global memory to shared memory */
	int refSeqLength2 = refSeqLength + 1;
	if (tx < refSeqLength)
	{
		int refSeqStartPos = blockIdx.x * (refSeqLength + 1);
		refSeqShared[tx] = toupper(refSeq[refSeqStartPos + tx]);
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
//			if (refSeqShared[refIndex - 1] == querySeqShared[queryIndex - 1])
			if (isMatch(refSeqShared[refIndex - 1],
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
 * Computes maximum alignment score and creates a matrix that can be used
 * for backtracking to create alignment sequences.
 *
 * @param maxQryLen Query sequence length.
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
__device__ void computeMaxScore2(int maxQryLen, int refSeqLength,
		char *topLeft, char *left, char *top, const char *querySeq,
		const char *refSeq, float match, float mismatch, float gapOpenPenalty,
		float gapExtPenalty, float *maxScore, int *cellMax)
{
	__shared__ char querySeqShared[MAX_QRY_SEQ_LENGTH];
	__shared__ char refSeqShared[REF_LENGTH];
	__shared__ float HShared[3 * (MAX_QRY_SEQ_LENGTH + 1)];
	__shared__ float scores[MAX_QRY_SEQ_LENGTH + 1];
	__shared__ int cell[MAX_QRY_SEQ_LENGTH + 1];

	int tx = threadIdx.x;

	/* Copy query sequence from global memory to shared memory */
	if (tx < maxQryLen)
	{
		int querySeqStartPos = blockIdx.x * (maxQryLen + 1);
		querySeqShared[tx] = toupper(querySeq[querySeqStartPos + tx]);
	}

	/* Copy reference sequence from global memory to shared memory */
	int refSeqLength2 = refSeqLength + 1;
	if (tx < refSeqLength)
	{
		int refSeqStartPos = blockIdx.x * (refSeqLength + 1);
		refSeqShared[tx] = toupper(refSeq[refSeqStartPos + tx]);
	}

	int maxQryLen2 = maxQryLen + 1;
	if (tx < maxQryLen2)
	{
		HShared[(0 * maxQryLen2) + tx] = 0.0f;
		HShared[(1 * maxQryLen2) + tx] = 0.0f;
		HShared[(2 * maxQryLen2) + tx] = 0.0f;
		scores[tx] = 0.0f;
		cell[tx] = 0;
	}
	__syncthreads();


	int numIterations = refSeqLength + maxQryLen + 1;
	int currentHOffset, prevHOffset, prevPrevHOffset;
	float topLeftCellScore, topCellScore, leftCellScore, currentCellScore;
	int hLength = 3 * maxQryLen2;
	prevPrevHOffset = 0 * maxQryLen2;
	prevHOffset = 1 * maxQryLen2;
	currentHOffset = 2 * maxQryLen2;
	int refIndex, queryIndex, i, currentCellIndex, topCellIndex, leftCellIndex;
	int byte = 0, bit = 0;
	int charSize = (int) (sizeof(char) * NUM_BITS_IN_A_BYTE);

	for (i = 2; i < numIterations; ++i)
	{
		refIndex = i - tx;
		queryIndex = tx;
		if (refIndex > 0 && refIndex <= refSeqLength && queryIndex > 0
				&& queryIndex <= maxQryLen)
		{
			if (isMatch(refSeqShared[refIndex - 1],
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
		currentHOffset = currentHOffset + maxQryLen2;
		if (currentHOffset >= hLength)
			currentHOffset = 0;
	}


	if (tx == 0)
	{
		float scoreMax = 0.0f;
		*cellMax = 0;
		for (i = 0; i < maxQryLen2; ++i)
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
 * This is a wrapper function that wraps @a isMatch function. It has been
 * added so that @a isMatch can be unit-tested.
 *
 * @param		a		First base.
 * @param		b		Second base.
 * @param[out]	retVal	Return value.
 */
__global__ void isMatch_wrap(char *a, char *b, int *retVal)
{
	*retVal = isMatch(*a, *b);
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
__device__ int isMatch(char a, char b)
{
	if (a == b || a == 'N' || b == 'N')
		return 0;
	else
		return -1;
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
int isMatch_cpu(char a, char b)
{
	if (a == b || a == 'N' || b == 'N')
		return 0;
	else
		return -1;
}


/**
 * This is a wrapper function that wraps backtrack function.
 *
 * @param qryConsensus Query consensus sequence.
 * @param refConsensus Reference consensus sequence.
 * @param startPos Alignment start position.
 * @param endPos Alignment end position.
 * @param alignLen Alignment length.
 * @param topLeft A bit value of '1' indicates that the cell was computed
 * using the top-left cell.
 * @param left Character array that contains a '1' bit if the corresponding
 * cell was computed using its left cell.
 * @param top Character array that contains a '1' bit if the corresponding
 * cell was computed using its top cell.
 * @param qryLen Query sequence length.
 * @param refLen Reference sequence length.
 * @param cellMax Index of the cell containing the max alignment score.
 * @param qrySeq Query sequence.
 * @param refSeq Reference sequence.
 */
__global__ void backtrack_wrapFunc(char *qryConsensus, char *refConsensus,
		int *startPos, int *endPos, int *alignLen, char *topLeft, char *left,
		char *top, int qryLen, int refLen, int cellMax, char *qrySeq,
		char *refSeq)
{
	backtrack(qryConsensus, refConsensus, startPos, endPos, alignLen, topLeft,
			left, top, qryLen, refLen, cellMax, qrySeq, refSeq);
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
__device__ void backtrack(char *queryConsensus, char *refConsensus,
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
 * @param maxQryLen Length of query sequence.
 * @param refSeqLength Length of reference sequence.
 * @param cellMax Cell containing the max score.
 * @param querySeq Query sequence.
 * @param refSeq Reference sequence.
 */
__device__ void backtrack2(char *queryConsensus, char *refConsensus,
		int *startPos, int *endPos, int *alignLength, char *topLeft, char *left,
		char *top, int maxQryLen, int refSeqLength, int cellMax, char *querySeq,
		char *refSeq)
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
	querySeq += (blockIdx.x * (maxQryLen + 1));
	refConsensus += (blockIdx.x * (refSeqLength + maxQryLen + 2));
	queryConsensus += (blockIdx.x * (refSeqLength + maxQryLen + 2));

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
__device__ int toupper(int ch)
{
	if ((unsigned int)(ch - 'a') < 26u)
		ch += 'A' - 'a';
	return ch;
}


/**
 * Returns the upper-case of the given character.
 *
 * @param ch Character which needs to be converted to upper-case.
 * @return Upper-case letter.
 *
 * Code from koders.com.
 */
int toupper_cpu(int ch)
{
	if ((unsigned int)(ch - 'a') < 26u)
		ch += 'A' - 'a';
	return ch;
}




