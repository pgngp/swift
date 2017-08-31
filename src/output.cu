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

#include <stdio.h>
#include <ctype.h>
#include "output.h"
#include "common.h"


/**
 * Prints the cigar string using the given reference and query consensus
 * sequences
 *
 * @param queryConsensus Character pointer to query consensus sequence
 * @param refConsensus Character pointer to reference consensus sequence
 * @param alignLength Length of consensus sequences
 * @param outputFilePtr Output file stream.
 */
void printCigar(
		char *queryConsensus,
		char *refConsensus,
		int alignLength,
		FILE *outputFilePtr)
{
	int operationLength = 0, z;
	char operation = ' ';
	for(z = alignLength - 1; z >= 0; --z)
	{
		/* Deletion from the reference */
		if (refConsensus[z] == '-')
		{
			if (operationLength > 0 && operation != 'D')
			{
				fprintf(outputFilePtr, "%d%c", operationLength, operation);
				operationLength = 0;
			}
			operation = 'D';
			++operationLength;
		}
		/* Insertion to the reference */
		else if (queryConsensus[z] == '-')
		{
			if (operationLength > 0 && operation != 'I')
			{
				fprintf(outputFilePtr, "%d%c", operationLength, operation);
				operationLength = 0;
			}
			operation = 'I';
			++operationLength;
		}
		/* Alignment match/mismatch */
		else
		{
			if (operationLength > 0 && operation != 'M')
			{
				fprintf(outputFilePtr, "%d%c", operationLength, operation);
				operationLength = 0;
			}
			operation = 'M';
			++operationLength;
		}
	}
	if (operationLength > 0)
		fprintf(outputFilePtr, "%d%c", operationLength, operation);
}


/**
 * Prints output in SAM format.
 *
 * @param numQueries Number of query sequences
 * @param queryName Character pointer to query name
 * @param refName Character pointer to reference name
 * @param score Floating point pointer to score array
 * @param alignStart Integer pointer to alignment start position array
 * @param alignEnd Integer pointer to alignment end position array
 * @param alignLength Integer pointer to alignment length array
 * @param combinedLength Combined length of reference and query sequence
 * @param queryConsensus Character pointer to query consensus sequence
 * @param refConsensus Character pointer to reference consensus sequence
 * @param outputFilePtr File descriptor of the file where alignments
 * will be saved.
 */
void printOutputSam(
		int numQueries,
		char *queryName,
		char *refName,
		float *score,
		int *alignStart,
		int *alignEnd,
		int *alignLength,
		int combinedLength,
		char *queryConsensus,
		char *refConsensus,
		FILE *outputFilePtr)
{
	int y, z;
	for (y = 0; y < numQueries; ++y)
	{
		fprintf(outputFilePtr,
				"%s\t*\t%s\t%d\t255\t",
				queryName + (y * MAX_SEQ_NAME_LENGTH),
				refName + (y * MAX_SEQ_NAME_LENGTH),
				alignStart[y]);
		printCigar(
				queryConsensus + (y * combinedLength),
				refConsensus + (y * combinedLength),
				(int) alignLength[y],
				outputFilePtr);
		fprintf(outputFilePtr, "\t*\t0\t0\t");
		for(z = (y * combinedLength) + alignLength[y] - 1;
				z >= (y * combinedLength); --z)
			fprintf(outputFilePtr, "%c", queryConsensus[z]);
		fprintf(outputFilePtr, "\t*\n");
	}
}


/**
 * Prints alignment output.
 *
 * Output includes query name, reference name, alignment score,
 * reference start position, reference end position, alignment length,
 * and the aligned reference and query sequences.
 *
 * @param numQueries Number of queries for which the output will be printed
 * @param queryName Array containing query names
 * @param refName Array containing reference names
 * @param score Array containing alignment scores
 * @param refDistance Array containing reference distances
 * @param alignStart Array containing alignment start position on the reference
 * @param alignEnd Array containing alignment end position on the reference
 * @param alignLength Array containing the alignment lengths
 * @param combinedLength Combined length of reference and query sequence
 * @param queryConsensus Array containing the aligned query sequence
 * @param refConsensus Array containing the aligned reference sequence
 * @param outputFilePtr File descriptor of the file where alignments
 * will be saved.
 */
void printOutput(int numQueries, char *queryName, char *refName, float *score,
		int *refDistance, int *alignStart, int *alignEnd, int *alignLength,
		int combinedLength, char *queryConsensus, char *refConsensus,
		FILE *outputFilePtr)
{
	int y, z;
	for (y = 0; y < numQueries; ++y)
	{
		if (score[y] < 0.0)
			continue;

		fprintf(outputFilePtr, "\n****************************************\n");
		fprintf(outputFilePtr,
				"Qry: %s | Ref: %s | Score: %.1f | AlignStart: %d | "
				"AlignEnd: %d | AlignLength: %d\n",
				queryName + (y * MAX_SEQ_NAME_LENGTH),
				refName + (y * MAX_SEQ_NAME_LENGTH),
				score[y],
				alignStart[y] + refDistance[y],
				alignEnd[y] + refDistance[y],
				alignLength[y]);
		fprintf(outputFilePtr, "R: ");
		for (z = (y * combinedLength) + alignLength[y] - 1;
				z >= (y * combinedLength); --z)
			fprintf(outputFilePtr, "%c", refConsensus[z]);
		fprintf(outputFilePtr, "\n");

		fprintf(outputFilePtr, "   ");
		for (z = (y * combinedLength) + alignLength[y] - 1;
				z >= (y * combinedLength); --z)
		{
			if (toupper(refConsensus[z]) == toupper(queryConsensus[z])
					|| refConsensus[z] == 'N' || queryConsensus[z] == 'N')
				fprintf(outputFilePtr, "|");
			else if (refConsensus[z] == '-' || queryConsensus[z] == '-')
				fprintf(outputFilePtr, "-");
			else
				fprintf(outputFilePtr, ".");
		}
		fprintf(outputFilePtr, "\n");

		fprintf(outputFilePtr, "Q: ");
		for (z = (y * combinedLength) + alignLength[y] - 1;
				z >= (y * combinedLength); --z)
			fprintf(outputFilePtr, "%c", queryConsensus[z]);
		fprintf(outputFilePtr, "\n");
	}
}


/**
 * Prints alignment output.
 *
 * Output includes query name, reference name, alignment score,
 * reference start position, reference end position, alignment length,
 * and the aligned reference and query sequences.
 *
 * @param numQueries Number of queries for which the output will be printed
 * @param queryName Array containing query names
 * @param refName Array containing reference names
 * @param score Array containing alignment scores
 * @param refDistance Array containing reference distances
 * @param alignStart Array containing alignment start position on the reference
 * @param alignEnd Array containing alignment end position on the reference
 * @param alignLength Array containing the alignment lengths
 * @param combinedLength Combined length of reference and query sequence
 * @param queryConsensus Array containing the aligned query sequence
 * @param refConsensus Array containing the aligned reference sequence
 * @param outputFilePtr File descriptor of the file where alignments
 * will be saved.
 */
void printOutput_paired(int numQueries, char *queryName, char *refName,
		float *score, int *refDistance, int *alignStart, int *alignEnd,
		int *alignLength, int combinedLength, char *queryConsensus,
		char *refConsensus, FILE *outputFilePtr)
{
	int y, z, count = 0;
	for (y = 0; y < numQueries; ++y)
	{
		if (score[y] == -1.0)
			continue;

		++count;
		if (count % 2 == 1)
			fprintf(outputFilePtr, "\n*************************************\n");
		else
			fprintf(outputFilePtr, "\n");
		fprintf(outputFilePtr,
				"Qry: %s | Ref: %s | Score: %.1f | AlignStart: %d |"
				" AlignEnd: %d | AlignLength: %d\n",
				queryName + (y * MAX_SEQ_NAME_LENGTH),
				refName + (y * MAX_SEQ_NAME_LENGTH),
				score[y],
				alignStart[y] + refDistance[y],
				alignEnd[y] + refDistance[y],
				alignLength[y]);
		fprintf(outputFilePtr, "R: ");
		for(z = (y * combinedLength) + alignLength[y] - 1;
				z >= (y * combinedLength); --z)
			fprintf(outputFilePtr, "%c", refConsensus[z]);
		fprintf(outputFilePtr, "\n");

		fprintf(outputFilePtr, "   ");
		for (z = (y * combinedLength) + alignLength[y] - 1;
				z >= (y * combinedLength); --z)
		{
			if (toupper(refConsensus[z]) == toupper(queryConsensus[z])
					|| refConsensus[z] == 'N' || queryConsensus[z] == 'N')
				fprintf(outputFilePtr, "|");
			else if (refConsensus[z] == '-' || queryConsensus[z] == '-')
				fprintf(outputFilePtr, "-");
			else
				fprintf(outputFilePtr, ".");
		}
		fprintf(outputFilePtr, "\n");

		fprintf(outputFilePtr, "Q: ");
		for(z = (y * combinedLength) + alignLength[y] - 1;
				z >= (y * combinedLength); --z)
			fprintf(outputFilePtr, "%c", queryConsensus[z]);
		fprintf(outputFilePtr, "\n");
	}
}


/**
 * Prints alignment output.
 *
 * Output includes query name, reference name, alignment score,
 * reference start position, reference end position, alignment length,
 * and the aligned reference and query sequences.
 *
 * @param numQueries Number of queries for which the output will be printed
 * @param queryName Array containing query names
 * @param refName Array containing reference names
 * @param score Array containing alignment scores
 * @param refDistance Array containing reference distances
 * @param alignStart Array containing alignment start position on the reference
 * @param alignEnd Array containing alignment end position on the reference
 * @param alignLength Array containing the alignment lengths
 * @param combinedLength Combined length of reference and query sequence
 * @param queryConsensus Array containing the aligned query sequence
 * @param refConsensus Array containing the aligned reference sequence
 * @param outputFilePtr File descriptor of the file where alignments
 * will be saved.
 */
void printOutput_paired2(int numQueries, int *qryIdx, char *queryName,
        char *refName, float *score, int *refDistance, int *alignStart,
        int *alignEnd, int *alignLength, int combinedLength,
        char *queryConsensus, char *refConsensus, FILE *outputFilePtr)
{
    static int y, z, oldQryIdx;
    oldQryIdx = -1;
    for (y = 0; y < numQueries; ++y)
    {
        if (score[y] == -1.0)
            continue;

        if (qryIdx[y] != oldQryIdx)
        {
            oldQryIdx = qryIdx[y];
            fprintf(outputFilePtr, "\n*************************************\n");
        }
        else
            fprintf(outputFilePtr, "\n");
        fprintf(outputFilePtr,
                "Qry: %s | Ref: %s | Score: %.1f | AlignStart: %d |"
                " AlignEnd: %d | AlignLength: %d\n",
                queryName + (y * MAX_SEQ_NAME_LENGTH),
                refName + (y * MAX_SEQ_NAME_LENGTH),
                score[y],
                alignStart[y] + refDistance[y],
                alignEnd[y] + refDistance[y],
                alignLength[y]);
        fprintf(outputFilePtr, "R: ");
        for(z = (y * combinedLength) + alignLength[y] - 1;
                z >= (y * combinedLength); --z)
            fprintf(outputFilePtr, "%c", refConsensus[z]);
        fprintf(outputFilePtr, "\n");

        fprintf(outputFilePtr, "   ");
        for (z = (y * combinedLength) + alignLength[y] - 1;
                z >= (y * combinedLength); --z)
        {
            if (toupper(refConsensus[z]) == toupper(queryConsensus[z])
                    || refConsensus[z] == 'N' || queryConsensus[z] == 'N')
                fprintf(outputFilePtr, "|");
            else if (refConsensus[z] == '-' || queryConsensus[z] == '-')
                fprintf(outputFilePtr, "-");
            else
                fprintf(outputFilePtr, ".");
        }
        fprintf(outputFilePtr, "\n");

        fprintf(outputFilePtr, "Q: ");
        for(z = (y * combinedLength) + alignLength[y] - 1;
                z >= (y * combinedLength); --z)
            fprintf(outputFilePtr, "%c", queryConsensus[z]);
        fprintf(outputFilePtr, "\n");
    }
}


/**
 * Prints alignment output.
 *
 * Output includes query name, reference name, alignment score,
 * reference start position, reference end position, alignment length,
 * and the aligned reference and query sequences.
 *
 * @param numQueries Number of queries for which the output will be printed
 * @param queryName Array containing query names
 * @param refName Array containing reference names
 * @param score Array containing alignment scores
 * @param refDistance Array containing reference distances
 * @param alignStart Array containing alignment start position on the reference
 * @param alignEnd Array containing alignment end position on the reference
 * @param alignLength Array containing the alignment lengths
 * @param combinedLength Combined length of reference and query sequence
 * @param queryConsensus Array containing the aligned query sequence
 * @param refConsensus Array containing the aligned reference sequence
 * @param outputFilePtr File descriptor of the file where alignments
 * will be saved.
 */
void printOutput_paired_cpu(int numQueries, char *queryName, char *refName,
		float *score, int *refDistance, int *alignStart, int *alignEnd,
		int *alignLength, int combinedLength, char *queryConsensus,
		char *refConsensus, FILE *outputFilePtr, int isFirstQry)
{
	int y, z;
	for (y = 0; y < numQueries; ++y)
	{
		if (isFirstQry == 1)
			fprintf(outputFilePtr, "\n*************************************\n");
		else
			fprintf(outputFilePtr, "\n");
		fprintf(outputFilePtr,
				"Ref: %s, Query: %s, Score: %f, Align start: %d, "
				"Align end: %d, Align length: %d\n",
				queryName + (y * MAX_SEQ_NAME_LENGTH),
				refName + (y * MAX_SEQ_NAME_LENGTH),
				score[y],
				alignStart[y] + refDistance[y] - NUM_NEW_LINES_BUFFER,
				alignEnd[y] + refDistance[y] - NUM_NEW_LINES_BUFFER,
				alignLength[y]);
		fprintf(outputFilePtr, "R: ");
		for(z = (y * combinedLength) + alignLength[y] - 1;
				z >= (y * combinedLength); --z)
			fprintf(outputFilePtr, "%c", refConsensus[z]);
		fprintf(outputFilePtr, "\n");
		fprintf(outputFilePtr, "Q: ");
		for(z = (y * combinedLength) + alignLength[y] - 1;
				z >= (y * combinedLength); --z)
			fprintf(outputFilePtr, "%c", queryConsensus[z]);
		fprintf(outputFilePtr, "\n");
	}
}


