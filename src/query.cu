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
#include "query.h"
#include "common.h"
#include <pthread.h>


char qryFileName[MAX_FILE_NAME_LENGTH];
FILE *qryFilePtr;
int qryLength = 0;
Query *qry;
static char line[MAX_LINE_LENGTH];
static char *line2;
static int lineLength;
pthread_mutex_t mutex;


/**
 * Initializes a query
 *
 * @param f Query file name
 * @param length Length of each query
 */
void qryInitialize(char *f, int length)
{
	strcpy(qryFileName, f);
	qryFilePtr = fopen(qryFileName, "r");
	qryLength = length;
	qry = (Query *) malloc(sizeof(Query));
	qry->name = (char *) malloc(MAX_QRY_NAME_LENGTH * sizeof(char));
	qry->seq = (char *) malloc((qryLength + 1) * sizeof(char));
	qry->revCompSeq = (char *) malloc((qryLength + 1) * sizeof(char));
}


/**
 * Cleans up: closes the query file and releases objects
 */
void qryClean()
{
	fclose(qryFilePtr);
	free(qry->name);
	free(qry->seq);
	free(qry->revCompSeq);
	free(qry);
}


/**
 * Creates a @a Query object using the given parameters and returns a pointer
 * to it.
 *
 * @param name Name of the query. A copy of this string will be created.
 * @param seq Query sequence. A copy of this string will be created.
 * @param revCompSeq Reverse complement of the query sequence. A copy of
 * this string will be created.
 * @param nameOffset Offset where the query name appears in the file.
 * @param seqOffset Offset where the query sequence appears in the file.
 * @return Query object.
 */
Query *qryCreate(const char *name, const char *seq, const char *revCompSeq,
		ulong nameOffset, ulong seqOffset)
{
	if (name == NULL || seq == NULL || revCompSeq == NULL
			|| nameOffset >= seqOffset)
		return NULL;

	Query *qry = (Query *) malloc(sizeof(Query));
	qry->name = (char *) malloc(sizeof(char) * strlen(name));
	qry->seq = (char *) malloc(sizeof(char) * (strlen(seq) + 1));
	qry->revCompSeq = (char *) malloc(sizeof(char) * (strlen(revCompSeq) + 1));
	qry->name = strcpy(qry->name, name);
	qry->seq = strcpy(qry->seq, seq);
	qry->revCompSeq = strcpy(qry->revCompSeq, revCompSeq);
	qry->nameOffset = nameOffset;
	qry->seqOffset = seqOffset;
	return qry;
}


/**
 * Creates a @a Query object with default values.
 *
 * @return A new @a Query object.
 */
Query *qryCreate2()
{
	Query *qry = (Query *) malloc(sizeof(Query));
	qry->name = (char *) malloc(sizeof(char) * MAX_QRY_NAME_LENGTH);
	qry->length = 0;
	qry->seq = (char *) malloc(sizeof(char) * MAX_QRY_SEQ_LENGTH);
	qry->revCompSeq = (char *) malloc(sizeof(char) * MAX_QRY_SEQ_LENGTH);
	qry->nameOffset = 0;
	qry->seqOffset = 0;
	return qry;
}


/**
 * Frees up resources occupied by the given @a Query object.
 */
void qryDelete(Query *qry)
{
	free(qry->name);
	free(qry->seq);
	free(qry->revCompSeq);
	free(qry);
}


/**
 * Returns a pointer to the next query sequence in the query file
 *
 * @return Pointer to the next query
 */
Query *qryGetNext()
{
	while (fgets(line, MAX_LINE_LENGTH, qryFilePtr) != NULL)
	{
		lineLength = strlen(line);
		if (line[lineLength - 1] == '\n')
		{
			line[lineLength - 1] = '\0';
			--lineLength;
		}

		/* Empty line. */
		if (lineLength == 0)
			continue;
		/* Line with query ID. */
		else if (line[0] == '>')
		{
			strcpy(qry->name, line + 1);
			qry->nameOffset = ftell(qryFilePtr) - lineLength - 1;
			qry->seqOffset = ftell(qryFilePtr);
		}
		/* Line with query sequence. */
		else
		{
			strncpy(qry->seq, line, lineLength);
			qry->seq[lineLength] = '\0';
			qryGetReverseComplement(qry->seq, qryLength, qry->revCompSeq);
			break;
		}
	}

	return qry;
}


/**
 * Returns a pointer to the next query sequence in the query file.
 *
 * @param qryFilePtr Pointer to the query sequence file.
 * @param qry @a Query object in which query data members will be stored.
 * @return Pointer to the next query.
 */
Query *qryGetNext2(FILE *qryFilePtr, Query *qry)
{
	line2 = fgets(line, MAX_LINE_LENGTH, qryFilePtr);

	while (line2 != NULL)
	{
		lineLength = strlen(line);
		if (line[lineLength - 1] == '\n')
		{
			line[lineLength - 1] = '\0';
			--lineLength;
		}

		/* Empty line. */
		if (lineLength == 0)
			continue;
		/* Line with query ID. */
		else if (line[0] == '>')
		{
			strcpy(qry->name, line + 1);
			qry->nameOffset = ftell(qryFilePtr) - lineLength - 1;
			qry->seqOffset = ftell(qryFilePtr);
		}
		/* Line with query sequence. */
		else
		{
			strncpy(qry->seq, line, lineLength);
			qry->seq[lineLength] = '\0';
			qry->length = strlen(qry->seq);
			qryGetReverseComplement(qry->seq, qry->length, qry->revCompSeq);
			break;
		}
		line2 = fgets(line, MAX_LINE_LENGTH, qryFilePtr);
	}

	if (line2 == NULL)
		return NULL;

	return qry;
}


/**
 * Returns a pointer to the next query sequence in the query file
 *
 * @return Pointer to the next query
 */
Query *qryGetNext3()
{
	pthread_mutex_lock(&mutex);
	while (fgets(line, MAX_LINE_LENGTH, qryFilePtr) != NULL)
	{
		lineLength = strlen(line);
		if (line[lineLength - 1] == '\n')
		{
			line[lineLength - 1] = '\0';
			--lineLength;
		}

		/* Empty line. */
		if (lineLength == 0)
			continue;
		/* Line with query ID. */
		else if (line[0] == '>')
		{
			strcpy(qry->name, line + 1);
			qry->nameOffset = ftell(qryFilePtr) - lineLength - 1;
			qry->seqOffset = ftell(qryFilePtr);
		}
		/* Line with query sequence. */
		else
		{
			strncpy(qry->seq, line, lineLength);
			qry->seq[lineLength] = '\0';
			qryGetReverseComplement(qry->seq, qryLength, qry->revCompSeq);
			break;
		}
	}
	pthread_mutex_unlock(&mutex);

	return qry;
}


/**
 * Returns the reverse complement of the given query sequence
 *
 * @param seq Query sequence
 * @param seqLength Length of sequence
 * @param[out] revComp Reverse complement of the given sequence
 */
void qryGetReverseComplement(char *seq, int seqLength, char *revComp)
{
	int i, index = seqLength;
	char c;
	for (i = 0; i < seqLength; ++i)
	{
		c = seq[--index];
		if (c == 'A')
			revComp[i] = 'T';
		else if (c == 'C')
			revComp[i] = 'G';
		else if (c == 'G')
			revComp[i] = 'C';
		else
			revComp[i] = 'A';
	}
	revComp[seqLength] = '\0';
}
