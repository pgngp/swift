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

#ifndef QUERY_H_
#define QUERY_H_

#include <stdio.h>


/**
 * @struct queryStruct
 *
 * Stores query data.
 */
typedef struct queryStruct
{
	char *name; /**< Name of the query. */
	int length; /**< Length of the query sequence. */
	char *seq; /**< Sequence of the query. */
	char *revCompSeq; /**< Reverse complement sequence of the query. */
	ulong nameOffset; /**< File offset of the query name. */
	ulong seqOffset; /**< File offset of the query sequence. */
} Query;


void qryInitialize(char *f, int length);
void qryClean();
Query *qryGetNext();
Query *qryGetNext2(FILE *qryFilePtr, Query *qry);
Query *qryGetNext3();
void qryGetReverseComplement(char *seq, int seqLength, char *revComp);
Query *qryCreate(const char *name, const char *seq, const char *revCompSeq,
		ulong nameOffset, ulong seqOffset);
Query *qryCreate2();
void qryDelete(Query *qry);

#endif /* QUERY_H_ */
