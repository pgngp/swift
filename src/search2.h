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

#ifndef SEARCH2_H_
#define SEARCH2_H_

#include "query.h"

static int getNumQrsPerIter(int seedLen, int numTotalTuples, int maxNumHits,
		int maxRefPosComposite);

void search2Queries(char *qryFile, char *refFile, char *matchFile, int seedLen,
		int maxNumHits, int *numMatches, int tupleIgnoreThres);

static void printHeuristicMatches(Query **qry, char *refIdx, int *refPos,
		int numQrsPerIter, int maxNumHits, int *numMatches, FILE *filePtr);

#endif /* SEARCH2_H_ */
