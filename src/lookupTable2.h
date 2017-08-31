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

#ifndef LOOKUPTABLE2_H_
#define LOOKUPTABLE2_H_

#include "hitList.h"
#include "query.h"

void lookupTable2Create_wrap(const char *refFile, int seedLen,
		int maxHitsPerQry, int **lookupTable, uint **refPos,
		int *numDistinctTuples, int **numRepeatsPerTuple,
		int tupleIgnoreThreshold);

void lookupTable2Create(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold);

void lookupTable2Delete();

HitList *lookupTable2MapQry(const Query *query, uint qryLen, uint isRevComp);

static HitList *getBestScoringPosition(HitList *list);

static int compare(const void *a, const void *b);

#endif /* LOOKUPTABLE2_H_ */
