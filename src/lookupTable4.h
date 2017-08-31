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

#ifndef LOOKUPTABLE4_H_
#define LOOKUPTABLE4_H_

#include "query.h"

void lookupTable4Create_wrap(const char *refFile, int seedLen,
		int maxHitsPerQry, int **lookupTable, char **refIdx, int **refPos,
		int *numDistinctTuples, int **numRepeatsPerTuple,
		int tupleIgnoreThreshold);

void lookupTable4Create(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold);

void lookupTable4Delete();

int lookupTable4MapQry(const Query *query, uint qryLen, uint isRevComp,
		char *refIdx, int *shift, int *refPos);

#endif /* LOOKUPTABLE4_H_ */
