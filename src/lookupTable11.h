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

#ifndef LOOKUPTABLE11_H_
#define LOOKUPTABLE11_H_

int lookupTable11MapQry(char *qrySeq, int qryLen, char *refIdx_bestHits,
		int *refPos_bestHits);

int lookupTable11MapQry2(char *qrySeq, int qryLen, char *refIdx_bestHits,
		int *refPos_bestHits);

void lookupTable11Delete();

void lookupTable11Reset();

void lookupTable11Create(const char *keysFile, const char *valsFile,
		const char *keysFile2, int maxHits, int seedLen, int tupleIgnoreThres);

void lookupTable11Create(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold);

void lookupTable11Create2(const char *refFile, int seedLen, int maxHitsPerQry,
		int tupleIgnoreThreshold);

static int getClustSize(char *qrySeq, int qryLen, char refIdx, int refPos,
		int *isContinous);

static int getClustSize2(char *qrySeq, int qryLen, char refIdx, int refPos,
		int *isContinous);

#endif /* LOOKUPTABLE11_H_ */
