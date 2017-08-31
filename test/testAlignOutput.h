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
 * This file contains code to validate the alignment output file.
 */

#ifndef TESTALIGNOUTPUT_H_
#define TESTALIGNOUTPUT_H_

static void validateOutput(const char *refFile, const char *qryFile,
		const char *alignFile, const char *outputFile);

static void validateMatchInfo();

static void validateRefSeq();

static void validateQrySeq();

const char *getUsage();

static void getRefSeq(char *refSeq);

#endif /* TESTALIGNOUTPUT_H_ */
