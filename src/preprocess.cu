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
#include "preprocess.h"
#include "common.h"


static int *_hashCodes = NULL; /* Array containing hash codes. */
static int *_powerVals = NULL; /* Array containing power values. */


/**
 * Creates and initializes data structures in this file.
 */
void preprocessCreate()
{
	_hashCodes = (int *) calloc(NUM_ASCII_CHARS, sizeof(int));
	_hashCodes['A'] = CODE_A;
	_hashCodes['a'] = CODE_A;
	_hashCodes['C'] = CODE_C;
	_hashCodes['c'] = CODE_C;
	_hashCodes['G'] = CODE_G;
	_hashCodes['g'] = CODE_G;
	_hashCodes['T'] = CODE_T;
	_hashCodes['t'] = CODE_T;
	_hashCodes['N'] = CODE_N;
	_hashCodes['n'] = CODE_N;

	_powerVals = (int *) calloc(MAX_SEED_LENGTH, sizeof(int));
	int i;
	for (i = 0; i < MAX_SEED_LENGTH; ++i)
		_powerVals[i] = (int) pow((float) 4, i);
}


/**
 * Deletes data structures in this file.
 */
void preprocessDelete()
{
	free(_hashCodes);
	_hashCodes = NULL;
	free(_powerVals);
	_powerVals = NULL;
}


/**
 * Returns the hash value of the given string
 *
 * @param s Given string
 * @param length Length of the given string
 * @return Hash value
 */
int getHash(char *s, int length)
{
	int i, sum;
	sum = 0;

	/* This section of code encodes the string into an integer.
	 * We use '4' as the base in the power function below because we
	 * have only 4 letters in our alphabet, i.e. A, C, G, and T.
	 *
	 * Example: Encoded value 'ACG' will be:
	 * = (4^0 * CODE_A) + (4^1 * CODE_C) + (4^2 * CODE_G)
	 * = (1 * CODE_A) + (4 * CODE_C) + (16 * CODE_G)
	 * (Assume CODE_A = 0, CODE_C = 1, and CODE_G = 2)
	 * = (1 * 0) + (4 * 1) + (16 * 2)
	 * = 0 + 4 + 32
	 * = 36
	 * That is, hash(ACG) = 36 */
	for (i = 0; i < length; ++i)
		sum += _powerVals[i] * _hashCodes[s[i]];

	return sum;
}


/**
 * Returns the hash value of the given string
 *
 * @param s Given string
 * @param length Length of the given string
 * @return Hash value
 */
long getLongHash(char *s, int length)
{
	int i;
	long sum = 0;

	/* This section of code encodes the string into an integer.
	 * We use '4' as the base in the power function below because we
	 * have only 4 letters in our alphabet, i.e. A, C, G, and T.
	 *
	 * Example: Encoded value 'ACG' will be:
	 * = (4^0 * CODE_A) + (4^1 * CODE_C) + (4^2 * CODE_G)
	 * = (1 * CODE_A) + (4 * CODE_C) + (16 * CODE_G)
	 * (Assume CODE_A = 0, CODE_C = 1, and CODE_G = 2)
	 * = (1 * 0) + (4 * 1) + (16 * 2)
	 * = 0 + 4 + 32
	 * = 36
	 * That is, hash(ACG) = 36 */
	for (i = 0; i < length; ++i)
		sum += _powerVals[i] * _hashCodes[s[i]];

	return sum;
}



///**
// * Returns the hash value of the given string
// *
// * @param s Given string
// * @param length Length of the given string
// * @return Hash value
// */
//int getHash(char *s, int length)
//{
//	static int i, sum;
//	sum = 0;
//
//	/* This section of code encodes the string into an integer.
//	 * We use '4' as the base in the power function below because we
//	 * have only 4 letters in our alphabet, i.e. A, C, G, and T.
//	 *
//	 * Example: Encoded value 'ACG' will be:
//	 * = (4^0 * CODE_A) + (4^1 * CODE_C) + (4^2 * CODE_G)
//	 * = (1 * CODE_A) + (4 * CODE_C) + (16 * CODE_G)
//	 * (Assume CODE_A = 0, CODE_C = 1, and CODE_G = 2)
//	 * = (1 * 0) + (4 * 1) + (16 * 2)
//	 * = 0 + 4 + 32
//	 * = 36
//	 * That is, hash(ACG) = 36 */
//	for (i = 0; i < length; ++i)
//	{
//		if (s[i] == 'A' || s[i] == 'a')
//			sum += (int) pow((float) 4, i) * CODE_A;
//		else if (s[i] == 'C' || s[i] == 'c')
//			sum += (int) pow((float) 4, i) * CODE_C;
//		else if (s[i] == 'G' || s[i] == 'g')
//			sum += (int) pow((float) 4, i) * CODE_G;
//		else if (s[i] == 'T' || s[i] == 't')
//			sum += (int) pow((float) 4, i) * CODE_T;
//		else if (s[i] == 'N' || s[i] == 'n')
//			sum += (int) pow((float) 4, i) * CODE_N;
//		else if (s[i] == '*')
//			sum += (int) pow((float) 4, i) * CODE_STAR;
//		else
//			sum += 0;
//	}
//
//	return sum;
//}


/**
 * Decodes the given hash
 *
 * @param hash Hash value to be decoded
 * @param[out] s Decoded sequence
 * @param seedLength Seed length
 */
void getUnhash(int hash, char *s, int seedLength)
{
	int quotient = hash, remainder, i = -1;

	while (quotient > 0)
	{
		remainder = quotient % 4;
		++i;
		if (remainder == 0)
			s[i] = 'A';
		else if (remainder == 1)
			s[i] = 'C';
		else if (remainder == 2)
			s[i] = 'G';
		else if (remainder == 3)
			s[i] = 'T';
		else
			s[i] = 'N';
		quotient = quotient / 4;
	}

	/* Pad with 'A's */
	while (i < (seedLength - 1))
	{
		++i;
		s[i] = 'A';
	}
	s[seedLength] = '\0';
}
