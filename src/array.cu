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

#include "array.h"


/**
 * Fetches the given number of random numbers between the given limits.
 *
 * @param n Number of random numbers to be created. This number should not
 * be greater than the range specified by lower and upper limit.
 * @param lowerLimit The lower limit of the random numbers.
 * @param upperLimit The upper limit of the random numbers.
 * @param[out] arr Array in which the random numbers will be stored. The size
 * of the array should be atleast @a n.
 */
void arrGetRandomNums2(int n, int lowerLimit, int upperLimit, int *arr)
{
	int range = upperLimit - lowerLimit + 1;
	if (n <= 0 || lowerLimit > upperLimit || arr == NULL || n > range)
		return;

	int i = 0, randNum;
	while (i < n)
	{
		randNum = (rand() % range) + lowerLimit;
		if (arrSearch((int *) arr, i, randNum) == -1)
		{
			arr[i] = randNum;
			++i;
		}
	}
}



/**
 * Fetches the given number of random numbers between the given limits.
 *
 * @param n Number of random numbers to be created. This number should not
 * be greater than the range specified by lower and upper limit.
 * @param lowerLimit The lower limit of the random numbers.
 * @param upperLimit The upper limit of the random numbers.
 * @param[out] arr Array in which the random numbers will be stored. The size
 * of the array should be atleast @a n.
 */
void arrGetRandomNums(uint n, uint lowerLimit, uint upperLimit, uint *arr)
{
	int range = upperLimit - lowerLimit + 1;
	if (n <= 0 || lowerLimit > upperLimit || arr == NULL || n > range)
		return;

	int i = 0, randNum;
	while (i < n)
	{
		randNum = (rand() % range) + lowerLimit;
		if (arrSearch((int *) arr, i, randNum) == -1)
		{
			arr[i] = randNum;
			++i;
		}
	}
}


///**
// * Fetches the given number of random numbers between the given limits.
// *
// * @param n Number of random numbers to be created. This number should not
// * be greater than the range specified by lower and upper limit.
// * @param lowerLimit The lower limit of the random numbers.
// * @param upperLimit The upper limit of the random numbers.
// * @param[out] arr Array in which the random numbers will be stored. The size
// * of the array should be atleast @a n.
// */
//__device__ void arrGetRandomNums_gpu(uint n, uint lowerLimit, uint upperLimit,
//		uint *arr)
//{
//	int range = upperLimit - lowerLimit + 1;
//	if (n <= 0 || lowerLimit > upperLimit || arr == NULL || n > range)
//		return;
//
//	int i = 0, randNum;
//	while (i < n)
//	{
//		randNum = (rand() % range) + lowerLimit;
//		if (arrSearch((int *) arr, i, randNum) == -1)
//		{
//			arr[i] = randNum;
//			++i;
//		}
//	}
//}


/**
 * Returns the index of the given number if it is already present in the
 * given array; otherwise, returns -1.
 *
 * @param arr Array in which the number will be searched.
 * @param arrSize Number of elements in the array.
 * @param num The number to be searched in the array.
 * @return Index of the searched number in the array; otherwise, -1.
 */
int arrSearch(int *arr, int arrSize, int num)
{
	if (arr == NULL || arrSize == 0)
		return -1;

	int i = 0;
	for (i = 0; i < arrSize; ++i)
	{
		if (arr[i] == num)
			return i;
	}

	return -1;
}

