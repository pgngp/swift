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
#include <check.h>
#include "testArray.h"
#include "../src/array.h"


/**
 * Tests arrGetRandomNums function.
 */
START_TEST(getRandomNums)
{
	/* Test behavior when given array is NULL. */
	{
		uint *arr = NULL;
		arrGetRandomNums(1, 1, 10, NULL);
		if (arr != NULL)
			fail("Incorrect behavior when array is NULL.\n");
	}

	/* Test behavior when given number of random numbers is less than 1. */
	{
		int n = 0, i = 0;
		uint arr[10];
		for (i = 0; i < 10; ++i)
			arr[i] = 0;
		arrGetRandomNums(n, 1, 10, arr);
		for (i = 0; i < 10; ++i)
		{
			if (arr[i] != 0)
				fail("Incorrect behaviour when number of random numbers "
						"is less than 1.\n");
		}
	}

	/* Test behavior when given number of random numbers is equal to 1. */
	{
		int n = 1;
		uint arr[n];
		arrGetRandomNums(n, 1, 10, arr);
		if (arr[0] < 1 || arr[0] > 10)
			fail("Incorrect behavior when number of random numbers equals to"
					" 1.\n");
	}

	/* Test behavior when lower limit is strictly greater than upper limit. */
	{
		int n = 2, i = 0;
		uint arr[n];
		for (i = 0; i < n; ++i)
			arr[i] = 0;
		arrGetRandomNums(n, 10, 1, arr);
		for (i = 0; i < n; ++i)
		{
			if (arr[i] > 0)
				fail("Incorrect behavior when lower limit is strictly "
						"greater than upper limit.\n");
		}
	}

	/* Test behavior when number of random numbers requested is greater
	 * than the range. */
	{
		int n = 10;
		int lowerLimit = 1;
		int upperLimit = 5;
		uint arr[n];
		arrGetRandomNums(n, lowerLimit, upperLimit, arr);
	}

	/* Test behavior when lower limit is equal to upper limit. */
	{
		int n = 1, i = 0;
		uint arr[n];
		for (i = 0; i < n; ++i)
			arr[i] = 0;
		arrGetRandomNums(n, 1, 1, arr);
		for (i = 0; i < n; ++i)
		{
			if (arr[i] != 1)
				fail("Incorrect behavior when lower limit is equal to upper "
						"limit.\n");
		}
	}

	/* Test whether there are duplicates in the returned array. */
	{
		int n = 10;
		uint arr[n];
		arrGetRandomNums(n, 1, 10, arr);
		int i, j;
		for (i = 0; i < n; ++i)
		{
			for (j = 0; j < i; ++j)
			{
				if (arr[i] == arr[j])
					fail("Duplicate entries in the random number array.\n");

			}
		}
	}
}
END_TEST


/**
 * Tests arrSearch function.
 */
START_TEST(search)
{
	/* Test case when number is not in array. */
	{
		int num = 8;
		int arr[] = {4, 6, 9, 10};
		int arrSize = 4;
		int index = arrSearch(arr, arrSize, num);
		if (index != -1)
			fail("Incorrect behavior when number is not in array.\n");
	}

	/* Test case when number is in array. */
	{
		int num = 8;
		int arr[] = {4, 6, 8, 10};
		int arrSize = 4;
		int index = arrSearch(arr, arrSize, num);
		if (index != 2)
			fail("Incorrect behavior when number is present in array.\n");
	}

	/* Test case when array has size 1 and number is not present. */
	{
		int num = 8;
		int arr[] = {4};
		int arrSize = 1;
		int index = arrSearch(arr, arrSize, num);
		if (index != -1)
			fail("Incorrect behavior when array size is 1 and number "
					"is not present in array.\n");
	}

	/* Test case when array has size 1 and number is present. */
	{
		int num = 8;
		int arr[] = {8};
		int arrSize = 1;
		int index = arrSearch(arr, arrSize, num);
		if (index != 0)
			fail("Incorrect behavior when array size is 1 and number is "
					"present.\n");
	}

	/* Test case when array is NULL. */
	{
		int num = 8;
		int *arr = NULL;
		int arrSize = 2;
		int index = arrSearch(arr, arrSize, num);
		if (index != -1)
			fail("Incorrect behavior when array is NULL.\n");
	}
}
END_TEST


/**
 * Creates test suite.
 */
Suite *arraySuite(void)
{
	Suite *s = suite_create("array");

	/* Core test case */
	TCase *testCaseCore = tcase_create("Core");
	tcase_add_test(testCaseCore, getRandomNums);
	tcase_add_test(testCaseCore, search);
	suite_add_tcase(s, testCaseCore);

	return s;
}
