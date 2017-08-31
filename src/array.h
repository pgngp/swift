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
 * Defines functions related to an array.
 */

#ifndef ARRAY_H_
#define ARRAY_H_

void arrGetRandomNums2(int n, int lowerLimit, int upperLimit, int *arr);
void arrGetRandomNums(uint n, uint lowerLimit, uint upperLimit, uint *arr);
int arrSearch(int *arr, int arrSize, int num);
//__device__ void arrGetRandomNums_gpu(uint n, uint lowerLimit, uint upperLimit,
//		uint *arr);
#endif /* ARRAY_H_ */
