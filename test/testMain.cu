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

#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include "testHitList.h"
#include "testRefPosList.h"
#include "testReference.h"
#include "testArray.h"
#include "testInput.h"
#include "testQuery.h"
#include "testSearch.h"
#include "testAlign.h"
#include "testSmithWaterman.h"
#include "testRefPosMap.h"
#include "testRefNameMap.h"
#include "testRefMap.h"
#include "testLookupTable.h"
#include "testLookupTable5.h"
#include "testLookupTable6.h"
#include "testLookupTable7.h"
#include "testMapHits.h"
#include "testMapHits2.h"
#include "testMapHits3.h"
#include "testMapHits4.h"
#include "testMapHits5.h"
#include "testMapHits6.h"


/**
 * Main function.
 *
 * Runs test cases.
 */
int main(void)
{
	srand(time(NULL));
	int numberFailed = 0;

	/* HitList test suite */
	Suite *hitListTestSuite = hitListSuite();
	SRunner *sr = srunner_create(hitListTestSuite);

	/* RefPosList test suite */
	Suite *refPosListTestSuite = refPosListSuite();
	srunner_add_suite(sr, refPosListTestSuite);

//	/* Reference test suite */
//	Suite *referenceTestSuite = referenceSuite();
//	srunner_add_suite(sr, referenceTestSuite);

	/* Array test suite */
	Suite *arrayTestSuite = arraySuite();
	srunner_add_suite(sr, arrayTestSuite);

	/* Input test suite */
	Suite *inputTestSuite = inputSuite();
	srunner_add_suite(sr, inputTestSuite);

	/* Query test suite */
	Suite *queryTestSuite = querySuite();
	srunner_add_suite(sr, queryTestSuite);

	/* Search test suite */
	Suite *searchTestSuite = searchSuite();
	srunner_add_suite(sr, searchTestSuite);

	/* Align test suite */
	Suite *alignTestSuite = alignSuite();
	srunner_add_suite(sr, alignTestSuite);

	/* SmithWaterman test suite */
	Suite *smithWatermanTestSuite = smithWatermanSuite();
	srunner_add_suite(sr, smithWatermanTestSuite);

	/* RefPosMap test suite */
//	Suite *refPosMapTestSuite = refPosMapSuite();
//	srunner_add_suite(sr, refPosMapTestSuite);

	/* RefNameMap test suite */
	Suite *refNameMapTestSuite = refNameMapSuite();
	srunner_add_suite(sr, refNameMapTestSuite);

	/* RefMap test suite */
	Suite *refMapTestSuite = refMapSuite();
	srunner_add_suite(sr, refMapTestSuite);

	/* Lookup table test suite */
	Suite *lookupTableTestSuite = lookupTableSuite();
	srunner_add_suite(sr, lookupTableTestSuite);

	/* LookupTable5 test suite */
	Suite *lookupTable5TestSuite = lookupTable5Suite();
	srunner_add_suite(sr, lookupTable5TestSuite);

	/* LookupTable6 test suite */
	Suite *lookupTable6TestSuite = lookupTable6Suite();
	srunner_add_suite(sr, lookupTable6TestSuite);

	/* LookupTable7 test suite */
	Suite *lookupTable7TestSuite = lookupTable7Suite();
	srunner_add_suite(sr, lookupTable7TestSuite);

	/* MapHits test suite */
	Suite *mapHitsTestSuite = mapHitsSuite();
	srunner_add_suite(sr, mapHitsTestSuite);

//	/* MapHits2 test suite */
//	Suite *mapHits2TestSuite = mapHits2Suite();
//	srunner_add_suite(sr, mapHits2TestSuite);

//	/* MapHits3 test suite */
//	Suite *mapHits3TestSuite = mapHits3Suite();
//	srunner_add_suite(sr, mapHits3TestSuite);

//	/* MapHits4 test suite */
//	Suite *mapHits4TestSuite = mapHits4Suite();
//	srunner_add_suite(sr, mapHits4TestSuite);

	/* MapHits5 test suite */
	Suite *mapHits5TestSuite = mapHits5Suite();
	srunner_add_suite(sr, mapHits5TestSuite);

	/* MapHits6 test suite */
	Suite *mapHits6TestSuite = mapHits6Suite();
	srunner_add_suite(sr, mapHits6TestSuite);

	srunner_run_all(sr, CK_NORMAL);
	numberFailed = srunner_ntests_failed(sr);
	srunner_free(sr);

	if (numberFailed == 0)
	{
		fprintf(stderr, "******************\n");
		fprintf(stderr, "PASSED!\n");
		fprintf(stderr, "******************\n");
	}
	else
	{
		fprintf(stderr, "******************\n");
		fprintf(stderr, "FAILED!\n");
		fprintf(stderr, "******************\n");
	}
	return (numberFailed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
