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

#include "preprocess.h"
#include "input.h"
#include "common.h"
#include "search.h"
#include "search2.h"
#include "align.h"
#include <time.h>
#include "refMap.h"
#include <limits.h>


static int getParameters(int numArgs, char *args[], char **queryFileName,
		char **queryFileName2, char **refFileName, int *querySeqLength,
		int *refSeqLength, int *numQueries, int *numReferences, float *match,
		float *mismatch, float *gapOpenPenalty, float *gapExtPenalty,
		int *outputFormat, char **outputFile, int *minFragSize,
		int *maxFragSize, int *useCpuOnly, int *tupleIgnoreThres, int *seedLen,
		const char *usage);
const char *getUsage();


/**
 * Main function.
 */
int main(int argc, char *argv[])
{
	char *queryFileName;
	char *queryFileName2 = NULL;
	char *refFileName;
	int querySeqLength = -1;
	int refSeqLength = -1;
	int numQueries = -1;
	int numReferences = -1;
	float match = DEFAULT_MATCH;
	float mismatch = DEFAULT_MISMATCH;
	float gapOpenPenalty = DEFAULT_GAP_OPEN_PENALTY;
	float gapExtPenalty = DEFAULT_GAP_EXT_PENALTY;
	int outputFormat = DEFAULT_OUTPUT_FORMAT;
	char *outputFile = NULL;
	int minFragSize = MIN_FRAG_SIZE, maxFragSize = MAX_FRAG_SIZE;
	int useCpuOnly = 0;
	time_t startTime, endTime;
	double diffTime1, diffTime2;
	int tupleIgnoreThres = TUPLE_IGNORE_THRES;
	int seedLen = SEED_LENGTH;

	srand(time(NULL));
	const char *usage = getUsage();
	preprocessCreate();

	int hasParameters = getParameters(argc, argv, &queryFileName,
			&queryFileName2, &refFileName, &querySeqLength, &refSeqLength,
			&numQueries, &numReferences, &match, &mismatch, &gapOpenPenalty,
			&gapExtPenalty, &outputFormat, &outputFile, &minFragSize,
			&maxFragSize, &useCpuOnly, &tupleIgnoreThres, &seedLen, usage);
	if (hasParameters == 0)
		return 0;
	if (numQueries == -1)
		numQueries = getSeqCount(queryFileName);
//	if (numReferences == -1)
//		numReferences = getSeqCount(refFileName);
	if (querySeqLength == -1)
		querySeqLength = getSeqSize(queryFileName);
//	if (refSeqLength == -1)
//		refSeqLength = getSeqSize(refFileName);

	char filteredResultsFile[MAX_FILE_NAME_LENGTH];
	sprintf(filteredResultsFile, "%s%s%s", TMP_DIR, PATH_SEPARATOR,
			FILTERED_RESULTS_FILE_NAME);

	/* Single-end query alignment. */
	if (queryFileName2 == NULL)
	{
		int numMatches;
		if (useCpuOnly == 1)
		{
			fprintf(stderr, "Running on CPU only\n");
			fprintf(stderr, "Creating reference map...");
			time(&startTime);
			refMapCreate(refFileName);
			time(&endTime);
			diffTime1 = difftime(endTime, startTime);
			fprintf(stderr, "done. (Time = %.2lf secs)\n", diffTime1);
//			searchQueries(queryFileName, numQueries, querySeqLength,
//					refFileName, seedLen, filteredResultsFile,
//					MAX_NUM_HITS_SINGLE_END, &numMatches);
			searchQueries3(queryFileName, numQueries, querySeqLength,
					refFileName, seedLen, filteredResultsFile,
					MAX_NUM_HITS_SINGLE_END, &numMatches, NUM_THREADS,
					tupleIgnoreThres);
//			alignQueries_cpu(queryFileName, numQueries, querySeqLength,
//					refFileName, filteredResultsFile, outputFormat, match,
//					mismatch, gapOpenPenalty, gapExtPenalty, outputFile);
			alignQueries_cpu_mt(queryFileName, numQueries, querySeqLength,
					refFileName, filteredResultsFile, outputFormat, match,
					mismatch, gapOpenPenalty, gapExtPenalty, outputFile);
		}
		else
		{
			fprintf(stderr, "******************\n");
			fprintf(stderr, "Parameters\n");
			fprintf(stderr, "******************\n");
			fprintf(stderr, "Tuple ignore threshold = %d\n", tupleIgnoreThres);
			fprintf(stderr, "Seed length = %d\n", seedLen);
			fprintf(stderr, "Max number of hits per query in Phase 1 = %d\n",
					MAX_NUM_HITS_SINGLE_END);

			fprintf(stderr, "\n******************\n");
			fprintf(stderr, "PHASE 1: Seeding\n");
			fprintf(stderr, "******************\n");
			fprintf(stderr, "Creating reference map...");
			time(&startTime);
			refMapCreate(refFileName);
			time(&endTime);
			diffTime1 = difftime(endTime, startTime);
			fprintf(stderr, "done. (Time = %.2lf secs)\n", diffTime1);

			time(&startTime);
//			searchQueries3(queryFileName, numQueries, querySeqLength,
//					refFileName, seedLen, filteredResultsFile,
//					MAX_NUM_HITS_SINGLE_END, &numMatches, NUM_THREADS,
//					tupleIgnoreThres);
			search2Queries(queryFileName, refFileName, filteredResultsFile,
					seedLen, MAX_NUM_HITS_SINGLE_END, &numMatches,
					tupleIgnoreThres);
			time(&endTime);
			diffTime1 = difftime(endTime, startTime);

			time(&startTime);
			fprintf(stderr, "\n******************\n");
			fprintf(stderr, "PHASE 2: Alignment\n");
			fprintf(stderr, "******************\n");
//			numMatches = 28610637;
//			alignQueries(queryFileName, querySeqLength, refFileName,
//					filteredResultsFile, outputFormat, match, mismatch,
//					gapOpenPenalty, gapExtPenalty, outputFile, numMatches);
//			alignQueries2(queryFileName, querySeqLength, refFileName,
//					filteredResultsFile, outputFormat, match, mismatch,
//					gapOpenPenalty, gapExtPenalty, outputFile, numMatches);
			alignQueries3(queryFileName, querySeqLength, refFileName,
					filteredResultsFile, outputFormat, match, mismatch,
					gapOpenPenalty, gapExtPenalty, outputFile, numMatches);
			time(&endTime);
			diffTime2 = difftime(endTime, startTime);

			fprintf(stderr, "Deleting reference map...");
			refMapFree();
			fprintf(stderr, "done.\n");

			fprintf(stderr, "Phase 1 time = %.2lf secs\n", diffTime1);
			fprintf(stderr, "Phase 2 time = %.2lf secs\n", diffTime2);
		}
	}
	/* Paired-end read alignment. */
	else
	{
		fprintf(stderr, "Creating reference map...");
		time(&startTime);
		refMapCreate(refFileName);
		time(&endTime);
		diffTime1 = difftime(endTime, startTime);
		fprintf(stderr, "done. (Time = %.2lf secs)\n", diffTime1);

		time(&startTime);
		int numMatches;
		fprintf(stderr, "******************\n");
		fprintf(stderr, "PHASE 1: Seeding\n");
		fprintf(stderr, "******************\n");
		searchQueries_paired2(queryFileName, queryFileName2, numQueries,
				querySeqLength, refFileName, seedLen, filteredResultsFile,
				MAX_NUM_HITS_PAIRED_END, &numMatches);
		time(&endTime);
		diffTime1 = difftime(endTime, startTime);

		time(&startTime);
		fprintf(stderr, "\n******************\n");
		fprintf(stderr, "PHASE 2: Alignment\n");
		fprintf(stderr, "******************\n");

		if (useCpuOnly == 1)
		{
			fprintf(stderr, "Running on CPU only\n");
			alignQueries_paired_cpu(queryFileName, queryFileName2, numQueries,
					querySeqLength, refFileName, numReferences, refSeqLength,
					filteredResultsFile, outputFormat, match, mismatch,
					gapOpenPenalty, gapExtPenalty, outputFile);
		}
		else
			alignQueries_paired2(queryFileName, queryFileName2, querySeqLength,
					refFileName, filteredResultsFile, outputFormat, match,
					mismatch, gapOpenPenalty, gapExtPenalty, outputFile,
					minFragSize, maxFragSize, numMatches);
		time(&endTime);
		diffTime2 = difftime(endTime, startTime);

		fprintf(stderr, "Deleting reference map...");
		refMapFree();
		fprintf(stderr, "done.\n");

		fprintf(stderr, "Phase 1 time = %.2lf secs\n", diffTime1);
		fprintf(stderr, "Phase 2 time = %.2lf secs\n", diffTime2);
	}
	preprocessDelete();

	return 0;
}


/**
 * Fetch parameters from the command-line. Return 1 on success and 0 on
 * failure.
 *
 * @param numArgs Number of arguments passed.
 * @param args Argument array.
 * @param queryFileName Query sequence file.
 * @param queryFileName2 Second query file (for paired-read alignment).
 * @param refFileName Reference sequence file.
 * @param querySeqLength Max query sequence length.
 * @param refSeqLength Max reference sequence length.
 * @param numQueries Number of queries in the query sequence file.
 * @param numReferences Number of references in the reference sequence file.
 * @param match Match score.
 * @param mismatch Mismatch score.
 * @param gapOpenPenalty Gap opening penalty.
 * @param gapExtPenalty Gap extension penalty.
 * @param outputFormat Output format.
 * @param outputFile Output file.
 * @param minFragSize Minimum DNA fragment size used in sequencing.
 * @param maxFragSize Maximum DNA fragment size used in sequencing.
 * @param useCpuOnly Indicates whether the program should be run on a CPU only.
 * A value of '1' indicates that the program must be run on CPU only.
 * @param tupleIgnoreThres Threshold value used to ignore reference tuples.
 * @param seedLen Seed length.
 * @param usage String containing help instructions.
 * @return Returns 1 on success and 0 on failure.
 */
static int getParameters(int numArgs, char *args[], char **queryFileName,
		char **queryFileName2, char **refFileName, int *querySeqLength,
		int *refSeqLength, int *numQueries, int *numReferences, float *match,
		float *mismatch, float *gapOpenPenalty, float *gapExtPenalty,
		int *outputFormat, char **outputFile, int *minFragSize,
		int *maxFragSize, int *useCpuOnly, int *tupleIgnoreThres, int *seedLen,
		const char *usage)
{
	if (numArgs == 1)
	{
		printf("%s\n", usage);
		return 0;
	}

	/* Extract command-line options */
	int i;
	for (i = 1; i < numArgs; ++i)
	{
		if (strcmp(args[i], "-q") == 0)
		{
			*queryFileName = args[i + 1];
			++i;
		}
		else if (strcmp(args[i], "-q2") == 0)
		{
			*queryFileName2 = args[i + 1];
			++i;
		}
		else if (strcmp(args[i], "-r") == 0)
		{
			*refFileName = args[i + 1];
			++i;
		}
		else if (strcmp(args[i], "-s") == 0)
		{
			*querySeqLength = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-n") == 0)
		{
			*numQueries = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-S") == 0)
		{
			*refSeqLength = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-N") == 0)
		{
			*numReferences = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-m") == 0)
		{
			*match = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-M") == 0)
		{
			*mismatch = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-O") == 0)
		{
			*gapOpenPenalty = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-E") == 0)
		{
			*gapExtPenalty = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-f") == 0)
		{
			*outputFormat = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-o") == 0)
		{
			*outputFile = args[i + 1];
			++i;
		}
		else if (strcmp(args[i], "-ms") == 0)
		{
			*minFragSize = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-MS") == 0)
		{
			*maxFragSize = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-cpu") == 0)
		{
			*useCpuOnly = 1;
		}
		else if (strcmp(args[i], "-t") == 0)
		{
			*tupleIgnoreThres = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-l") == 0)
		{
			*seedLen = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-v") == 0)
		{
			printf("Version: %s\n", PROG_VERSION);
			return 0;
		}
		else if (strcmp(args[i], "-h") == 0 || strcmp(args[i], "--help") == 0)
		{
			printf("%s\n", usage);
			return 0;
		}
		else
		{
			printf("%s\n", usage);
			return 0;
		}
	}

	/* If reference, query, or output files are not entered,
	 * print usage. */
	if (*queryFileName == NULL || *refFileName == NULL || *outputFile == NULL)
	{
		printf("%s\n", usage);
		return 0;
	}

	return 1;
}


/**
 * Returns the program usage string.
 *
 * @return Pointer to a character string where the program usage will be
 * stored.
 */
const char *getUsage()
{
	const char *usage = "\nUSAGE:\n"
		"   /path/to/swift -q <query fasta file> -r <reference fasta file> "
		"-o <output file> [optional parameters]\n\n"
		"DESCRIPTION:\n"
		"   Aligns multiple query sequences to multiple reference sequences "
		"using the Smith-Waterman algorithm.\n\n"
		"PARAMETERS:\n"
		"   -q     Query fasta file (required)\n"
		"   -q2    Paired query fasta file\n"
		"   -r     Reference fasta file (required)\n"
		"   -o     Output file (required)\n"
		"   -s     Query sequence size\n"
		"   -n     Number of queries\n"
		"   -l     Length of a seed (Default: 12)\n"
		"   -t     Threshold value used to ignore reference tuples\n"
		"   -m     Match score (Default: 2)\n"
		"   -M     Mismatch score (Default: -1)\n"
		"   -O     Gap open penalty (Default: -10)\n"
		"   -E     Gap extension penalty (Default: -1)\n"
		"   -cpu   Run program on CPU only\n"
		"   -f     Output format\n"
		"          0 - Default program output format. "
		"Output includes alignment, score, positions, and "
		"length.\n"
		"          1 - SAM format\n"
		"   -v     Print program version\n"
		"   -h     Print usage\n\n"
		"AUTHOR:\n"
		"   Pankaj Gupta (pankaj.gupta@stjude.org)\n"
		"   John Obenauer (john.obenauer@stjude.org)\n\n"
		"ORGANIZATION:\n"
		"   St. Jude Children's Research Hospital (http://www.stjude.org)\n\n";
	return usage;
}
