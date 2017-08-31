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

#include "testAlignOutput.h"
#include <stdio.h>
#include "../src/common.h"
#include <ctype.h>


static FILE *_refFilePtr;
static FILE *_qryFilePtr;
static FILE *_alignFilePtr;
static FILE *_outputFilePtr;
static char _line[MAX_LINE_LENGTH];
static char _refName[MAX_REF_NAME_LENGTH], _qryName[MAX_QRY_NAME_LENGTH];
static char _refSeq[MAX_LINE_LENGTH], _qrySeq[MAX_LINE_LENGTH];
static float _alignScore;
static long _alignStart, _alignEnd;
static int _alignLen;
static char _message[1000];


/**
 * Main function.
 */
int main(int argc, char *argv[])
{
	char *refFile = NULL, *qryFile = NULL, *alignFile = NULL, *outputFile = NULL;

	if (argc == 1)
	{
		fprintf(stderr, "%s", getUsage());
		return 0;
	}

	int i;
	for (i = 1; i < argc; ++i)
	{
		if (strcmp(argv[i], "-r") == 0)
		{
			refFile = argv[i + 1];
			++i;
		}
		else if (strcmp(argv[i], "-q") == 0)
		{
			qryFile = argv[i + 1];
			++i;
		}
		else if (strcmp(argv[i], "-a") == 0)
		{
			alignFile = argv[i + 1];
			++i;
		}
		else if (strcmp(argv[i], "-o") == 0)
		{
			outputFile = argv[i + 1];
			++i;
		}
		else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
		{
			fprintf(stderr, "%s", getUsage());
			return 0;
		}
		else
		{
			fprintf(stderr, "%s", getUsage());
			return 0;
		}
	}

	if (refFile == NULL || qryFile == NULL || alignFile == NULL
			|| outputFile == NULL)
	{
		fprintf(stderr, "%s", getUsage());
		return 0;
	}

	validateOutput(refFile, qryFile, alignFile, outputFile);
}


/**
 * Validates the alignment output file (default output format). This function
 * does the following things in the alignment file:
 * (1) Make sure score is correct.
 * (2) Make sure alignment start position is correct.
 * (3) Make sure alignment end position is correct.
 * (4) Make sure alignment length is correct.
 * (5) Make sure the reference sequence is indeed between the specified start
 * and end positions.
 *
 * Other things that are important to check but are not currently
 * checked by this function are:
 * (1) The match in the output file is indeed the best match and there are
 * no better matches that have not been included in the output file.
 * (2) Extending the reference and the query sequence on either side
 * will not increase the alignment score.
 *
 * @param	refFile		Reference file.
 * @param	qryFile		Query file.
 * @param	alignFile	File containing alignment output.
 * @param	outputFile	Output file where error messages will be printed out.
 */
static void validateOutput(const char *refFile, const char *qryFile,
		const char *alignFile, const char *outputFile)
{
	_refFilePtr = fopen(refFile, "r");
	_qryFilePtr = fopen(qryFile, "r");
	_alignFilePtr = fopen(alignFile, "r");
	_outputFilePtr = fopen(outputFile, "w");

	while (fgets(_line, MAX_LINE_LENGTH, _alignFilePtr) != NULL)
	{
		if (_line[0] == '\n')
			continue;
		else if (strncmp(_line, "*****", 5) == 0)
			continue;
		else if (strncmp(_line, "Qry:", 4) == 0)
		{
			_message[0] = '\0';
			validateMatchInfo();
		}
		else if (strncmp(_line, "R:", 2) == 0)
		{
			validateRefSeq();
		}
		else if (strncmp(_line, "Q:", 2) == 0)
		{
			validateQrySeq();
			if (_message[0] != '\0')
			{
				fprintf(_outputFilePtr, "\n********************************\n");
				fprintf(_outputFilePtr, "%s vs %s:\n", _qryName, _refName);
				fprintf(_outputFilePtr, "%s", _message);
			}
		}
	}

	fclose(_refFilePtr);
	fclose(_qryFilePtr);
	fclose(_alignFilePtr);
	fclose(_outputFilePtr);
}


/**
 * Validates match information.
 */
static void validateMatchInfo()
{
	sscanf(_line, "%*s %s %*s %*s %s %*s %*s %f %*s %*s %ld %*s %*s %ld %*s "
			"%*s %d", _qryName, _refName, &_alignScore, &_alignStart,
			&_alignEnd, &_alignLen);

	/* Make sure alignment score is positive. */
	if (_alignScore < 0.0)
		strcat(_message, "- Score is negative.\n");

	/* Make sure alignment score is not greater than the max possible score. */
	static float maxScore;
	maxScore = _alignLen * DEFAULT_MATCH;
	if (_alignScore > maxScore)
		strcat(_message, "- Score is greater than max possible score.\n");

	/* Make sure alignment start position is positive. */
	if (_alignStart < 0)
		strcat(_message, "- Align start position is negative.\n");

	/* Make sure alignment end position is positive. */
	if (_alignEnd < 0)
		strcat(_message, "- Align end position is negative.\n");

	/* Make sure alignment start position is less than or equal to alignment
	 * end position. */
	if (_alignStart > _alignEnd)
		strcat(_message, "- Align start position is greater than "
				"alignment end position.\n");
}


/**
 * Validates reference sequence.
 */
static void validateRefSeq()
{
	sscanf(_line, "%*s %s", _refSeq);

	/* Compare reference sequence length to alignment length. */
	if (strlen(_refSeq) != _alignLen)
		strcat(_message, "- Reference sequence length is not equal "
				"to alignment length.\n");

	/* Make sure the reference sequence is indeed between the start and
	 * end positions. */
	static char refSeq[MAX_LINE_LENGTH];
	getRefSeq(refSeq);
	static char str[1000];
	static int i, j;
	j = 0;
	for (i = 0; i < _alignLen; ++i)
	{
		if (_refSeq[i] == '-')
			continue;
		else if (_refSeq[i] != refSeq[j])
		{
			strcat(_message, "- Reference is not between the start and "
					"end positions.\n");
			sprintf(str, "start = %ld, end = %ld\n", _alignStart, _alignEnd);
			strcat(_message, str);
			strcat(_message, "_rf: ");
			strcat(_message, _refSeq);
			strcat(_message, "\nref: ");
			strcat(_message, refSeq);
			strcat(_message, "\n");
			break;
		}
		else
			++j;
	}
}


/**
 * Extract reference sequence from the reference sequence file between
 * @a _alignStart and @a _alignEnd positions.
 *
 * @param[out]	refSeq	Reference sequence that was extracted from the
 * reference sequence file.
 */
static void getRefSeq(char *refSeq)
{
	static char line[MAX_LINE_LENGTH], refName[MAX_REF_NAME_LENGTH];
	static int isCorrectRef, isFragFound, lineLength, sum, skip;
	static int seqLen, refSeqLen;
	seqLen = _alignEnd - _alignStart + 1;
	isCorrectRef = FALSE;
	isFragFound = FALSE;
	sum = 0;
	refSeqLen = 0;
	rewind(_refFilePtr);
	while (fgets(line, MAX_LINE_LENGTH, _refFilePtr) != NULL)
	{
		lineLength = strlen(line);
		if (line[lineLength - 1] == '\n')
		{
			line[lineLength - 1] = '\0';
			--lineLength;
		}

		if (line[0] == '>')
		{
			strcpy(refName, ((char *) &line) + 1);
			if (strcmp(refName, _refName) == 0)
				isCorrectRef = TRUE;
		}
		else if (isCorrectRef == TRUE)
		{
			sum += lineLength;
			if (sum >= _alignStart && isFragFound == FALSE)
			{
				skip = lineLength - (sum - _alignStart) - 1;
				strcpy(refSeq, ((char *) &line) + skip);
				refSeqLen = lineLength - skip;
				isFragFound = TRUE;
			}
			else if (isFragFound == TRUE)
			{
				strcat(refSeq, line);
				refSeqLen += lineLength;
			}

			if (refSeqLen > seqLen)
			{
				refSeq[refSeqLen - (refSeqLen - seqLen)] = '\0';
				break;
			}
		}
	}
}


/**
 * Validates query sequence.
 */
static void validateQrySeq()
{
	sscanf(_line, "%*s %s", _qrySeq);

	/* Compare query sequence length to alignment length. */
	if (strlen(_qrySeq) != _alignLen)
		strcat(_message, "- Query sequence length is not equal "
				"to alignment length.\n");

	/* Compare alignment score with the score computed manually by comparing
	 * the reference and query strings. */
	static int i, isRefFirstGap, isQryFirstGap;
	static float sum;
	sum = 0.0f;
	isRefFirstGap = TRUE;
	isQryFirstGap = TRUE;
	for (i = 0; i < _alignLen; ++i)
	{
		if (_refSeq[i] == '-' && _qrySeq[i] == '-')
		{
			strcat(_message, "- Both reference and query have a gap "
					"at the same position.\n");
		}
//		else if (toupper(_refSeq[i]) == toupper(_qrySeq[i]))
		else if (toupper(_refSeq[i]) == toupper(_qrySeq[i])
				|| _refSeq[i] == 'N' || _qrySeq[i] == 'N')
		{
			sum += DEFAULT_MATCH;
			isRefFirstGap = TRUE;
			isQryFirstGap = TRUE;
		}
		else if (_refSeq[i] == '-')
		{
			if (isRefFirstGap == TRUE)
			{
				sum += DEFAULT_GAP_OPEN_PENALTY;
				isRefFirstGap = FALSE;
			}
			else
				sum += DEFAULT_GAP_EXT_PENALTY;
		}
		else if (_qrySeq[i] == '-')
		{
			if (isQryFirstGap == TRUE)
			{
				sum += DEFAULT_GAP_OPEN_PENALTY;
				isQryFirstGap = FALSE;
			}
			else
				sum += DEFAULT_GAP_EXT_PENALTY;
		}
		else
		{
			sum += DEFAULT_MISMATCH;
			isRefFirstGap = TRUE;
			isQryFirstGap = TRUE;
		}
	}
	if (sum != _alignScore)
		strcat(_message, "- Alignment score does not match "
				"the score obtained manually by comparing reference and "
				"query strings.\n");
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
		"   /path/to/testOutput -q <query fasta file> -r <reference fasta file> "
		"-a <alignment file> -o <output file>\n\n"
		"DESCRIPTION:\n"
		"   Tests the output alignment file for the following things:\n"
		"   (1) Make sure score is correct.\n"
		"   (2) Make sure alignment start position is correct.\n"
		"   (3) Make sure alignment end position is correct.\n"
		"   (4) Make sure alignment length is correct.\n"
		"   (5) Make sure the reference sequence is indeed between the \n"
		"       specified start and end positions.\n\n"
		"   Other things that are important to check but are not currently \n"
		"   checked by this function are:\n"
		"   (1) The match in the output file is indeed the best match and \n"
		"       there are no better matches that have not been included in \n"
		"       the output file.\n"
		"   (2) Extending the reference and the query sequence on either \n"
		"       side will not increase the alignment score.\n\n"
		"PARAMETERS:\n"
		"   -q     Query fasta file (required)\n"
		"   -r     Reference fasta file (required)\n"
		"   -a     Alignment file (required)\n"
		"   -o     Output file where error messages will be output (required)\n"
		"   -h     Print usage\n\n"
		"AUTHOR:\n"
		"   Pankaj Gupta (pankaj.gupta@stjude.org)\n"
		"   John Obenauer (john.obenauer@stjude.org)\n\n"
		"ORGANIZATION:\n"
		"   St. Jude Children's Research Hospital (http://www.stjude.org)\n\n";
	return usage;
}
