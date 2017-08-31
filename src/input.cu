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
#include "common.h"
#include "input.h"
#include "assert.h"

char line[MAX_LINE_LENGTH];
int lineLength;
int newSeqFlag;
int charCount;


/**
 * Fetches the sequence name from the file pointed to by the given file
 * descriptor and saves the string in the given character pointer.
 *
 * @param filePtr File descriptor of the sequence file
 * @param fileOffset File offset from where parsing should begin
 * @param name Character pointer where the sequence name should be stored
 */
void getSeqName(FILE *filePtr, ulong fileOffset, char *name)
{
	if (fseek(filePtr, fileOffset, SEEK_SET) != 0)
	{
		fprintf(stderr, "Error: %s:%s:%d: fseek has error\n", __FILE__,
				__FUNCTION__, __LINE__);
		fprintf(stderr, "file offset = %lu\n", fileOffset);
		exit(FAILURE);
	}
	fgets(line, MAX_LINE_LENGTH, filePtr);
	lineLength = strlen(line);
	if (line[lineLength - 1] == '\n')
		line[lineLength - 1] = '\0';
	strcpy(name, ((char *) &line) + 1);
}


/**
 * Fetches offsets stored in the file pointed to by the given file pointer.
 *
 * @param filePtr File descriptor of the file containing name and sequence
 * offsets of queries and references. This must not be NULL.
 * @param[out] queryNameOffset Long integer where query name offset will
 * be stored. This must not be NULL.
 * @param[out] querySeqOffset Long integer where query sequence offset will
 * be stored. This must not be NULL.
 * @param[out] refNameOffset Long integer where reference name offset will
 * be stored. This must not be NULL.
 * @param[out] refSeqOffset Long integer where reference sequence offset will
 * be stored. This must not be NULL.
 * @param[out] refSeqFragOffset Long integer where reference sequence
 * fragment offset will be stored. This must not be NULL.
 * @param[out] refDistance Integer where distance from the beginning of
 * the reference will be stored. This must not be NULL.
 * @param[out] isReverseComplement Indicates whether the offsets belong to
 * the reverse complement of a query. This must not be NULL.
 * @param[out] refIdx Reference index.
 */
void getOffsets(FILE *filePtr, uint *queryNameOffset, uint *querySeqOffset,
		uint *refNameOffset, uint *refSeqOffset, uint *refSeqFragOffset,
		int *refDistance, uint *isReverseComplement, int *refIdx)
{
	fgets(line, MAX_LINE_LENGTH, filePtr);
	sscanf(line, "%*s\t%d\t%d\t%*d\t%d\t%u\t%u", isReverseComplement, refIdx,
			refDistance, queryNameOffset, querySeqOffset);
}


/**
 * Fetches offsets stored in the file pointed to by the given file pointer.
 *
 * @param		filePtr			File descriptor of the file containing name
 * and sequence offsets of queries and references. This must not be NULL.
 * @param[out]	queryNameOffset	Long integer where query name offset will
 * be stored. This must not be NULL.
 * @param[out]	querySeqOffset	Long integer where query sequence offset will
 * be stored. This must not be NULL.
 * @param[out]	mateIdx			Variable where mate index will be stored.
 * @param[out]	refNameOffset	Long integer where reference name offset will
 * be stored. This must not be NULL.
 * @param[out]	refSeqOffset	Long integer where reference sequence offset will
 * be stored. This must not be NULL.
 * @param[out]	refSeqFragOffset	Long integer where reference sequence
 * fragment offset will be stored. This must not be NULL.
 * @param[out]	refDistance		Integer where distance from the beginning of
 * the reference will be stored. This must not be NULL.
 * @param[out]	isRevComp		Indicates whether the offsets belong to
 * the reverse complement of a query. This must not be NULL.
 * @param[out]  refIdx          Variable where reference index will be stored.
 */
void getOffsets_paired(FILE *filePtr, ulong *queryNameOffset,
		ulong *querySeqOffset, int *mateIdx, ulong *refNameOffset,
		ulong *refSeqOffset, ulong *refSeqFragOffset, int *refDistance,
		uint *isRevComp, int *refIdx)
{
    fgets(line, MAX_LINE_LENGTH, filePtr);
    sscanf(line, "%*d\t%*s\t%d\t%d\t%d\t%*d\t%d\t%lu\t%lu",
            isRevComp, mateIdx, refIdx, refDistance, queryNameOffset,
            querySeqOffset);
}


/**
 * Fetches offsets stored in the file pointed to by the given file pointer.
 *
 * @param       filePtr         File descriptor of the file containing name
 * and sequence offsets of queries and references. This must not be NULL.
 * @param[out]  queryNameOffset Long integer where query name offset will
 * be stored. This must not be NULL.
 * @param[out]  querySeqOffset  Long integer where query sequence offset will
 * be stored. This must not be NULL.
 * @param[out]  mateIdx         Variable where mate index will be stored.
 * @param[out]  refNameOffset   Long integer where reference name offset will
 * be stored. This must not be NULL.
 * @param[out]  refSeqOffset    Long integer where reference sequence offset will
 * be stored. This must not be NULL.
 * @param[out]  refSeqFragOffset    Long integer where reference sequence
 * fragment offset will be stored. This must not be NULL.
 * @param[out]  refDistance     Integer where distance from the beginning of
 * the reference will be stored. This must not be NULL.
 * @param[out]  isRevComp       Indicates whether the offsets belong to
 * the reverse complement of a query. This must not be NULL.
 * @param[out]  refIdx          Variable where reference index will be stored.
 * @param[out]  qryIdx          Variable where query index will be stored.
 */
void getOffsets_paired2(FILE *filePtr, uint *queryNameOffset,
        uint *querySeqOffset, int *mateIdx, uint *refNameOffset,
        uint *refSeqOffset, uint *refSeqFragOffset, int *refDistance,
        uint *isRevComp, int *refIdx, int *qryIdx)
{
    fgets(line, MAX_LINE_LENGTH, filePtr);
    sscanf(line, "%d\t%*s\t%d\t%d\t%d\t%*d\t%d\t%u\t%u",
            qryIdx, isRevComp, mateIdx, refIdx, refDistance, queryNameOffset,
            querySeqOffset);
}


/**
 * Fetch sequence from the file pointed to by the given file descriptor
 * and store the sequence in the given character pointer.
 *
 * @param filePtr File descriptor of the sequence file
 * @param fileOffset File offset from where parsing should begin
 * @param seq Character pointer where the parsed sequence should be stored
 * @param maxLength Length of the sequence that should be parsed
 */
void getSeq(FILE *filePtr, ulong fileOffset, char *seq, int maxLength)
{
	if (fseek(filePtr, fileOffset, SEEK_SET) != 0)
	{
		fprintf(stderr, "Error: %s:%s:%d: fseek has error\n", __FILE__,
				__FUNCTION__, __LINE__);
		fprintf(stderr, "file offset = %lu\n", fileOffset);
		exit(FAILURE);
	}
	newSeqFlag = 1;
	charCount = 0;
	while (fgets(line, MAX_LINE_LENGTH, filePtr) != NULL)
	{
		lineLength = strlen(line);
		if (line[0] == '>')
			break;
		else if (line[lineLength - 1] == '\n')
		{
			line[lineLength - 1] = '\0';
			--lineLength;
		}
		charCount += lineLength;
		if (charCount > maxLength)
		{
			lineLength -= (charCount - maxLength);
			line[lineLength] = '\0';
		}

		if (newSeqFlag == 1)
		{
			strcpy(seq, ((char *) &line));
			newSeqFlag = 0;
		}
		else
			strcat(seq, line);

		/* If parsed sequence length has gone beyond the max length, then
		 * break out of the loop */
		if (charCount >= maxLength)
			break;
	}
}


/**
 * Fetch sequence from the file pointed to by the given file descriptor
 * and store the sequence in the given character pointer.
 *
 * @param filePtr File descriptor of the sequence file
 * @param fileOffset File offset from where parsing should begin
 * @param seq Character pointer where the parsed sequence should be stored
 * @param maxLength Length of the sequence that should be parsed
 */
void getSeq2(FILE *filePtr, ulong fileOffset, char *seq, int maxLength)
{
	if (fseek(filePtr, fileOffset, SEEK_SET) != 0)
	{
		fprintf(stderr, "Error: %s:%s:%d: fseek has error\n", __FILE__,
				__FUNCTION__, __LINE__);
		fprintf(stderr, "file offset = %lu\n", fileOffset);
		exit(FAILURE);
	}
	newSeqFlag = 1;
	charCount = 0;
	while (fgets(line, MAX_LINE_LENGTH, filePtr) != NULL)
	{
		lineLength = strlen(line);
		if (line[0] == '>')
			break;
		else if (line[lineLength - 1] == '\n')
		{
			line[lineLength - 1] = '\0';
			--lineLength;
		}
		charCount += lineLength;
		if (charCount > maxLength)
		{
			lineLength -= (charCount - maxLength);
			line[lineLength] = '\0';
		}

		if (newSeqFlag == 1)
		{
			strcpy(seq, ((char *) &line));
			newSeqFlag = 0;
		}
		else
			strcat(seq, line);

		/* If parsed sequence length has gone beyond the max length, then
		 * break out of the loop */
		if (charCount >= maxLength)
			break;
	}
	strcat(seq, "\0");
}


/**
 * Fetch sequence from the file pointed to by the given file descriptor
 * and store the sequence in the given character pointer.
 *
 * @param filePtr File descriptor of the sequence file
 * @param fileOffset File offset from where parsing should begin
 * @param seq Character pointer where the parsed sequence should be stored
 * @param maxLength Length of the sequence that should be parsed
 */
void getRefSeq(FILE *filePtr, ulong fileOffset, char *seq, int maxLength)
{
	if (fseek(filePtr, fileOffset, SEEK_SET) != 0)
	{
		fprintf(stderr, "Error: %s:%s:%d: fseek has error\n", __FILE__,
				__FUNCTION__, __LINE__);
		fprintf(stderr, "file offset = %lu\n", fileOffset);
		exit(FAILURE);
	}
	newSeqFlag = 1;
	charCount = 0;
	char isNFine;
	char *line2;
	isNFine = 0;
	while (fgets(line, MAX_LINE_LENGTH, filePtr) != NULL)
	{
		lineLength = strlen(line);
		if (line[0] == '>')
			break;
		else if (line[lineLength - 1] == '\n')
		{
			line[lineLength - 1] = '\0';
			--lineLength;
		}
		charCount += lineLength;
		if (charCount > maxLength)
		{
			lineLength -= (charCount - maxLength);
			line[lineLength] = '\0';
		}

		/* If there are any 'n's in the beginning, skip them. */
		line2 = line;
		if (isNFine == 0)
		{
			while (lineLength > 0)
			{
				if (line2[0] == 'n' || line2[0] == 'N')
				{
					++line2;
					--lineLength;
				}
				else
				{
					isNFine = 1;
					break;
				}
			}
		}

		if (lineLength > 0)
		{
			if (newSeqFlag == 1)
			{
				strcpy(seq, line2);
				newSeqFlag = 0;
			}
			else
				strcat(seq, line2);
		}

		/* If parsed sequence length has gone beyond the max length, then
		 * break out of the loop */
		if (charCount >= maxLength)
			break;
	}
}


/**
 * Fetch sequence from the file pointed to by the given file descriptor
 * and store the sequence in the given character pointer.
 *
 * @param filePtr File descriptor of the sequence file
 * @param fileOffset File offset from where parsing should begin
 * @param seq Character pointer where the parsed sequence should be stored
 * @param maxLength Length of the sequence that should be parsed
 */
void getRefSeq2(FILE *filePtr, ulong fileOffset, char *seq, int maxLength)
{
	if (fseek(filePtr, fileOffset, SEEK_SET) != 0)
	{
		fprintf(stderr, "Error: %s:%s:%d: fseek has error\n", __FILE__,
				__FUNCTION__, __LINE__);
		fprintf(stderr, "file offset = %lu\n", fileOffset);
		exit(FAILURE);
	}
	newSeqFlag = 1;
	charCount = 0;
	char tmpSeq[maxLength + 1];
	while (fgets(line, MAX_LINE_LENGTH, filePtr) != NULL)
	{
		lineLength = strlen(line);
		if (line[0] == '>')
			break;
		else if (line[lineLength - 1] == '\n')
		{
			line[lineLength - 1] = '\0';
			--lineLength;
		}
		charCount += lineLength;
		if (charCount > maxLength)
		{
			lineLength -= (charCount - maxLength);
			line[lineLength] = '\0';
		}

		if (lineLength > 0)
		{
			if (newSeqFlag == 1)
			{
				strcpy(tmpSeq, line);
				newSeqFlag = 0;
			}
			else
				strcat(tmpSeq, line);
		}

		/* If parsed sequence length has gone beyond the max length, then
		 * break out of the loop */
		if (charCount >= maxLength)
			break;
	}

	/* Remove 'n's from the beginning and the end. */
	char *tmpSeq2 = tmpSeq;
	while (tmpSeq2[0] == 'n' || tmpSeq2[0] == 'N')
			++tmpSeq2;
	int len = strlen(tmpSeq);
	while (tmpSeq2[len - 1] == 'n' || tmpSeq2[len - 1] == 'N')
		--len;
	strncpy(seq, tmpSeq2, len);
	seq[len] = '\0';
}


/**
 * Fetch sequence from the file pointed to by the given file descriptor
 * and store the sequence in the given character pointer.
 *
 * @param			filePtr		File descriptor of the sequence file
 * @param[in,out]	fileOffset	File offset from where parsing should begin
 * @param[in,out]	distance	Position of the first base of the sequence
 * fragment.
 * @param			shift		Number of positions the sequence should be
 * shifted.
 * @param			seq			Character pointer where the parsed sequence
 * should be stored
 * @param			maxLength	Length of the sequence that should be parsed
 */
//void getRefSeq(FILE *filePtr, ulong *fileOffset, int *distance, int shift,
//		char *seq, int maxLength)
//{
//	if (fseek(filePtr, fileOffset, SEEK_SET) != 0)
//	{
//		fprintf(stderr, "Error: %s:%s:%d: fseek has error\n", __FILE__,
//				__FUNCTION__, __LINE__);
//		fprintf(stderr, "file offset = %lu\n", fileOffset);
//		exit(FAILURE);
//	}
//	newSeqFlag = 1;
//	charCount = 0;
//	while (fgets(line, MAX_LINE_LENGTH, filePtr) != NULL)
//	{
//		lineLength = strlen(line);
//		if (line[0] == '>')
//			break;
//		else if (line[lineLength - 1] == '\n')
//		{
//			line[lineLength - 1] = '\0';
//			--lineLength;
//		}
//		charCount += lineLength;
//		if (charCount > maxLength)
//		{
//			lineLength -= (charCount - maxLength);
//			line[lineLength] = '\0';
//		}
//
//		if (newSeqFlag == 1)
//		{
//			strcpy(seq, ((char *) &line));
//			newSeqFlag = 0;
//		}
//		else
//			strcat(seq, line);
//
//		/* If parsed sequence length has gone beyond the max length, then
//		 * break out of the loop */
//		if (charCount >= maxLength)
//			break;
//	}
//}


/**
 * Returns the size of each sequence in the given file.
 *
 * @note The size of the first sequence is assumed to be the size of every
 * sequence.
 *
 * @param fileName Pointer to the file name
 * @return Size of each sequence. The size of the first sequence is assumed
 * to be the size of every sequence.
 */
int getSeqSize(char *fileName)
{
	FILE *filePtr = fopen(fileName, "r");
	if (filePtr == NULL)
	{
		fprintf(stderr, "Error: Could not open file %s\n", fileName);
		return 0;
	}

	int nameLineCount = 0, seqSize = 0;
	while (fgets(line, MAX_LINE_LENGTH, filePtr) != NULL)
	{
		if (line[0] == '>' && nameLineCount >= 1)
			break;
		else if (line[0] == '>' && nameLineCount < 1)
			nameLineCount++;
		else /* Sequence line */
			seqSize += strlen(line) - 1;
	}
	fclose(filePtr);
	return seqSize;
}


/**
 * Returns the number of sequences in the given file
 *
 * @param fileName Pointer to the file name
 * @return Number of sequences in the file
 */
int getSeqCount(char *fileName)
{
	FILE *filePtr = fopen(fileName, "r");
	if (filePtr == NULL)
	{
		fprintf(stderr, "Error: Could not open file %s\n", fileName);
		return 0;
	}

	int maxLength = 2; /* We are only concerned with the first character */
	char line[maxLength];
	int count = 0;
	while (fgets(line, 2, filePtr) != NULL)
	{
		if (line[0] == '>')
			++count;
	}
	fclose(filePtr);
	return count;
}
