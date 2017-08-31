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

#ifndef COMMON_H_
#define COMMON_H_

#define	PROG_VERSION		"0.11.2"

#define CODE_A				0
#define CODE_C				1
#define	CODE_G				2
#define	CODE_T				3
#define	CODE_N				CODE_A
#define	CODE_STAR			CODE_A
#define SEED_LENGTH			12 /* Must be greater than or equal to 8. */
#define	DNA_ALPHABET_SIZE	4 /* Do not change this. */

#define CUDA_FAILURE		1
#define FAILURE				1

#define NUM_H_MATRIX_ROWS 			2
#define MAX_FILE_NAME_LENGTH		100
#define MAX_QRY_NAME_LENGTH			50
#define	MAX_REF_NAME_LENGTH			MAX_QRY_NAME_LENGTH
#define MAX_SEQ_NAME_LENGTH			MAX_QRY_NAME_LENGTH
#define MAX_QRY_SEQ_LENGTH			200 /* Do not change this. */
#define	MAX_LINE_LENGTH				100000
#define SELECTED_DEVICE_INDEX		2
#define	GLOBAL_MEM_MARGIN			100000000
#define NUM_REGS_USED				48
#define	DEFAULT_MATCH				2
#define	DEFAULT_MISMATCH			-1
#define DEFAULT_GAP_OPEN_PENALTY	-10
#define	DEFAULT_GAP_EXT_PENALTY		-1
#define TMP_DIR						"/home/pgupta/data/gpusw" /**< Path to the directory
where files will be stored temporarily. Make sure this directory already
exists. */
#define PATH_SEPARATOR				"/" /**< File path separator. */
#define FILTERED_RESULTS_FILE_NAME	"matchFile.txt"
#define	TEMP_REF_FILE				"tmpRefFile.txt"
#define SEQ_PADDING_CHAR			"*"
#define QUERY_LENGTH				100
#define	REF_LENGTH					300 /* Do not change this. */
#define NUM_BITS_IN_A_BYTE			8
#define	MAX_NUM_HITS_SINGLE_END		3 /**< Max number of hits per query for
single-end alignment. */
//#define	MAX_NUM_HITS_PAIRED_END		5 /**< Max number of hits per query for
//paired-end alignment. */
#define	MAX_NUM_HITS_PAIRED_END		3 /**< Max number of hits per query for
paired-end alignment. */
#define NUM_NEW_LINES_BUFFER		2
#define	MIN_FRAG_SIZE				99 /**< Minimum DNA fragment size used
while sequencing. */
#define MAX_FRAG_SIZE				600 /**< Maximum DNA fragment size used
while sequencing. */

#define REF_IDX_BITS				5 /**< Number of bits used to represent
reference index of a tuple. */
#define	POS_BITS					27 /**< Number of bits used to represent
tuple position within reference. */
#define	TUPLE_IGNORE_THRES			1000	/**< Tuples that have more than this
number of matching hits in the seeding step will be ignored. */

#define	REF_NUM_BASES_PER_LINE		50	/**< Number of bases per line
in the reference sequence. */

#define MAX_GENOMIC_SIZE			4000000000 /* Max number of bases
in the genome. */

//#define	REF_IDX_BITS2				8	/* Number of bits for reference index. */
#define SHIFT_BITS					28	/* Number of bits for shift. */
#define	REF_POS_BITS				28	/* Number of bits for reference position. */
#define MAX_BASES_IN_CHR			268435455	/* Max number of bases in
any human chromosome. (2^28 - 1 = 268435455) */
#define REF_IDX_SHIFT_BITS			56	/* Number of bits reference index
needs to be shifted to the right. (SHIFT_BITS + REF_POS_BITS = 56) */
#define SHIFT_SHIFT_BITS			28	/* Number of bits "shift" needs to be
shifted to the right. */
#define SHIFT_MASK					72057593769492480	/* This is a bit-mask
that will be used to extract "shift" out of the composite reference position. */
#define	REF_POS_MASK				268435455	/* This is a bit-mask that will
be used to extract "reference position" out of the composite reference position. */
#define	NUM_THREADS					16 /* Number of CPU threads to use
in Phase 1 (indexing phase). */
#define	WARP_SIZE					32

#define NUM_ASCII_CHARS				128
#define MAX_SEED_LENGTH				15
#define	MAX_NUM_REFS				25
#define	NUM_BLOCKS_PHASE1			55000

#define	CHR1_MAX_SIZE				249250621
#define	CHR2_MAX_SIZE				243199373
#define	CHR3_MAX_SIZE				198022430
#define	CHR4_MAX_SIZE				191154264
#define	CHR5_MAX_SIZE				180915260
#define	CHR6_MAX_SIZE				171115067
#define	CHR7_MAX_SIZE				159138663
#define	CHR8_MAX_SIZE				146364022
#define	CHR9_MAX_SIZE				141213431
#define	CHR10_MAX_SIZE				135534747
#define	CHR11_MAX_SIZE				135006516
#define	CHR12_MAX_SIZE				133851895
#define	CHR13_MAX_SIZE				115169878
#define	CHR14_MAX_SIZE				107349540
#define	CHR15_MAX_SIZE				102531392
#define	CHR16_MAX_SIZE				90354753
#define	CHR17_MAX_SIZE				81195210
#define CHR18_MAX_SIZE				78077248
#define	CHR19_MAX_SIZE				59128983
#define	CHR20_MAX_SIZE				63025520
#define CHR21_MAX_SIZE				48129895
#define	CHR22_MAX_SIZE				51304566
#define	CHRX_MAX_SIZE				155270550
#define CHRY_MAX_SIZE				59373566
#define CHRMT_MAX_SIZE				16569

#define CLUST_SPAN					3000
#define	CLUST_OVERLAP				MAX_QRY_SEQ_LENGTH


typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned short ushort;
typedef unsigned char uchar;

/**
 * Enumerations for output format.
 */
enum OutputFormat {DEFAULT_OUTPUT_FORMAT, SAM_OUTPUT_FORMAT};

/**
 * Enumerations for boolean values.
 */
enum Boolean {FALSE, TRUE};


/**
 * Macro to print CUDA error.
 */
#define PRINT_CUDA_ERROR() { \
							cudaError_t cudaError = cudaGetLastError(); \
							if (cudaError != cudaSuccess) \
							{ \
								fprintf(stderr, \
										"Error: %s:%s:%d: %s\n", \
										__FILE__, \
										__FUNCTION__, \
										__LINE__, \
										cudaGetErrorString(cudaError)); \
										exit(CUDA_FAILURE); \
							} \
}

#endif /* COMMON_H_ */
