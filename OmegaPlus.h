
/*  
 *  OmegaPlus: A Parallel Tool for Rapid & Scalable Detection of 
 *	       Selective Sweeps in Genome Datasets
 *
 *  Copyright February 2012 by Nikolaos Alachiotis and Pavlos Pavlidis
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other enquiries send an email to
 *  Pavlos Pavlidis (pavlidisp@gmail.com) or
 *  Nikolaos Alachiotis (n.alachiotis@gmail.com)
 *  
 */


#ifndef _OMEGAPLUS_H
#define _OMEGAPLUS_H

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <limits.h>
#include <ctype.h>
// qLD ADDED
#include <stdint.h>
#include "gpu_kernel/blislike.h"

#ifdef _USE_PTHREADS
#include <pthread.h>
#include <unistd.h>
#endif

#ifdef _USE_OPENMP_GENERIC
#include <omp.h>
#include <unistd.h>
#endif

#define ENTRIES_PER_INT 32
#define DENOMINATOR_OFFSET 0.00001
#define REGISTER_WIDTH 32
#define MAXSIZE 10000
#define INFILENAMESIZE 1000
#define SEQNAMESIZE 1000
#define DECREASE 0.99
#define MINSNPPERWINDOW 5
#define MINSNPS_THRESHOLD 4
#define STRINGLENGTH 1000
#define VCF_HLENGTH 9 // number of fields in the VCF header line
#define MAX_CHROM_NAME_VCF 100
#define MAX_STATES_VCF 5
#define MAXAFLENGTH 9
#define RESULTS_DEFAULT 0
#define RESULTS_ALL 1

#define RSQUARE 0
#define DOM 1
#define ABSDOM 2
#define JUSTD 3
#define ABSD 4
#define ABSDOM2 5

#define ZERO '0'
#define ONE  '1'
#define GAP '-'
#define AD 'A'
#define CY 'C'
#define GU 'G'
#define TH 'T'
#define UN 'N'
#define ad 'a'
#define cy 'c'
#define gu 'g'
#define th 't'

#define BINARY 2
#define BINARY_WITH_GAPS 3
#define DNA 4
#define DNA_WITH_GAPS 5
#define STATESALL 8

#define OTHER_FORMAT -1
#define MS_FORMAT 0
#define FASTA_FORMAT 1
#define MACS_FORMAT 2
#define VCF_FORMAT 3

#ifdef _USE_PTHREADS
#define EXIT 127
#define BUSYWAIT 0
#define PAIRWISECORRELATION 1
#define COMPUTEOMEGAS 2
#endif

#ifdef _USE_OPENMP_GENERIC
#define SNP_GROUP_SIZE 32//8//128//32//256//4//128
#define SNP_GROUP_SIZE_POW2 1024//64//16384//1024//65536//16//16384
#endif

char bits_in_16bits [0x1u << 16];

char VCF_alignment_name [MAX_CHROM_NAME_VCF];
char VCF_alignment_name_cur [MAX_CHROM_NAME_VCF];

int VCF_header_lines;
int VCF_first_SNP;

int nxtVCFalignment;

int linkage_disequilibrium;

int borderTol;

char runName[INFILENAMESIZE];

typedef float cor_t;

cor_t ABS (cor_t input);

double mainTime0;

typedef struct 
{
	int              states;
	int		 length;
	int 		 segsites;
	float 	       * positions;
	int 	       * positionsInd;
	int 		 sequences;
	char 	      ** seqtable;
	char	      ** poltable;
	int		 siteSize;
	int              compressedSize;
	unsigned int  ** compressedArrays;
	cor_t         ** correlationMatrix;
	int 		 VCFsamples;
	int 	       * VCFsample_valid;

} alignment_struct;

typedef struct
{
	int 		 valid;
	float 		 omegaRealPos; // position in the alignment
	int		 omegaPos; // SNP position of the very left SNP
	int		 leftIndex; // Leftmost SNP position for specific maxw
	int 		 rightIndex; // Rightmost SNP position for specigic maxw
	int		 leftminIndex; //leftmost SNP position for specific minw
	int		 rightminIndex; //rightmost SNP position for specific minw
	float		 maxValue; // Maximum omega value
	int		 maxLeftIndex; // Leftmost SNP position for maxValue
	int		 maxRightIndex; // Rightmot SNP position for maxValue
 
} omega_struct;


#ifdef _USE_PTHREADS

pthread_t * workerThreadL;

typedef struct
{
	int limitLeft;
	int limitRight;
		
} gridPartition_t;

typedef struct
{
	alignment_struct * alignment;
	omega_struct * omega;
	int omegaIndex;
	int firstRow;
}threadArgPWC_t;

typedef struct
{
	alignment_struct * alignment;
	omega_struct * omega;
	int omegaIndex;

	float maxValue;
	int maxLeftIndex;
	int maxRightIndex;

}threadArgCO_t;

typedef struct
{
	alignment_struct * alignment;
	omega_struct * omega;
	int lvw_i;
	int cvw_i;
	int firstRowToCopy;
} threadArgSMV_t;


#ifdef _USE_PTHREADS_MULTI
typedef struct
{
	int size;
	int * coarsePartDone;
	int * myTurn;
	int * threadAv;
	int * ttid;
	int * mtid;
	int * tthreads;
	int * hold;
	int * load;
}multi_sync_t;

typedef struct
{
	alignment_struct * alignment;
	omega_struct * omega;
	int cvw_i;
	int firstRowToCompute;
	float ** myCorrelationMatrix;
	char * bits_in_16bitsLocal;
}threadArgPWC_t_MULTI;


threadArgPWC_t_MULTI * threadArgPWC_MULTI;

float * maxOmegaValueThreadsMULTI;
int * maxOmegaLeftIndexThreadsMULTI;
int * maxOmegaRightIndexThreadsMULTI;

#endif

typedef struct
{
	int threadID;
	int threadTOTAL;

	int threadBARRIER;
	int threadOPERATION;

#ifdef _USE_PTHREADS_MEMINT
	gridPartition_t * gridPartition;
	alignment_struct * alignment;
	omega_struct * omega;
	FILE * fpReport;
	int matrixSizeMax;

#ifdef _USE_PTHREADS_MULTI
	multi_sync_t * multi_sync;
#endif
#else
	threadArgPWC_t * threadArgPWC; // PairWise Correlation
	threadArgCO_t * threadArgCO; // Compute Omegas	
#endif	

}threadData_t;

#endif

cor_t ** createCorrelationMatrix(cor_t ** correlationMatrix, int matrixSize);
int findOmegaBounds (alignment_struct * alignment, omega_struct * omega, int grid, int * maxw, int minw, int minsnps);
int findNextValidOmega(omega_struct *omega,int lvw_i, int grid);
int validGridP(int cvw_i, int grid);
void checkSNIPPositions (FILE* fp, alignment_struct * alignment, int index);
void computeOmegas (alignment_struct * alignment, omega_struct * omega, int omegaIndex, void * threadData, cor_t ** correlationMatrix);
void appendOmegaResultToFile (alignment_struct * alignment, omega_struct * omega, int omegaIndex, int gridIndex, FILE * fpOut, int resultType);
void maxOmegaResultReport (float * maxomegaRealPos, float * maxomegamaxValue, int * maxomegamaxLeftIndex, int * maxomegamaxRightIndex, int grid, alignment_struct * alignment, omega_struct * omega, FILE * fpInfo);
void initializeGlobalPointers(alignment_struct* alignment);
#ifdef _SHARED
void compute_bits_in_16bits(void);
#else
void compute_bits_in_16bitsLocal(char * bits_in_16bitsLocal);
#endif
void printHeading (FILE * fp);
void computeCorrelationMatrixPairwise(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRowIndex, void * threadData, cor_t ** myCorrelationMatrix, char * lookuptable);
void applyCorrelationMatrixAdditions (omega_struct * omega, int omegaIndex, int firstRowIndex, cor_t ** correlationMatrix);
void overlapCorrelationMatrixAdditions (alignment_struct * alignment, omega_struct * omega, int lvw_i, int cvw_i, int * firstRowToCopy, int * firstRowToCompute, int * firstRowToAdd);
void shiftCorrelationMatrixValues (omega_struct * omega, int lvw_i, int cvw_i, int firstRowToCopy, cor_t ** correlationMatrix);
void commandLineParser(int argc, char** argv, char* infile, int* grid, int* length, int* minw, int* maxw, char** recfile, int* minsnps, int * imputeN, int * imputeG, int * binary, unsigned int * seed, int * fileFormat, int * threads,int * results,int * ld,int *borderTol, int *filterOut, int* noSeparator, char * samplefile_i, int * generateVCFsamplelist, int * memLimit, int * reports, double * maf, int * fileFormatMBS);
int findFirstAlignment(alignment_struct *alignment, FILE *fp, FILE *fpInfo,int format, FILE *fpVCFsamples, int generateVCFsamplelist, char * vcf_samples_filename);
int findNextAlignment(FILE *fp, int fileFormat);
void freeAlignment(alignment_struct *alignment, int matrixSizeMax);
int readAlignment(FILE *fp, alignment_struct *alignment, int imputeG, int imputeN, int binary, int format, FILE * fpInfo, int filterOut, double maf, int readAlignment);
void compressAlignment(alignment_struct *alignment);
#ifdef _USE_PTHREADS
void startThreadOperations(threadData_t * threadData, int operation);
void correlationThread(threadData_t * currentThread);
void omegasThread(threadData_t * currentThread);
#ifdef _USE_PTHREADS_MULTI
void computeCorrelationsMULTI(alignment_struct * alignment, omega_struct * omega, int cvw_i, int firstRowToCompute, float ** myCorrelationMatrixN, char * bits_in_16bitsLocal, int ttid, int tthreads);
void initmaxOmegaValueThreadsMULTI(int size);
void getmaxOmegaValueThreadsMULTI(omega_struct * omega, int cvw_i, int size);
void computeOmegasMULTI (alignment_struct * alignment, omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, int ttid, int tthreads, int tid);
#endif
#endif

#endif

#ifdef _USE_OPENMP_GENERIC


double * total_dp_init_time;
double * total_dp_update_time;
double * total_omega_values_time;


int compute_genLoad(int startIndex, int finishIndex, int prev_totalLoad, int prev_startIndex, int prev_finishIndex);
int get_SNP_groupID(int SNP_index);
float get_Mem_GroupSize();
void update_workgroup_map_ptr(float *** workgroup_map_ptr, int start, int finish, int first_group_index);
void update_workgroup_map_partial_ptr(float *** workgroup_map_ptr, int start, int finish, int prev_start, int prev_finish, int first_group_index);
void dp_on_tiles_overlap_ptr (int first_DP_tile, int last_DP_tile, float *** workgroup_map_ptr, float *** overlap_workgroup_map_ptr, int overlap, int first_group_index, alignment_struct * alignment, int leftSNPindex, int rightSNPindex, uint32_t * qLD_res);
unsigned int precomputed16_bitcount (unsigned int n);
cor_t computeCorrelationValueBIN(int sequences, unsigned int * accumXvec);
cor_t computeCorrelationValueDNA(int sequences, cor_t pairwiseCorrelationMatrix[4][4], unsigned int * valid);
void count01Combs (int total, unsigned int inputL, unsigned int inputR, unsigned int * accumXvec);
int count01GAPCombs (unsigned int inputL, unsigned int inputR, unsigned int inputLVld, unsigned int inputRVld, unsigned int * accumXvec);
void computeOmegas_generic (omega_struct * omega, int omegaIndex, float *** workgroup_map_ptr,  int first_group_index);
int max(int a, int b);
int min(int a, int b);
float computeOmega (float LS, float RS, float TS, int k, int ksel2, int m, int msel2);
#endif

// General ADDED
double gettime(void);

// qLD ADDED
enum gemm_block_sizes_e
{
    BLOCK_NC=4032,
    BLOCK_KC=256,
    BLOCK_MC=72,
    BLOCK_NR=6,
    BLOCK_MR=8
};

typedef uint64_t inputDataType_x64;
typedef uint32_t inputDataType_x32;
void Pack_A(uint32_t *A,
            unsigned int lda,
            uint32_t *A_pack,
            unsigned int m,
            unsigned int k);

void Pack_B(uint32_t *B,
            unsigned int ldb,
            uint32_t *B_pack,
            unsigned int k,
            unsigned int n);

void dgemm_ref(int k,
               int mr_alg,
               int nr_alg,
               uint32_t alpha,
               uint32_t* restrict a,
               uint32_t* restrict b,
               uint32_t* restrict c,
               int rs_c,
               int cs_c);

void gemm(unsigned int m,
          unsigned int n,
          unsigned int k,
		  uint32_t alphap,
          uint32_t * A,
          unsigned int lda,
          uint32_t * B,
          unsigned int ldb,
          uint32_t * C,
          unsigned int ldc,
          void * Ac_pack_v,
          void * Bc_pack_v);

void mlt(unsigned int m,
         unsigned int k,
         uint32_t* A,
         uint32_t* tableA);

uint32_t * correlate(uint32_t* tableA,
               int tableAsize,
               int compressed_snp_size,
			   int group_size);

// -- GENERAL OPENCL STUFF -- //
// macro define needed for deprecated warning
#define CL_TARGET_OPENCL_VERSION 110
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <CL/cl.h>

// -- qLD OPENCL STUFF -- //
#define PROGRAM_FILE "gpu_kernel/blislike.cl"
#define KERNEL_NAME "blis_like8x4"
#define RESULTS_PART_SIZE_GPU 4 //2*2
#define DOUBLE unsigned int

cl_platform_id *platforms;
cl_device_id *devices;
cl_context context;
cl_program program;
cl_mem    a_buffers[2];
cl_mem    b_buffers[2];
cl_mem    c_buffers[2];
cl_kernel kernels[2];
cl_command_queue io_queue;
cl_command_queue compute_queue;
cl_event events[4*2];
unsigned int rs_c;
// NOTE: since we're writing into a smaller buffer and not the full output
// matrix, we mult by MC to get to the next row instead of the full "m".
unsigned int cs_c;

// -- OMEGA OPENCL STUFF -- //
// #define OMEGA_NAME "omega"
// #define OMEGA_NAME "omega2"
// #define OMEGA_NAME "omega3"
#define OMEGA_NAME "omega4"
// #define OMEGA_NAME "omega5"

cl_mem omega_buffer;
cl_mem LS_buffer;
cl_mem RS_buffer;
cl_mem TS_buffer;
cl_mem k_buffer;
cl_mem m_buffer;
cl_mem LRkm_buffer;
cl_mem omega_im;
cl_mem index_im;
cl_mem index_buffer;
cl_kernel omega_kernel;

cl_uint comp_units;
size_t group_size;
cl_uint work_items;

void gpu_init(void);

void gpu_release(void);

uint32_t * correlate_gpu(uint32_t* tableA,
               int tableAsize,
               int compressed_snp_size,
			   int group_size);

void printCLErr(cl_int err,int line, char* file);

void GPU_Pack_A(inputDataType_x32 *A,
                unsigned int lda,
                DOUBLE *A_pack,
                unsigned int m,
                unsigned int k);

void GPU_Pack_B(inputDataType_x32 *B,
                unsigned int ldb,
                DOUBLE *B_pack,
                unsigned int k,
                unsigned int n);

void gpu_gemm(unsigned int m,
              unsigned int n,
              unsigned int k,
              inputDataType_x32 *A,
              unsigned int lda,
              inputDataType_x32 *B,
              unsigned int ldb,
              inputDataType_x32 * C,
              unsigned int ldc,
              void * Ac_pack_v,
              void * Bc_pack_v);

void mlt_gpu(unsigned int m,
         unsigned int k,
         inputDataType_x32 *A,
         inputDataType_x32* tableA);

void computeOmegas_gpu (alignment_struct * alignment, omega_struct * omega, int omegaIndex, void * threadData, cor_t ** correlationMatrix);

void computeOmegaValues_gpu (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData);

void computeOmega_gpu(float * omegas, float * LSs, float * RSs, float * TSs, int * ks, int * ms, int outer_cnt, int inner_cnt, unsigned int total);

void computeOmegaValues_gpu2 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData);

void computeOmega_gpu2(float * omegas, float * LSs, float * RSs, float * TSs, int * ks, int * ms, unsigned int total);
/*
void computeOmegaValues_gpu3 (omega_struct * omega, cor_t ** correlationMatrix, void * threadData, unsigned int * indexes, unsigned int cnt);

void computeOmegaValues_gpu4 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData, float * omegas, float * LSs, float * RSs, float * TSs, int * ks, int * ms);
*/
void computeOmegaValues_gpu5 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData);

void computeOmega_gpu5 (float * omegas, float * LRkm, float * TSs, int outer_cnt, int inner_cnt, unsigned int total);

void computeOmegaValues_gpu6 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData);

void computeOmega_gpu6(float * maxW, unsigned int * maxI, float * LRkm, float * TSs, int outer_cnt, int inner_cnt, unsigned int total);

void computeOmegaValues_gpu7 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData);

void computeOmega_gpu7(float * maxW, unsigned int * maxI, float * LSs, float * RSs, float * TSs, int * ks, int * ms, unsigned int total);