
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


#include "OmegaPlus.h"

inline int min(int a, int b)
{
  if(a < b)
    return a;
  return b;
}

inline int max(int a, int b)
{
  if(a > b)
    return a;
  return b;
}

int moveright(int start, int* intsnps, int snpsize, int pos)
{
	int i = start;

	while(i<snpsize && intsnps[i] <= pos)
		++i;

	return i-1;
}

int moveleft(int start, int* intsnps, int snpsize, int pos)
{
	int i = start;

	while(i>=0 && intsnps[i] >= pos)
		--i;

	return i+1;
}

void boundaries(omega_struct* omstruct, 
		int grid, 
		int* intsnps, 
		int snpsize, 
		int minw, 
		int minsnps, /* the minimum number of snps a window is allowed to have (usually 2)*/
		int maxw, 
		int *changemaxw, 
		int *maxSizeMatrix)
{
	int startsnp = intsnps[0];

	int endsnp = intsnps[snpsize - 1];

	assert(grid>1);

	
	int i, left_boundary, right_boundary, left_min_boundary, right_min_boundary,
	    maxsizematrix = 0, dif;

	float step = (float)(endsnp - startsnp)/(grid-1), omega_position = (float)startsnp;

	if(step < 1)
	{
		fprintf(stdout, "\n\n WARNING: Gridsize is too large (%d) for the region between the first and last SNPs (%d)\n\n", grid, endsnp-startsnp+1);
		
		/* assert(endsnp-startsnp >= grid); */
		
	}

	int checksnp=0;

	for(i=0; i<grid; ++i)
	{
		omstruct[i].omegaRealPos = omega_position;

		left_boundary = omega_position - maxw;

		right_boundary = omega_position + maxw;

		left_min_boundary = omega_position - minw;

		right_min_boundary = omega_position + minw;

		omstruct[i].omegaPos = moveright(checksnp, intsnps, snpsize, omega_position);

		omstruct[i].leftIndex = moveleft(omstruct[i].omegaPos, intsnps, snpsize, left_boundary);

		omstruct[i].rightIndex = moveright(omstruct[i].omegaPos + 1, intsnps, snpsize, right_boundary);

		omstruct[i].leftminIndex = moveleft(omstruct[i].omegaPos, intsnps, snpsize, left_min_boundary);

		omstruct[i].rightminIndex = moveright(omstruct[i].omegaPos, intsnps, snpsize, right_min_boundary);

		omstruct[i].maxValue = 0.0;

		omstruct[i].maxLeftIndex = 0;
		omstruct[i].maxRightIndex = 0;

		/* change the minw if it is too small */
		while(omstruct[i].omegaPos - omstruct[i].leftminIndex + 1  < minsnps)
			--omstruct[i].leftminIndex;

		while( omstruct[i].rightminIndex - omstruct[i].omegaPos  < minsnps)
			++omstruct[i].rightminIndex;

		omstruct[i].valid = 1;

		/* do some checks for the minimum window boundaries */
		if(omstruct[i].leftminIndex < omstruct[i].leftIndex || omstruct[i].rightminIndex > omstruct[i].rightIndex)
			omstruct[i].valid = 0;

		/* do some checks for the maximum window boundaries */
		if(omstruct[i].omegaPos + 1 >= omstruct[i].rightIndex || omstruct[i].omegaPos <= omstruct[i].leftIndex)
			omstruct[i].valid = 0;

		dif = omstruct[i].rightIndex - omstruct[i].leftIndex + 1;

		if(omstruct[i].valid == 1 && dif  > MAXSIZE)
		{
			*changemaxw = 1;
			return;
		}

		if(dif > maxsizematrix)
			maxsizematrix = dif;

		omega_position += step;
	}

	*maxSizeMatrix = maxsizematrix;
}

int findOmegaBounds (alignment_struct * alignment, omega_struct * omega, int grid, int * maxw, int minw, int minsnps)
{
	int changemaxw = 1,
            matrixSizeMax=0,
            maxwc = *maxw;

	while (changemaxw)
	{
		changemaxw = 0;

		boundaries(omega, grid, alignment-> positionsInd, alignment->segsites, 	minw, minsnps, maxwc, &changemaxw, &matrixSizeMax);
	
		if(changemaxw)
			maxwc = (int)(DECREASE * maxwc);
	}

	*maxw = maxwc;

	return matrixSizeMax;	
}

int findNextValidOmega(omega_struct *omega, int lvw_i, int grid)
{
	int i=lvw_i+1;
	
	while(i<grid && omega[i].valid!=1)
		i++;
	
	return i;		
}

int validGridP(int cvw_i, int grid)
{
	if(cvw_i>=0 && cvw_i<grid)
		return 1;
	
	return 0;
}

float computeOmega (float LS, float RS, float TS, int k, int ksel2, int m, int msel2)
{
	float   numerator = (LS + RS) / (ksel2 + msel2);

	float denominator = (TS - LS - RS) / (k*m) + DENOMINATOR_OFFSET;

	float omega =  numerator / denominator;

	return omega;
}

void computeOmegaValues (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	float LS, RS, TS, tmpW = 0.0, maxW=0.0;

	int i, j, ksel2, msel2, k, m, maxLeftIndex=0, maxRightIndex=0,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex,

	rightMinIndexORIG = rightMinIndex,

	rightMaxIndexORIG = rightMaxIndex;

#ifdef _UNROLL
	int vw = 2;
	int iter, iterations = (rightMaxIndex-rightMinIndex+1) / vw;
	int finaliterations = (rightMaxIndex-rightMinIndex+1) % vw;
	int omegaSNIPIndexPlusOne = omegaSNIPIndex + 1;
	float maxW_0=0.0, maxW_1=0.0;
	int maxLeftIndex_0=0, maxRightIndex_0=0, maxLeftIndex_1=0, maxRightIndex_1=0;
#endif

#ifdef _USE_PTHREADS
#ifndef _USE_PTHREADS_MEMINT

	threadData_t * threadDataL = (threadData_t *) threadData;

	int t=0, 	

	tid = threadDataL->threadID,
	
	threads = threadDataL->threadTOTAL;
#endif
#endif	
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{

#ifdef _USE_PTHREADS

#ifndef _USE_PTHREADS_MEMINT
	  
	if(t%threads==tid)
		{
#endif
#endif	
			LS = correlationMatrix[omegaSNIPIndex][i];

			k = omegaSNIPIndex - i + 1;
		
			ksel2 = (k * (k-1)) / 2;

			if(borderTol > 0)
			{
				rightMinIndex = rightMinIndexORIG;

				rightMaxIndex = rightMaxIndexORIG;

			    //fprintf(stderr, "---------------------\nrightMinIndex: %d, rightMaxIndex: %d\n", rightMinIndex, rightMaxIndex);

				int leftSNPs = omegaSNIPIndex - i + 1;
				int equalRightPosition = omegaSNIPIndex + leftSNPs;
  
				rightMinIndex = max(rightMinIndex, equalRightPosition - borderTol);
				rightMaxIndex = min(rightMaxIndex, equalRightPosition + borderTol);
			    
			}
#ifdef _UNROLL
			j = rightMinIndex;

			for(iter=0;iter<iterations;iter++)
			{

				int j_0 = j;
				int j_1 = j_0 + 1;

				j=j+vw;

				float RS_0 = correlationMatrix[j_0][omegaSNIPIndexPlusOne];
				float RS_1 = correlationMatrix[j_1][omegaSNIPIndexPlusOne];

				int m_0 = j_0 - omegaSNIPIndex;
				int m_1 = j_1 - omegaSNIPIndex;
	
				int mmin1_0 = m_0 - 1;
				int mmin1_1 = m_1 - 1;

				int mmin1multm_0 = mmin1_0 * m_0;
				int mmin1multm_1 = mmin1_1 * m_1;

				int msel2_0 = mmin1multm_0 / 2;
				int msel2_1 = mmin1multm_1 / 2;
					
				float TS_0 = correlationMatrix[j_0][i];
				float TS_1 = correlationMatrix[j_1][i];

				int ksel2plusmsel2_0 = ksel2 + msel2_0;
				int ksel2plusmsel2_1 = ksel2 + msel2_1;

				float LSplusRS_0 = LS+RS_0;
				float LSplusRS_1 = LS+RS_1;


				float numerator_0 = LSplusRS_0 / ksel2plusmsel2_0;
				float numerator_1 = LSplusRS_1 / ksel2plusmsel2_1;

				float TSminusLS_0 = TS_0 - LS;
				float TSminusLS_1 = TS_1 - LS;

				float TSminusLSminusRS_0 = TSminusLS_0 - RS_0;
				float TSminusLSminusRS_1 = TSminusLS_1 - RS_1;

				int kmultm_0 = k * m_0;
				int kmultm_1 = k * m_1;

				float denominator_0 = TSminusLSminusRS_0 / kmultm_0;
				float denominator_1 = TSminusLSminusRS_1 / kmultm_1;

				float denominatorplusOFFSET_0 = denominator_0 + DENOMINATOR_OFFSET;
				float denominatorplusOFFSET_1 = denominator_1 + DENOMINATOR_OFFSET;

				float tmpW_0 =  numerator_0 / denominatorplusOFFSET_0;
				float tmpW_1 =  numerator_1 / denominatorplusOFFSET_1;
	
				if(tmpW_0>maxW_0)
				{
					maxW_0 = tmpW_0;
					maxLeftIndex_0 = i + omega[omegaIndex].leftIndex;
					maxRightIndex_0 = j_0 + omega[omegaIndex].leftIndex;
				}				
	
				if(tmpW_1>maxW_1)
				{
					maxW_1 = tmpW_1;
					maxLeftIndex_1 = i + omega[omegaIndex].leftIndex;
					maxRightIndex_1 = j_1 + omega[omegaIndex].leftIndex;
				}			

			}

			if(maxW_0>maxW_1)
			{
				maxW = maxW_0;
				maxLeftIndex = maxLeftIndex_0;
				maxRightIndex = maxRightIndex_0;
			}
			else
			{
				maxW = maxW_1;
				maxLeftIndex = maxLeftIndex_1;
				maxRightIndex = maxRightIndex_1;
			}


			if(finaliterations!=0)			
			{
			
				RS = correlationMatrix[j][omegaSNIPIndexPlusOne];
				m = j - omegaSNIPIndex;	
				int mmin1 = m - 1;
				int mmin1multm = mmin1 * m;
				msel2 = mmin1multm / 2;					
				TS = correlationMatrix[j][i];
				int ksel2plusmsel2 = ksel2 + msel2;
				float LSplusRS = LS+RS;
				float numerator = LSplusRS / ksel2plusmsel2;
				float TSminusLS = TS - LS;
				float TSminusLSminusRS = TSminusLS - RS;
				int kmultm = k * m;
				float denominator = TSminusLSminusRS / kmultm;
				float denominatorplusOFFSET = denominator + DENOMINATOR_OFFSET; 
				tmpW =  numerator / denominatorplusOFFSET;
	
				if(tmpW>maxW)
				{
					maxW = tmpW;
					maxLeftIndex = i + omega[omegaIndex].leftIndex;
					maxRightIndex = j + omega[omegaIndex].leftIndex;
				}

			}

			maxW_0 = maxW;
			maxLeftIndex_0 = maxLeftIndex;
			maxRightIndex_0 = maxRightIndex;
	
			maxW_1 = maxW;
			maxLeftIndex_1 = maxLeftIndex;
			maxRightIndex_1 = maxRightIndex;		

#else

		
			for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
			{
				RS = correlationMatrix[j][omegaSNIPIndex+1];

				m = j - omegaSNIPIndex;
	
				msel2 = (m * (m-1)) / 2;
					
				TS = correlationMatrix[j][i];

				tmpW = computeOmega(LS, RS, TS, k, ksel2, m, msel2);
	
				if(tmpW>maxW)
				{
					maxW = tmpW;
					maxLeftIndex = i + omega[omegaIndex].leftIndex;
					maxRightIndex = j + omega[omegaIndex].leftIndex;
				}
			}

#endif

#ifdef _USE_PTHREADS
#ifndef _USE_PTHREADS_MEMINT
		}
		t++;
#endif
#endif
	}
#ifdef _USE_PTHREADS
#ifndef _USE_PTHREADS_MEMINT
	threadDataL->threadArgCO->maxValue = maxW;
	threadDataL->threadArgCO->maxLeftIndex = maxLeftIndex;
	threadDataL->threadArgCO->maxRightIndex = maxRightIndex;	
#else
	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;
#endif
#else
	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;
#endif
}

void computeOmega_gpu(float * omegas, float * LSs, float * RSs, float * TSs, int * ks, int * ms, int outer_cnt, int inner_cnt, unsigned int total){
	// static cl_ulong p_start, p_end, p_total=0;
	int err=0;

	//set kernel arguments
	err = clSetKernelArg(omega_kernel, 6, sizeof(int), &inner_cnt);
	printCLErr(err,__LINE__,__FILE__);

	// size_t local = LOCAL_2;
	size_t global = total;// + (local-(total%local));

	// write values to GPU buffers
	// LS
	err=clEnqueueWriteBuffer(
			io_queue, LS_buffer, CL_FALSE, 0,
			outer_cnt*sizeof(float), LSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// RS
	err=clEnqueueWriteBuffer(
			io_queue, RS_buffer, CL_FALSE, 0,
			inner_cnt*sizeof(float), RSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// TS
	err=clEnqueueWriteBuffer(
			io_queue, TS_buffer, CL_FALSE, 0,
			total*sizeof(float), TSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);
	
	// k
	err=clEnqueueWriteBuffer(
			io_queue, k_buffer, CL_FALSE, 0,
			outer_cnt*sizeof(int), ks,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// m
	err=clEnqueueWriteBuffer(
			io_queue, m_buffer, CL_FALSE, 0,
			inner_cnt*sizeof(int), ms,
			0, NULL, &events[0]
			);
	printCLErr(err,__LINE__,__FILE__);

	//deploy kernel to execute program
	err=clEnqueueNDRangeKernel(
			io_queue, omega_kernel, 1, NULL, &global, 0,
			1, &events[0], &events[1]
			);
	printCLErr(err,__LINE__,__FILE__);
	// //deploy kernel to execute program

	// clWaitForEvents(1, &events[1]);

    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
    //                         &p_start, NULL);
	// printCLErr(err,__LINE__,__FILE__);
    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
    //                         &p_end, NULL);
	// printCLErr(err,__LINE__,__FILE__);

	// p_total += p_end - p_start;

	// printf("Compute: %lu\n",p_total);

	//read back omega values in omega buffer
	err=clEnqueueReadBuffer(
			io_queue, omega_buffer, CL_FALSE, 0,
			total*sizeof(float), omegas,
			1, &events[1], NULL
			);
	printCLErr(err,__LINE__,__FILE__);
}

void computeOmegaValues_gpu (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	// static double mtime0, mtime1, mtimetot = 0;
	float tmpW, maxW=0.0;
	static float * omegas = NULL, * LSs = NULL, * RSs = NULL, * TSs = NULL;

	unsigned int total;
	static int * ks = NULL, * ms = NULL;
	int i, j, maxLeftIndex=0, maxRightIndex=0, 
	
	outer_cnt=0, inner_cnt=0, outer_i=0, inner_i=0,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	outer_cnt = leftMinIndex-leftMaxIndex+1;
	inner_cnt = rightMaxIndex-rightMinIndex+1;
	total = outer_cnt * inner_cnt;
	
	omegas = malloc(sizeof(*omegas)*total);
	LSs = malloc(sizeof(*LSs)*outer_cnt);
	RSs = malloc(sizeof(*RSs)*inner_cnt);
	TSs = malloc(sizeof(*TSs)*total);
	ks = malloc(sizeof(*ks)*outer_cnt);
	ms = malloc(sizeof(*ms)*inner_cnt);

	if(omegas==NULL || LSs==NULL || RSs==NULL || TSs==NULL || ks==NULL || ms==NULL)
		printf("MALLOC error, %u, %u, %u\n", (leftMinIndex-leftMaxIndex+1), inner_cnt, total);
	
	// Doesn't make a difference in execution time
	// omegas = aligned_alloc (4096, GPU_BLOCK_MC*LOCAL_0*sizeof(float));
	// LSs = aligned_alloc (4096, GPU_BLOCK_KC*sizeof(float));
	// RSs = aligned_alloc (4096, GPU_BLOCK_KC*sizeof(float));
	// TSs = aligned_alloc (4096, GPU_BLOCK_MC*LOCAL_0*sizeof(float));
	// ks = aligned_alloc (4096, GPU_BLOCK_KC*sizeof(int));
	// ms = aligned_alloc (4096, GPU_BLOCK_KC*sizeof(int));
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{
		outer_i = leftMinIndex - i;

		LSs[outer_i] = correlationMatrix[omegaSNIPIndex][i];

		ks[outer_i] = omegaSNIPIndex - i + 1;

		if(borderTol > 0)	// Not implemented
		{
			int leftSNPs = omegaSNIPIndex - i + 1;
			printf("Left: %d\n",leftSNPs);
		}
		for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
		{
			inner_i = j - rightMinIndex;
			if(!outer_i)
			{
				RSs[inner_i] = correlationMatrix[j][omegaSNIPIndex+1];

				ms[inner_i] = j - omegaSNIPIndex;
			}
			
			TSs[(outer_i*inner_cnt)+inner_i] = correlationMatrix[j][i];
		}
	}
	// mtime0 = gettime();
	computeOmega_gpu(omegas, LSs, RSs, TSs, ks, ms, outer_cnt, inner_cnt, total);
	// mtime1 = gettime();
	// mtimetot += mtime1 - mtime0;

	// printf("GPU: %f\n",mtimetot);
	
	for(i=0;i<total;i++){
		tmpW = omegas[i];
		if(tmpW>maxW)
			{
				maxW = tmpW;
				maxLeftIndex = (leftMinIndex - (int)(i/inner_cnt)) + omega[omegaIndex].leftIndex;
				maxRightIndex = (rightMinIndex + (int)(i%inner_cnt)) + omega[omegaIndex].leftIndex;
			}
	}

	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;
	
	free(omegas);
	free(LSs);
	free(RSs);
	free(TSs);
	free(ks);
	free(ms);
}

void computeOmega_gpu2(float * omegas, float * LSs, float * RSs, float * TSs, int * ks, int * ms, unsigned int total){
	// static cl_ulong p_start, p_end, p_total=0;
	int err=0;

	// size_t local = LOCAL_2;
	size_t global = total;// + (local-(total%local));

	// write values to GPU buffers
	// LS
	err=clEnqueueWriteBuffer(
			io_queue, LS_buffer, CL_FALSE, 0,
			total*sizeof(float), LSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// RS
	err=clEnqueueWriteBuffer(
			io_queue, RS_buffer, CL_FALSE, 0,
			total*sizeof(float), RSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// TS
	err=clEnqueueWriteBuffer(
			io_queue, TS_buffer, CL_FALSE, 0,
			total*sizeof(float), TSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);
	
	// k
	err=clEnqueueWriteBuffer(
			io_queue, k_buffer, CL_FALSE, 0,
			total*sizeof(int), ks,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// m
	err=clEnqueueWriteBuffer(
			io_queue, m_buffer, CL_FALSE, 0,
			total*sizeof(int), ms,
			0, NULL, &events[0]
			);
	printCLErr(err,__LINE__,__FILE__);

	//deploy kernel to execute program
	err=clEnqueueNDRangeKernel(
			io_queue, omega_kernel, 1, NULL, &global, 0,
			1, &events[0], &events[1]
			);
	printCLErr(err,__LINE__,__FILE__);

	// clWaitForEvents(1, &events[1]);

    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
    //                         &p_start, NULL);
	// printCLErr(err,__LINE__,__FILE__);
    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
    //                         &p_end, NULL);
	// printCLErr(err,__LINE__,__FILE__);

	// p_total += p_end - p_start;

	// printf("Compute: %lu\n",p_total);

	//read back omega values in omega buffer
	err=clEnqueueReadBuffer(
			io_queue, omega_buffer, CL_TRUE, 0,
			total*sizeof(float), omegas,
			1, &events[1], NULL
			);
	printCLErr(err,__LINE__,__FILE__);
}

void computeOmegaValues_gpu2 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	// static double mtime0, mtime1, mtimetot = 0;
	float tmpW, maxW=0.0;
	static float * omegas = NULL, * LSs = NULL, * RSs = NULL, * TSs = NULL;

	unsigned int inner_cnt, total, index = 0;
	static int * ks = NULL, * ms = NULL;
	int i, j, maxLeftIndex=0, maxRightIndex=0,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	inner_cnt = rightMaxIndex-rightMinIndex+1;
	total = (leftMinIndex-leftMaxIndex+1) * inner_cnt;

	omegas = malloc(sizeof(*omegas)*total);
	LSs = malloc(sizeof(*LSs)*total);
	RSs = malloc(sizeof(*RSs)*total);
	TSs = malloc(sizeof(*TSs)*total);
	ks = malloc(sizeof(*ks)*total);
	ms = malloc(sizeof(*ms)*total);

	// int size = (total*sizeof(float))+(4096-((total*sizeof(float))%4096));

	// omegas = aligned_alloc (4096, size);
	// LSs = aligned_alloc (4096, size);
	// RSs = aligned_alloc (4096, size);
	// TSs = aligned_alloc (4096, size);
	// ks = aligned_alloc (4096, size);
	// ms = aligned_alloc (4096, size);

	// void * omegas = NULL, * LSs = NULL, * RSs = NULL, * TSs = NULL;
	// void * ks = NULL, * ms = NULL;
	// int pm;
	// pm=posix_memalign(&(omegas), 4096, size);
    // assert(!pm);
	// pm=posix_memalign(&(LSs), 4096, size);
    // assert(!pm);
	// pm=posix_memalign(&(RSs), 4096, size);
    // assert(!pm);
	// pm=posix_memalign(&(TSs), 4096, size);
    // assert(!pm);
	// pm=posix_memalign(&(ks), 4096, size);
    // assert(!pm);
	// pm=posix_memalign(&(ms), 4096, size);
    // assert(!pm);

	if(omegas==NULL || LSs==NULL || RSs==NULL || TSs==NULL || ks==NULL || ms==NULL)
	{
		printf("MALLOC error, %u, %u, %u\n", (leftMinIndex-leftMaxIndex+1), inner_cnt, total);
		exit(1);
	}
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{
		if(borderTol > 0)	// Not implemented
		{
			int leftSNPs = omegaSNIPIndex - i + 1;
			printf("Left: %d\n",leftSNPs);
		}
		for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
		{
			LSs[index] = correlationMatrix[omegaSNIPIndex][i];
			ks[index] = omegaSNIPIndex - i + 1;
			RSs[index] = correlationMatrix[j][omegaSNIPIndex+1];
			ms[index] = j - omegaSNIPIndex;
			TSs[index] = correlationMatrix[j][i];
			index++;
		}
	}

	// mtime0 = gettime();
	computeOmega_gpu2(omegas, LSs, RSs, TSs, ks, ms, total);
	// mtime1 = gettime();
	// mtimetot += mtime1 - mtime0;

	// printf("GPU: %f\n",mtimetot);

	for(i=0;i<total;i++){
		tmpW = omegas[i];
		if(tmpW>maxW)
		{
			maxW = tmpW;
			maxLeftIndex = (leftMinIndex - (int)(i/inner_cnt)) + omega[omegaIndex].leftIndex;
			maxRightIndex = (rightMinIndex + (int)(i%inner_cnt)) + omega[omegaIndex].leftIndex;
		}
	}

	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;
	
	free(omegas);
	free(LSs);
	free(RSs);
	free(TSs);
	free(ks);
	free(ms);
}

void computeOmega_gpu3(float * omegas, float * LR, int * km, float * TSs, int in_out_cnt, int inner_cnt, unsigned int total){
	static cl_ulong p_start, p_end, p_total=0;
	int err=0;

	// //set kernel arguments
	err = clSetKernelArg(omega_kernel, 4, sizeof(int), &inner_cnt);
	printCLErr(err,__LINE__,__FILE__);

	const size_t local = LOCAL_2;
	const size_t global = total + (LOCAL_2-(total%LOCAL_2));

	// write values to GPU buffers
	// LRkm
	err=clEnqueueWriteBuffer(
			io_queue, LR_buffer, CL_FALSE, 0,
			in_out_cnt*sizeof(float), LR,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// TS
	err=clEnqueueWriteBuffer(
			io_queue, TS_buffer, CL_FALSE, 0,
			total*sizeof(float), TSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// km
	err=clEnqueueWriteBuffer(
			io_queue, km_buffer, CL_FALSE, 0,
			in_out_cnt*sizeof(int), km,
			0, NULL, &events[0]
			);
	printCLErr(err,__LINE__,__FILE__);

	//deploy kernel to execute program
	err=clEnqueueNDRangeKernel(
			io_queue, omega_kernel, 1, NULL, &global, &local, //&local
			1, &events[0], &events[1]
			);
	printCLErr(err,__LINE__,__FILE__);

	clWaitForEvents(1, &events[1]);

    err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                            &p_start, NULL);
	printCLErr(err,__LINE__,__FILE__);
    err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
                            &p_end, NULL);
	printCLErr(err,__LINE__,__FILE__);

	p_total += p_end - p_start;

	// printf("Compute: %lu\n",p_total);
	printf("Compute: %lu\n",p_end - p_start);

	//read back omega values in omega buffer
	err=clEnqueueReadBuffer(
			io_queue, omega_buffer, CL_FALSE, 0,
			total*sizeof(float), omegas,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);
}

void computeOmegaValues_gpu3 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	// static double mtime0, mtime1, mtimetot = 0;
	float tmpW, maxW=0.0;
	static float * omegas = NULL, * LR = NULL, * TSs = NULL;

	unsigned int total;
	static int * km = NULL;
	int i, j, maxLeftIndex=0, maxRightIndex=0, 
	
	outer_cnt, inner_cnt, in_out_cnt, outer_i, inner_i, Lk_i, Rm_i,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	outer_cnt = leftMinIndex - leftMaxIndex + 1;
	inner_cnt = rightMaxIndex - rightMinIndex + 1;
	in_out_cnt = outer_cnt + inner_cnt;
	total = outer_cnt * inner_cnt;

	omegas = malloc(sizeof(*omegas)*total);
	LR = malloc(sizeof(*LR)*in_out_cnt);
	km = malloc(sizeof(*km)*in_out_cnt);
	TSs = malloc(sizeof(*TSs)*total);

	if(omegas==NULL || LR==NULL || km==NULL || TSs==NULL)
		printf("MALLOC error, %u, %u, %u\n", (leftMinIndex-leftMaxIndex+1), inner_cnt, total);
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{
		outer_i = leftMinIndex - i;
		Lk_i = inner_cnt + outer_i;

		LR[Lk_i] = correlationMatrix[omegaSNIPIndex][i];				// LSs

		km[Lk_i] = omegaSNIPIndex - i + 1;							// ks

		if(borderTol > 0)	// Not implemented
		{
			int leftSNPs = omegaSNIPIndex - i + 1;
			printf("Left: %d\n",leftSNPs);
		}
		for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
		{
			inner_i = j - rightMinIndex;
			Rm_i = inner_i;		// Compiler optimisation?!
			if(!outer_i)
			{
				LR[Rm_i] = correlationMatrix[j][omegaSNIPIndex+1];	// RSs

				km[Rm_i] = j - omegaSNIPIndex;						// ms
			}
			
			TSs[(outer_i*inner_cnt)+inner_i] = correlationMatrix[j][i];
		}
	}

	// mtime0 = gettime();
	computeOmega_gpu3(omegas, LR, km, TSs, in_out_cnt, inner_cnt, total);
	// mtime1 = gettime();
	// mtimetot += mtime1 - mtime0;

	// printf("GPU: %f\n",mtimetot);

	for(i=0;i<total;i++)
	{
		tmpW = omegas[i];
		// Validate
		// if(tmpW != omegasVal[i])
		// 	printf("NOT EQUAL!! %d, %f, %f\n",i,tmpW,omegasVal[i]);
		if(tmpW>maxW)
		{
			maxW = tmpW;
			maxLeftIndex = (leftMinIndex - (int)(i/inner_cnt)) + omega[omegaIndex].leftIndex;
			maxRightIndex = (rightMinIndex + (int)(i%inner_cnt)) + omega[omegaIndex].leftIndex;
		}
	}

	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;

	free(omegas);
	free(LR);
	free(km);
	free(TSs);
}

void computeOmega_gpu4(float * maxW, unsigned int * maxI, float * LRkm, float * TSs, int in_out_cnt, int inner_cnt, unsigned int total){
	// static cl_ulong p_start, p_end, p_total=0;

	int err = 0;
	const size_t local = group_size;
	const size_t global = work_items; 
	const cl_uint iter = total / global;

	// set kernel arguments
	err |= clSetKernelArg(omega_kernel, 6, sizeof(int), &inner_cnt);
	err |= clSetKernelArg(omega_kernel, 7, sizeof(cl_uint), &iter);
	printCLErr(err,__LINE__,__FILE__);

	// write values to GPU buffers
	// LRkm
	err=clEnqueueWriteBuffer(
			io_queue, LRkm_buffer, CL_FALSE, 0,
			in_out_cnt*sizeof(float), LRkm,
			0, NULL, &events[2]
			);
	printCLErr(err,__LINE__,__FILE__);

	// TS
	err=clEnqueueWriteBuffer(
			io_queue, TS_buffer, CL_FALSE, 0,
			total*sizeof(float), TSs,	//TSs must be of length: total + (global-(total%global))
			0, NULL, &events[0]
			);
	printCLErr(err,__LINE__,__FILE__);

	//deploy kernel to execute program
	err=clEnqueueNDRangeKernel(
			io_queue, omega_kernel, 1, NULL, &global, &local,
			1, &events[0], &events[1]
			);
	printCLErr(err,__LINE__,__FILE__);

	// clWaitForEvents(1, &events[1]);

    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
    //                         &p_start, NULL);
	// printCLErr(err,__LINE__,__FILE__);
    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
    //                         &p_end, NULL);
	// printCLErr(err,__LINE__,__FILE__);

	// p_total += p_end - p_start;

	// printf("%lu\n",p_total);

	err=clEnqueueReadBuffer(
			io_queue, omega_im, CL_FALSE, 0,
			sizeof(float), maxW,
			1, &events[1], NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	err=clEnqueueReadBuffer(
			io_queue, index_im, CL_FALSE, 0,
			sizeof(unsigned int), maxI,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);
}

void computeOmegaValues_gpu4 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	// static double mtime0, mtime1, mtimetot = 0;
	float maxW, * LRkm = NULL, * TSs = NULL;

	unsigned int maxI, total = 0;
	int i, j, maxLeftIndex=0, maxRightIndex=0, 
	
	outer_cnt, inner_cnt, in_out_cnt, outer_i, inner_i, Lk_i, Rm_i,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	outer_cnt = leftMinIndex - leftMaxIndex + 1;
	inner_cnt = rightMaxIndex - rightMinIndex + 1;
	// in_out_cnt = (outer_cnt + inner_cnt) * 2;
	total = outer_cnt * inner_cnt;

	const size_t global = work_items;
	unsigned int work_total = ((total + (global - 1)) / global) * global;	//Faster??

	unsigned int outer_work = ((work_total + (inner_cnt - 1)) / inner_cnt);

	in_out_cnt = (outer_work + inner_cnt) * 2;
	
	LRkm = malloc(sizeof(*LRkm)*in_out_cnt);
	TSs = malloc(sizeof(*TSs)*work_total);

	if(LRkm==NULL || TSs==NULL)
		printf("MALLOC error\n");

	// Check int to float cast inaccuracy LRkm[]
	if((abs(omegaSNIPIndex - leftMaxIndex + 1) > 16777216) || (abs(rightMaxIndex - omegaSNIPIndex) > 16777216)){
		printf("Error: inaccurate int to float cast");
		exit(1);
	}
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{
		outer_i = leftMinIndex - i;
		Lk_i = (outer_i * 2) + (inner_cnt * 2);

		LRkm[Lk_i] = correlationMatrix[omegaSNIPIndex][i];				// LSs

		LRkm[Lk_i+1] = omegaSNIPIndex - i + 1;							// ks

		if(borderTol > 0)	// Not implemented
		{
			int leftSNPs = omegaSNIPIndex - i + 1;
			printf("Left: %d\n",leftSNPs);
		}
		for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
		{
			inner_i = j - rightMinIndex;
			Rm_i = inner_i * 2;
			if(!outer_i)
			{
				LRkm[Rm_i] = correlationMatrix[j][omegaSNIPIndex+1];	// RSs

				LRkm[Rm_i+1] = j - omegaSNIPIndex;						// ms
			}
			
			TSs[(outer_i*inner_cnt)+inner_i] = correlationMatrix[j][i];
		}
	}

	for (i=total;i<work_total;i++){
		TSs[i] = 0;
	}

	for (i=(outer_cnt+inner_cnt)*2;i<in_out_cnt;i++){
		LRkm[i] = 1;
	}

	// mtime0 = gettime();
	computeOmega_gpu4(&maxW, &maxI, LRkm, TSs, in_out_cnt, inner_cnt, work_total);
	// mtime1 = gettime();
	// mtimetot += mtime1 - mtime0;

	// printf("GPU: %f\n",mtimetot);

	maxLeftIndex = (leftMinIndex - (int)(maxI/inner_cnt)) + omega[omegaIndex].leftIndex;
	maxRightIndex = (rightMinIndex + (int)(maxI%inner_cnt)) + omega[omegaIndex].leftIndex;

	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;
	
	free(LRkm);
	free(TSs);
}

void computeOmega_gpu5(float * maxW, unsigned int * maxI, float * LSs, float * RSs, float * TSs, int * ks, int * ms, unsigned int total){
	// static cl_ulong p_start, p_end, p_total=0;

	int err = 0;
	const size_t local = group_size;
	const size_t global = work_items;
	const cl_uint iter = total  / global;

	// set kernel arguments
	err |= clSetKernelArg(omega_kernel, 9, sizeof(cl_uint), &iter);
	printCLErr(err,__LINE__,__FILE__);

	// write values to GPU buffers
	err=clEnqueueWriteBuffer(
			io_queue, LS_buffer, CL_FALSE, 0,
			total*sizeof(float), LSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// RS
	err=clEnqueueWriteBuffer(
			io_queue, RS_buffer, CL_FALSE, 0,
			total*sizeof(float), RSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// TS
	err=clEnqueueWriteBuffer(
			io_queue, TS_buffer, CL_FALSE, 0,
			total*sizeof(float), TSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);
	
	// k
	err=clEnqueueWriteBuffer(
			io_queue, k_buffer, CL_FALSE, 0,
			total*sizeof(int), ks,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// m
	err=clEnqueueWriteBuffer(
			io_queue, m_buffer, CL_FALSE, 0,
			total*sizeof(int), ms,
			0, NULL, &events[0]
			);
	printCLErr(err,__LINE__,__FILE__);

	//deploy kernel to execute program
	err=clEnqueueNDRangeKernel(
			io_queue, omega_kernel, 1, NULL, &global, &local,
			1, &events[0], &events[1]
			);
	printCLErr(err,__LINE__,__FILE__);

	// clWaitForEvents(1, &events[1]);

    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
    //                         &p_start, NULL);
	// printCLErr(err,__LINE__,__FILE__);
    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
    //                         &p_end, NULL);
	// printCLErr(err,__LINE__,__FILE__);

	// p_total += p_end - p_start;

	// printf("Compute: %lu\n",p_total);

	err=clEnqueueReadBuffer(
			io_queue, omega_im, CL_FALSE, 0,
			sizeof(float), maxW,
			1, &events[1], NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	err=clEnqueueReadBuffer(
			io_queue, index_im, CL_FALSE, 0,
			sizeof(unsigned int), maxI,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);
}

void computeOmegaValues_gpu5 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	// static double mtime0, mtime1, mtimetot = 0;
	float maxW;
	static float * LSs = NULL, * RSs = NULL, * TSs = NULL;

	unsigned int inner_cnt, maxI, total, index = 0;
	static int * ks = NULL, * ms = NULL;
	int i, j, maxLeftIndex=0, maxRightIndex=0,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	inner_cnt = rightMaxIndex-rightMinIndex+1;
	total = (leftMinIndex-leftMaxIndex+1) * inner_cnt;

	const size_t global = work_items;
	unsigned int work_total = ((total + (global - 1)) / global) * global;	//Faster??

	LSs = malloc(sizeof(*LSs)*work_total);
	RSs = malloc(sizeof(*RSs)*work_total);
	TSs = malloc(sizeof(*TSs)*work_total);
	ks = malloc(sizeof(*ks)*work_total);
	ms = malloc(sizeof(*ms)*work_total);

	if(LSs==NULL || RSs==NULL || TSs==NULL || ks==NULL || ms==NULL)
	{
		printf("MALLOC error\n");
		exit(1);
	}
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{
		if(borderTol > 0)	// Not implemented
		{
			int leftSNPs = omegaSNIPIndex - i + 1;
			printf("Left: %d\n",leftSNPs);
		}
		for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
		{
			LSs[index] = correlationMatrix[omegaSNIPIndex][i];
			ks[index] = omegaSNIPIndex - i + 1;
			RSs[index] = correlationMatrix[j][omegaSNIPIndex+1];
			ms[index] = j - omegaSNIPIndex;
			TSs[index] = correlationMatrix[j][i];
			index++;
		}
	}

	for (i=total;i<work_total;i++){
		LSs[i] = 0;
		ks[i] = 0;
		RSs[i] = 0;
		ms[i] = 0;
		TSs[i] = 0;
	}

	// mtime0 = gettime();
	computeOmega_gpu5(&maxW, &maxI, LSs, RSs, TSs, ks, ms, work_total);
	// mtime1 = gettime();
	// mtimetot += mtime1 - mtime0;

	// printf("GPU: %f\n",mtimetot);

	maxLeftIndex = (leftMinIndex - (int)(maxI/inner_cnt)) + omega[omegaIndex].leftIndex;
	maxRightIndex = (rightMinIndex + (int)(maxI%inner_cnt)) + omega[omegaIndex].leftIndex;

	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;

	free(LSs);
	free(RSs);
	free(TSs);
	free(ks);
	free(ms);
}

void computeOmegaValues_gpu6 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	// static double mtime0, mtime1, mtimetot = 0;
	float tmpW, maxW=0.0;
	static float * omegas = NULL, * LR = NULL, * TSs = NULL;

	unsigned int total, * indexes = NULL, index = 0;
	static int * km = NULL;
	int i, j, maxLeftIndex=0, maxRightIndex=0, 
	
	outer_cnt, inner_cnt, in_out_cnt, outer_i, inner_i, Lk_i, Rm_i,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	outer_cnt = leftMinIndex - leftMaxIndex + 1;
	inner_cnt = rightMaxIndex - rightMinIndex + 1;
	total = outer_cnt * inner_cnt;

	const size_t global = work_items;
	unsigned int work_total = ((total + (global - 1)) / global) * global;	//Faster??

	unsigned int outer_work = ((work_total + (inner_cnt - 1)) / inner_cnt);

	in_out_cnt = outer_work + inner_cnt;

	omegas = malloc(sizeof(*omegas)*global);
	indexes = malloc(sizeof(*indexes)*global);
	LR = malloc(sizeof(*LR)*in_out_cnt);
	km = malloc(sizeof(*km)*in_out_cnt);
	TSs = malloc(sizeof(*TSs)*work_total);

	if(omegas==NULL || LR==NULL || km==NULL || TSs==NULL)
		printf("MALLOC error, %u, %u, %u\n", (leftMinIndex-leftMaxIndex+1), inner_cnt, total);
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{
		outer_i = leftMinIndex - i;
		Lk_i = inner_cnt + outer_i;

		LR[Lk_i] = correlationMatrix[omegaSNIPIndex][i];				// LSs

		km[Lk_i] = omegaSNIPIndex - i + 1;							// ks

		if(borderTol > 0)	// Not implemented
		{
			int leftSNPs = omegaSNIPIndex - i + 1;
			printf("Left: %d\n",leftSNPs);
		}
		for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
		{
			inner_i = j - rightMinIndex;
			Rm_i = inner_i;		// Compiler optimisation?!
			if(!outer_i)
			{
				LR[Rm_i] = correlationMatrix[j][omegaSNIPIndex+1];	// RSs

				km[Rm_i] = j - omegaSNIPIndex;						// ms
			}
			
			TSs[(outer_i*inner_cnt)+inner_i] = correlationMatrix[j][i];
		}
	}

	for (i=total;i<work_total;i++){
		TSs[i] = 0;
	}

	for (i=outer_cnt+inner_cnt;i<in_out_cnt;i++){
		LR[i] = 1;
		km[i] = 1;
	}

	// mtime0 = gettime();
	computeOmega_gpu7(omegas, indexes, LR, km, TSs, in_out_cnt, inner_cnt, work_total);
	// mtime1 = gettime();
	// mtimetot += mtime1 - mtime0;

	// printf("GPU: %f\n",mtimetot);

	for(i=0;i<global;i++)
	{
		tmpW = omegas[i];
		if(tmpW>maxW)
		{
			maxW = tmpW;
			index = indexes[i];
		}
	}
	
	maxLeftIndex = (leftMinIndex - (int)(index/inner_cnt)) + omega[omegaIndex].leftIndex;
	maxRightIndex = (rightMinIndex + (int)(index%inner_cnt)) + omega[omegaIndex].leftIndex;

	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;

	free(omegas);
	free(indexes);
	free(LR);
	free(km);
	free(TSs);
}

void computeOmega_gpu7(float * omegas, unsigned int * indexes, float * LR, int * km, float * TSs, int in_out_cnt, int inner_cnt, unsigned int total){
	static cl_ulong p_start, p_end, p_total=0;

	int err=0;
	const size_t local = group_size;
	const size_t global = work_items; 
	const cl_uint iter = total / global;

	// //set kernel arguments
	err |= clSetKernelArg(omega_kernel, 5, sizeof(int), &inner_cnt);
	err |= clSetKernelArg(omega_kernel, 6, sizeof(cl_uint), &iter);
	printCLErr(err,__LINE__,__FILE__);

	// write values to GPU buffers
	// LR
	err=clEnqueueWriteBuffer(
			io_queue, LR_buffer, CL_FALSE, 0,
			in_out_cnt*sizeof(float), LR,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// TS
	err=clEnqueueWriteBuffer(
			io_queue, TS_buffer, CL_FALSE, 0,
			total*sizeof(float), TSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// km
	err=clEnqueueWriteBuffer(
			io_queue, km_buffer, CL_FALSE, 0,
			in_out_cnt*sizeof(int), km,
			0, NULL, &events[0]
			);
	printCLErr(err,__LINE__,__FILE__);

	//deploy kernel to execute program
	err=clEnqueueNDRangeKernel(
			io_queue, omega_kernel, 1, NULL, &global, &local,
			1, &events[0], &events[1]
			);
	printCLErr(err,__LINE__,__FILE__);

	clWaitForEvents(1, &events[1]);

    err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                            &p_start, NULL);
	printCLErr(err,__LINE__,__FILE__);
    err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
                            &p_end, NULL);
	printCLErr(err,__LINE__,__FILE__);

	p_total += p_end - p_start;

	printf("%lu\n",p_total);

	//read back omega values in omega buffer
	err=clEnqueueReadBuffer(
			io_queue, omega_buffer, CL_FALSE, 0,
			global*sizeof(float), omegas,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	//read back indexes values in indexes buffer
	err=clEnqueueReadBuffer(
			io_queue, index_buffer, CL_FALSE, 0,
			global*sizeof(unsigned int), indexes,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);
}

void computeOmegaValues_gpu7 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	// static double mtime0, mtime1, mtimetot = 0;
	float tmpW, maxW=0.0;
	static float * omegas = NULL, * LR = NULL, * TSs = NULL;
	static unsigned int total_om = 0;
	unsigned int total, * indexes = NULL, index = 0;
	static int * km = NULL;
	int i, j, maxLeftIndex=0, maxRightIndex=0, 
	
	outer_cnt, inner_cnt, in_out_cnt, outer_i, inner_i, Lk_i, Rm_i,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	outer_cnt = leftMinIndex - leftMaxIndex + 1;
	inner_cnt = rightMaxIndex - rightMinIndex + 1;
	total = outer_cnt * inner_cnt;
	total_om += total;
	printf("total: %u\n",total_om);
	
	const size_t global = work_items;
	unsigned int work_total = ((total + (global - 1)) / global) * global;	//Faster??

	unsigned int outer_work = ((work_total + (inner_cnt - 1)) / inner_cnt);

	in_out_cnt = outer_work + inner_cnt;

	omegas = malloc(sizeof(*omegas)*global);
	indexes = malloc(sizeof(*indexes)*global);
	LR = malloc(sizeof(*LR)*in_out_cnt);
	km = malloc(sizeof(*km)*in_out_cnt);
	TSs = malloc(sizeof(*TSs)*work_total);

	if(omegas==NULL || LR==NULL || km==NULL || TSs==NULL)
		printf("MALLOC error\n");
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{
		outer_i = leftMinIndex - i;
		Lk_i = inner_cnt + outer_i;

		LR[Lk_i] = correlationMatrix[omegaSNIPIndex][i];				// LSs

		km[Lk_i] = omegaSNIPIndex - i + 1;							// ks

		if(borderTol > 0)	// Not implemented
		{
			int leftSNPs = omegaSNIPIndex - i + 1;
			printf("Left: %d\n",leftSNPs);
		}
		for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
		{
			inner_i = j - rightMinIndex;
			Rm_i = inner_i;		// Compiler optimisation?!
			if(!outer_i)
			{
				LR[Rm_i] = correlationMatrix[j][omegaSNIPIndex+1];	// RSs

				km[Rm_i] = j - omegaSNIPIndex;						// ms
			}
			
			TSs[(outer_i*inner_cnt)+inner_i] = correlationMatrix[j][i];
		}
	}

	for (i=total;i<work_total;i++){
		TSs[i] = 0;
	}

	for (i=outer_cnt+inner_cnt;i<in_out_cnt;i++){
		LR[i] = 1;
		km[i] = 1;
	}

	// mtime0 = gettime();
	computeOmega_gpu7(omegas, indexes, LR, km, TSs, in_out_cnt, inner_cnt, work_total);
	// mtime1 = gettime();
	// mtimetot += mtime1 - mtime0;

	// printf("GPU: %f\n",mtimetot);

	for(i=0;i<global;i++)
	{
		tmpW = omegas[i];
		if(tmpW>maxW)
		{
			maxW = tmpW;
			index = indexes[i];
		}
	}

	maxLeftIndex = (leftMinIndex - (int)(index/inner_cnt)) + omega[omegaIndex].leftIndex;
	maxRightIndex = (rightMinIndex + (int)(index%inner_cnt)) + omega[omegaIndex].leftIndex;

	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;

	free(omegas);
	free(indexes);
	free(LR);
	free(km);
	free(TSs);
}

void computeOmega_gpu8(float * omegas, unsigned int * indexes, float * LR, int * km, float * TSs, int in_out_cnt, int inner_cnt, unsigned int total, unsigned long int global){
	static cl_ulong p_start, p_end, p_total=0;

	int err=0;
	const size_t local = group_size;

	// //set kernel arguments
	err |= clSetKernelArg(omega_kernel, 5, sizeof(int), &inner_cnt);
	printCLErr(err,__LINE__,__FILE__);

	// write values to GPU buffers
	// LR
	err=clEnqueueWriteBuffer(
			io_queue, LR_buffer, CL_FALSE, 0,
			in_out_cnt*sizeof(float), LR,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// TS
	err=clEnqueueWriteBuffer(
			io_queue, TS_buffer, CL_FALSE, 0,
			total*sizeof(float), TSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// km
	err=clEnqueueWriteBuffer(
			io_queue, km_buffer, CL_FALSE, 0,
			in_out_cnt*sizeof(int), km,
			0, NULL, &events[0]
			);
	printCLErr(err,__LINE__,__FILE__);

	//deploy kernel to execute program
	err=clEnqueueNDRangeKernel(
			io_queue, omega_kernel, 1, NULL, &global, &local,
			1, &events[0], &events[1]
			);
	printCLErr(err,__LINE__,__FILE__);

	clWaitForEvents(1, &events[1]);

    err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                            &p_start, NULL);
	printCLErr(err,__LINE__,__FILE__);
    err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
                            &p_end, NULL);
	printCLErr(err,__LINE__,__FILE__);

	p_total += p_end - p_start;

	printf("%lu\n",p_total);

	//read back omega values in omega buffer
	err=clEnqueueReadBuffer(
			io_queue, omega_buffer, CL_FALSE, 0,
			global*sizeof(float), omegas,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	//read back indexes values in indexes buffer
	err=clEnqueueReadBuffer(
			io_queue, index_buffer, CL_FALSE, 0,
			global*sizeof(unsigned int), indexes,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);
}

void computeOmegaValues_gpu8 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	// static double mtime0, mtime1, mtimetot = 0;
	float tmpW, maxW=0.0;
	static float * omegas = NULL, * LR = NULL, * TSs = NULL;

	unsigned int work_total, * indexes = NULL, index = 0;
	static int * km = NULL;
	int i, j, maxLeftIndex=0, maxRightIndex=0,
	
	outer_work, outer_cnt, inner_cnt, in_out_cnt, outer_i, inner_i, Lk_i, Rm_i,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	outer_cnt = leftMinIndex - leftMaxIndex + 1;
	inner_cnt = rightMaxIndex - rightMinIndex + 1;

	const size_t local = group_size;
	outer_work = ((outer_cnt + (local - 1)) / local) * local;
	// work_total = outer_work * inner_cnt;
	work_total = outer_cnt * inner_cnt;
	// in_out_cnt = outer_work + inner_cnt;
	in_out_cnt = outer_cnt + inner_cnt;

	omegas = malloc(sizeof(*omegas)*outer_work);
	indexes = malloc(sizeof(*indexes)*outer_work);
	LR = malloc(sizeof(*LR)*in_out_cnt);
	km = malloc(sizeof(*km)*in_out_cnt);
	TSs = malloc(sizeof(*TSs)*work_total);

	if(omegas==NULL || LR==NULL || km==NULL || TSs==NULL)
		printf("MALLOC error\n");
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{
		outer_i = leftMinIndex - i;
		Lk_i = inner_cnt + outer_i;

		LR[Lk_i] = correlationMatrix[omegaSNIPIndex][i];				// LSs

		km[Lk_i] = omegaSNIPIndex - i + 1;							// ks

		if(borderTol > 0)	// Not implemented
		{
			int leftSNPs = omegaSNIPIndex - i + 1;
			printf("Left: %d\n",leftSNPs);
		}
		for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
		{
			inner_i = j - rightMinIndex;
			Rm_i = inner_i;		// Compiler optimisation?!
			if(!outer_i)		// Could be taken out of the for loop!
			{
				LR[Rm_i] = correlationMatrix[j][omegaSNIPIndex+1];	// RSs

				km[Rm_i] = j - omegaSNIPIndex;						// ms
			}
			
			TSs[(outer_i*inner_cnt)+inner_i] = correlationMatrix[j][i];
		}
	}

	// mtime0 = gettime();
	computeOmega_gpu8(omegas, indexes, LR, km, TSs, in_out_cnt, inner_cnt, work_total, outer_work);
	// mtime1 = gettime();
	// mtimetot += mtime1 - mtime0;

	// printf("GPU: %f\n",mtimetot);

	for(i=0;i<outer_cnt;i++)
	{
		tmpW = omegas[i];
		if(tmpW>maxW)
		{
			maxW = tmpW;
			index = indexes[i];
		}
	}

	maxLeftIndex = (leftMinIndex - (int)(index/inner_cnt)) + omega[omegaIndex].leftIndex;
	maxRightIndex = (rightMinIndex + (int)(index%inner_cnt)) + omega[omegaIndex].leftIndex;

	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;

	free(omegas);
	free(indexes);
	free(LR);
	free(km);
	free(TSs);
}

void computeOmega_gpu10(float * omegas, unsigned int * indexes, float * LR, int * km, float * TSs, int in_out_cnt, int inner_cnt, unsigned int total, unsigned long int global){
	static cl_ulong p_start, p_end, p_total=0;

	int err=0;
	const size_t local = group_size;

	// //set kernel arguments
	err |= clSetKernelArg(omega_kernel2, 5, sizeof(int), &inner_cnt);
	printCLErr(err,__LINE__,__FILE__);

	// write values to GPU buffers
	// LR
	err=clEnqueueWriteBuffer(
			io_queue, LR_buffer, CL_FALSE, 0,
			in_out_cnt*sizeof(float), LR,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// TS
	err=clEnqueueWriteBuffer(
			io_queue, TS_buffer, CL_FALSE, 0,
			total*sizeof(float), TSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// km
	err=clEnqueueWriteBuffer(
			io_queue, km_buffer, CL_FALSE, 0,
			in_out_cnt*sizeof(int), km,
			0, NULL, &events[0]
			);
	printCLErr(err,__LINE__,__FILE__);

	//deploy kernel to execute program
	err=clEnqueueNDRangeKernel(
			io_queue, omega_kernel2, 1, NULL, &global, &local,
			1, &events[0], &events[1]
			);
	printCLErr(err,__LINE__,__FILE__);

	clWaitForEvents(1, &events[1]);

    err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                            &p_start, NULL);
	printCLErr(err,__LINE__,__FILE__);
    err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
                            &p_end, NULL);
	printCLErr(err,__LINE__,__FILE__);

	p_total += p_end - p_start;

	// printf("Compute: %lu\n",p_total);
	printf("Compute: %lu\n",p_end - p_start);

	//read back omega values in omega buffer
	err=clEnqueueReadBuffer(
			io_queue, omega_buffer, CL_FALSE, 0,
			global*sizeof(float), omegas,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	//read back indexes values in indexes buffer
	err=clEnqueueReadBuffer(
			io_queue, index_buffer, CL_FALSE, 0,
			global*sizeof(unsigned int), indexes,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);
}

void computeOmegaValues_gpu10 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	// static double mtime0, mtime1, mtimetot = 0;
	float tmpW, maxW=0.0;
	static float * omegas = NULL, * LR = NULL, * TSs = NULL;

	unsigned int total, * indexes = NULL, index = 0;
	static int * km = NULL;
	int i, j, maxLeftIndex=0, maxRightIndex=0,
	
	outer_work, outer_cnt, inner_cnt, in_out_cnt, outer_i, inner_i, Lk_i, Rm_i,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	outer_cnt = leftMinIndex - leftMaxIndex + 1;
	inner_cnt = rightMaxIndex - rightMinIndex + 1;
	total = outer_cnt * inner_cnt;
	in_out_cnt = outer_cnt + inner_cnt;

	LR = malloc(sizeof(*LR)*in_out_cnt);
	km = malloc(sizeof(*km)*in_out_cnt);
	TSs = malloc(sizeof(*TSs)*total);

	if(LR==NULL || km==NULL || TSs==NULL)
		printf("MALLOC error\n");
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{
		outer_i = leftMinIndex - i;
		Lk_i = inner_cnt + outer_i;

		LR[Lk_i] = correlationMatrix[omegaSNIPIndex][i];				// LSs

		km[Lk_i] = omegaSNIPIndex - i + 1;							// ks

		if(borderTol > 0)	// Not implemented
		{
			int leftSNPs = omegaSNIPIndex - i + 1;
			printf("Left: %d\n",leftSNPs);
		}
		for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
		{
			inner_i = j - rightMinIndex;
			Rm_i = inner_i;		// Compiler optimisation?!
			if(!outer_i)
			{
				LR[Rm_i] = correlationMatrix[j][omegaSNIPIndex+1];	// RSs

				km[Rm_i] = j - omegaSNIPIndex;						// ms
			}
			
			TSs[(outer_i*inner_cnt)+inner_i] = correlationMatrix[j][i];
		}
	}

	if(total < 20000){	// or outer_cnt < 20
		// Do work items = total ("omega3")
		omegas = malloc(sizeof(*omegas)*total);

		// mtime0 = gettime();
		computeOmega_gpu3(omegas, LR, km, TSs, in_out_cnt, inner_cnt, total);
		// mtime1 = gettime();
		// mtimetot += mtime1 - mtime0;

		for(i=0;i<total;i++)
		{
			tmpW = omegas[i];
			if(tmpW>maxW)
			{
				maxW = tmpW;
				index = i;
			}
		}
	}
	else{
		// Do use for loop ("omega10") or ("omega9/8")
		const size_t local = group_size;
		outer_work = ((outer_cnt + (local - 1)) / local) * local;

		omegas = malloc(sizeof(*omegas)*outer_work);
		indexes = malloc(sizeof(*indexes)*outer_work);

		// mtime0 = gettime();
		computeOmega_gpu10(omegas, indexes, LR, km, TSs, in_out_cnt, inner_cnt, total, outer_work);
		// mtime1 = gettime();
		// mtimetot += mtime1 - mtime0;

		for(i=0;i<outer_cnt;i++)
		{
			tmpW = omegas[i];
			if(tmpW>maxW)
			{
				maxW = tmpW;
				index = indexes[i];
			}
		}
		free(indexes);
	}

	// printf("GPU: %f\n",mtimetot);

	maxLeftIndex = (leftMinIndex - (int)(index/inner_cnt)) + omega[omegaIndex].leftIndex;
	maxRightIndex = (rightMinIndex + (int)(index%inner_cnt)) + omega[omegaIndex].leftIndex;

	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;

	free(omegas);
	free(LR);
	free(km);
	free(TSs);
}

void computeOmega_gpu11(float * omegas, unsigned int * indexes, float * LR, int * km, float * TSs, int in_out_cnt, int inner_cnt, unsigned int total, unsigned long int global){
	cl_ulong p_start, p_end, p_total=0;

	int err=0;
	const size_t local = in_out_cnt;

	// printf("Global: %u, Local: %u, Groups: %u, Values: %u, Iterations: %u\n", global, local, global/local, total, inner_cnt);
	
	// //set kernel arguments
	err |= clSetKernelArg(omega_kernel, 7, sizeof(int), &inner_cnt);
	printCLErr(err,__LINE__,__FILE__);

	// write values to GPU buffers
	// TS
	err=clEnqueueWriteBuffer(
			io_queue, TS_buffer, CL_FALSE, 0,
			total*sizeof(float), TSs,
			0, NULL, &events[0]
			);
	printCLErr(err,__LINE__,__FILE__);
	
	//deploy kernel to execute program
	err=clEnqueueNDRangeKernel(
			io_queue, omega_kernel, 1, NULL, &global, &local,
			1, &events[0], &events[1]
			);
	printCLErr(err,__LINE__,__FILE__);

	clWaitForEvents(1, &events[1]);

    err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                            &p_start, NULL);
	printCLErr(err,__LINE__,__FILE__);
    err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
                            &p_end, NULL);
	printCLErr(err,__LINE__,__FILE__);

	p_total = p_end - p_start;

	// printf("%lu\n",p_total);

	*indexes = p_total;

	//read back omega values in omega buffer
	err=clEnqueueReadBuffer(
			io_queue, omega_buffer, CL_FALSE, 0,
			global*sizeof(float), omegas,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);
}

void computeOmegaValues_gpu11 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	static float * omegas = NULL, * LR = NULL, * TSs = NULL;

	unsigned int time;
	static int * km = NULL;
	int i;

	unsigned int optGroupSize 	= 0;
	unsigned int optGroups 		= 0;
	unsigned int tmpTime		= 0;
	unsigned int minTime		= UINT32_MAX;
	unsigned int values			= 1000000;
	unsigned int groups 		= 500;
	unsigned int group_size 	= 256;

	for(group_size = 128; group_size <= 256; group_size += 128){
		for(groups = 6; groups < 100; groups += 1){
			unsigned int global_size	= groups * group_size;
			unsigned int iterations		= (values + (global_size - 1)) / global_size;
			unsigned int work_values	= iterations * global_size;

			// global_size = (values + (group_size - 1)) / group_size * group_size;
			// work_values = global_size;

			omegas = malloc(sizeof(*omegas)*global_size);
			TSs = malloc(sizeof(*TSs)*work_values);

			for(i=0;i<work_values;i++){
				TSs[i] = 1;
			}

			for(i=0;i<100;i++){
				computeOmega_gpu11(omegas, &time, LR, km, TSs, group_size, iterations, work_values, global_size);
				tmpTime += time;
			}
			tmpTime /= 100;

			printf("%u;%u;%u;%u;%u;%u\n", global_size, group_size, groups, work_values, iterations, tmpTime);

			if(tmpTime < minTime){
				minTime = tmpTime;
				optGroupSize = group_size;
				optGroups = groups;
			}

			for(i=0;i<global_size;i++){
				if(omegas[i] != iterations)
					printf("FAULT ANSWER\n");
			}

			tmpTime = 0;

			free(omegas);
			free(TSs);
		}
		printf("Optimal group size: %u, Optimal groups: %u, Time: %u\n", optGroupSize, optGroups, minTime);
		minTime = UINT32_MAX;
	}

	free(LR);
	free(km);
}

void computeOmega_gpu13(float * omegas, unsigned int * indexes, float * LR, int * km, float * TSs, int in_out_cnt, int mult, int iter, int inner, unsigned int total){
	// static cl_ulong p_start, p_end, p_total=0;

	int err=0;
	const size_t local = group_size;
	const size_t global = work_items; 

	// //set kernel arguments
	err |= clSetKernelArg(omega_kernel, 5, sizeof(cl_int), &mult);
	err |= clSetKernelArg(omega_kernel, 6, sizeof(cl_int), &iter);
	err |= clSetKernelArg(omega_kernel, 7, sizeof(cl_int), &inner);
	printCLErr(err,__LINE__,__FILE__);

	// write values to GPU buffers
	// LR
	err=clEnqueueWriteBuffer(
			io_queue, LR_buffer, CL_FALSE, 0,
			in_out_cnt*sizeof(float), LR,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// TS
	err=clEnqueueWriteBuffer(
			io_queue, TS_buffer, CL_FALSE, 0,
			total*sizeof(float), TSs,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	// km
	err=clEnqueueWriteBuffer(
			io_queue, km_buffer, CL_FALSE, 0,
			in_out_cnt*sizeof(int), km,
			0, NULL, &events[0]
			);
	printCLErr(err,__LINE__,__FILE__);

	//deploy kernel to execute program
	err=clEnqueueNDRangeKernel(
			io_queue, omega_kernel, 1, NULL, &global, &local,
			1, &events[0], &events[1]
			);
	printCLErr(err,__LINE__,__FILE__);

	// clWaitForEvents(1, &events[1]);

    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
    //                         &p_start, NULL);
	// printCLErr(err,__LINE__,__FILE__);
    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
    //                         &p_end, NULL);
	// printCLErr(err,__LINE__,__FILE__);

	// p_total += p_end - p_start;

	// printf("%lu\n",p_total);

	//read back omega values in omega buffer
	err=clEnqueueReadBuffer(
			io_queue, omega_buffer, CL_FALSE, 0,
			global*sizeof(float), omegas,
			1, &events[1], NULL
			);
	printCLErr(err,__LINE__,__FILE__);

	//read back indexes values in indexes buffer
	err=clEnqueueReadBuffer(
			io_queue, index_buffer, CL_FALSE, 0,
			global*sizeof(unsigned int), indexes,
			0, NULL, NULL
			);
	printCLErr(err,__LINE__,__FILE__);
}

void computeOmegaValues_gpu13 (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	// static double mtime0, mtime1, mtimetot = 0;
	float tmpW, maxW=0.0;
	static float * omegas = NULL, * LR = NULL, * TSs = NULL;

	unsigned int * indexes = NULL, index = 0;
	static int * km = NULL;
	int i, j, maxLeftIndex=0, maxRightIndex=0, mult, iter, inner_work, outer_work, 
	
	work_total, outer_cnt, inner_cnt, in_out_cnt, inner_i=0, Lk_i=0, Rm_i, T_i=0,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	outer_cnt = leftMinIndex - leftMaxIndex + 1;
	inner_cnt = rightMaxIndex - rightMinIndex + 1;

	const size_t global = work_items;
	mult = global / outer_cnt;
	iter = (inner_cnt + (mult - 1)) / mult;
	inner_work = mult * iter;
	outer_work = global / mult;
	// work_total = inner_work * outer_work;
	work_total = iter * global;
	in_out_cnt = outer_work + inner_work;

	omegas = malloc(sizeof(*omegas)*global);
	indexes = malloc(sizeof(*indexes)*global);
	LR = malloc(sizeof(*LR)*in_out_cnt);
	km = malloc(sizeof(*km)*in_out_cnt);
	TSs = malloc(sizeof(*TSs)*work_total);

	if(omegas==NULL || LR==NULL || km==NULL || TSs==NULL)
		printf("MALLOC error\n");

	Rm_i = outer_work;

	// VALIDATE
	// float maxW2 = 0.0f;
	// int maxLeftIndex2;
	// int maxRightIndex2;

	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{
		LR[Lk_i] = correlationMatrix[omegaSNIPIndex][i];				// LSs

		km[Lk_i] = omegaSNIPIndex - i + 1;							// ks

		if(borderTol > 0)	// Not implemented
		{
			int leftSNPs = omegaSNIPIndex - i + 1;
			printf("Left: %d\n",leftSNPs);
		}
		for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
		{
			if(!Lk_i)
			{
				LR[Rm_i] = correlationMatrix[j][omegaSNIPIndex+1];	// RSs

				km[Rm_i] = j - omegaSNIPIndex;						// ms

				Rm_i++;
			}
			
			TSs[T_i] = correlationMatrix[j][i];

			inner_i++;

			if(inner_i % iter == 0)
				T_i = (Lk_i * mult) + (inner_i / iter);
			else
				T_i+=global;

			// VALIDATE
			// float LS = correlationMatrix[omegaSNIPIndex][i];
			// int k = omegaSNIPIndex - i + 1;
			// int ksel2 = (k * (k-1)) / 2;
			// float RS = correlationMatrix[j][omegaSNIPIndex+1];
			// int m = j - omegaSNIPIndex;
			// int msel2 = (m * (m-1)) / 2;
			// float TS = correlationMatrix[j][i];
			// tmpW = computeOmega(LS, RS, TS, k, ksel2, m, msel2);
			// if(tmpW>maxW2)
			// {
			// 	maxW2 = tmpW;
			// 	maxLeftIndex2 = i + omega[omegaIndex].leftIndex;
			// 	maxRightIndex2 = j + omega[omegaIndex].leftIndex;
			// }
		}
		inner_i = 0;
		Lk_i++;
		T_i = Lk_i * mult;
	}

	for(i=outer_work+inner_cnt;i<outer_work+inner_work;i++){
		LR[i] = 1;
		km[i] = 0;
	}

	// mtime0 = gettime();
	computeOmega_gpu13(omegas, indexes, LR, km, TSs, in_out_cnt, mult, iter, inner_cnt, work_total);
	// mtime1 = gettime();
	// mtimetot += mtime1 - mtime0;
	// printf("%f\n",mtimetot);

	for(i=0;i<outer_cnt*mult;i++)
	{
		tmpW = omegas[i];
		if(tmpW>maxW)
		{
			maxW = tmpW;
			index = indexes[i];
		}
	}

	maxLeftIndex = (leftMinIndex - (int)(index/inner_cnt)) + omega[omegaIndex].leftIndex;
	maxRightIndex = (rightMinIndex + (int)(index%inner_cnt)) + omega[omegaIndex].leftIndex;

	// VALIDATE
	// if(maxW2 != maxW)
	// 	printf("Max fault\n");

	// if(maxLeftIndex2 != maxLeftIndex)
	// 	printf("Left fault\n");

	// if(maxRightIndex2 != maxRightIndex)
	// 	printf("Right fault\n");

	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;

	free(omegas);
	free(indexes);
	free(LR);
	free(km);
	free(TSs);
}

#ifdef _USE_PTHREADS
#ifdef _USE_PTHREADS_MEMINT
void computeOmegas (alignment_struct * alignment, omega_struct * omega, int omegaIndex, void * threadData, cor_t ** correlationMatrix)
{
	computeOmegaValues (omega, omegaIndex, correlationMatrix, NULL);
}
#ifdef _USE_PTHREADS_MULTI
void initmaxOmegaValueThreadsMULTI(int size)
{
	int j;

	for(j=0;j<size;j++)
		maxOmegaValueThreadsMULTI[j]=-1.0;
}

void getmaxOmegaValueThreadsMULTI(omega_struct * omega, int cvw_i, int size)
{

	int j;

	omega[cvw_i].maxValue = maxOmegaValueThreadsMULTI[0];
	omega[cvw_i].maxLeftIndex  = maxOmegaLeftIndexThreadsMULTI[0];
	omega[cvw_i].maxRightIndex = maxOmegaRightIndexThreadsMULTI[0];

	for(j=1;j<size;j++)
		if(maxOmegaValueThreadsMULTI[j]>omega[cvw_i].maxValue)
		{
			omega[cvw_i].maxValue = maxOmegaValueThreadsMULTI[j];
			omega[cvw_i].maxLeftIndex  = maxOmegaLeftIndexThreadsMULTI[j];
			omega[cvw_i].maxRightIndex = maxOmegaRightIndexThreadsMULTI[j];
		}
}

void computeOmegaValuesMULTI (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, int ttid, int tthreads, int tid)
{
	float LS, RS, TS, tmpW = 0, maxW=0;

	int i, j, ksel2, msel2, k, m, maxLeftIndex=0, maxRightIndex=0,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex,

	rightMinIndexORIG = rightMinIndex,

	rightMaxIndexORIG = rightMaxIndex;

	int t =0;

#ifdef _UNROLL
	int vw = 2;
	int iter, iterations = (rightMaxIndex-rightMinIndex+1) / vw;
	int finaliterations = (rightMaxIndex-rightMinIndex+1) % vw;
	int omegaSNIPIndexPlusOne = omegaSNIPIndex + 1;
	float maxW_0=0.0, maxW_1=0.0;
	int maxLeftIndex_0=0, maxRightIndex_0=0, maxLeftIndex_1=0, maxRightIndex_1=0;
#endif
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{

		if(t%tthreads==ttid)
		{
			LS = correlationMatrix[omegaSNIPIndex][i];

			k = omegaSNIPIndex - i + 1;
		
			ksel2 = (k * (k-1)) / 2;

			if(borderTol > 0)
			{
				rightMinIndex = rightMinIndexORIG;

				rightMaxIndex = rightMaxIndexORIG;

			    //fprintf(stderr, "---------------------\nrightMinIndex: %d, rightMaxIndex: %d\n", rightMinIndex, rightMaxIndex);

				int leftSNPs = omegaSNIPIndex - i + 1;
				int equalRightPosition = omegaSNIPIndex + leftSNPs;
  
				rightMinIndex = max(rightMinIndex, equalRightPosition - borderTol);
				rightMaxIndex = min(rightMaxIndex, equalRightPosition + borderTol);
			    
			}
		
			#ifdef _UNROLL
			j = rightMinIndex;

			for(iter=0;iter<iterations;iter++)
			{

				int j_0 = j;
				int j_1 = j_0 + 1;

				j=j+vw;

				float RS_0 = correlationMatrix[j_0][omegaSNIPIndexPlusOne];
				float RS_1 = correlationMatrix[j_1][omegaSNIPIndexPlusOne];

				int m_0 = j_0 - omegaSNIPIndex;
				int m_1 = j_1 - omegaSNIPIndex;
	
				int mmin1_0 = m_0 - 1;
				int mmin1_1 = m_1 - 1;

				int mmin1multm_0 = mmin1_0 * m_0;
				int mmin1multm_1 = mmin1_1 * m_1;

				int msel2_0 = mmin1multm_0 / 2;
				int msel2_1 = mmin1multm_1 / 2;
					
				float TS_0 = correlationMatrix[j_0][i];
				float TS_1 = correlationMatrix[j_1][i];

				int ksel2plusmsel2_0 = ksel2 + msel2_0;
				int ksel2plusmsel2_1 = ksel2 + msel2_1;

				float LSplusRS_0 = LS+RS_0;
				float LSplusRS_1 = LS+RS_1;


				float numerator_0 = LSplusRS_0 / ksel2plusmsel2_0;
				float numerator_1 = LSplusRS_1 / ksel2plusmsel2_1;

				float TSminusLS_0 = TS_0 - LS;
				float TSminusLS_1 = TS_1 - LS;

				float TSminusLSminusRS_0 = TSminusLS_0 - RS_0;
				float TSminusLSminusRS_1 = TSminusLS_1 - RS_1;

				int kmultm_0 = k * m_0;
				int kmultm_1 = k * m_1;

				float denominator_0 = TSminusLSminusRS_0 / kmultm_0;
				float denominator_1 = TSminusLSminusRS_1 / kmultm_1;

				float denominatorplusOFFSET_0 = denominator_0 + DENOMINATOR_OFFSET;
				float denominatorplusOFFSET_1 = denominator_1 + DENOMINATOR_OFFSET;

				float tmpW_0 =  numerator_0 / denominatorplusOFFSET_0;
				float tmpW_1 =  numerator_1 / denominatorplusOFFSET_1;
	
				if(tmpW_0>maxW_0)
				{
					maxW_0 = tmpW_0;
					maxLeftIndex_0 = i + omega[omegaIndex].leftIndex;
					maxRightIndex_0 = j_0 + omega[omegaIndex].leftIndex;
				}				
	
				if(tmpW_1>maxW_1)
				{
					maxW_1 = tmpW_1;
					maxLeftIndex_1 = i + omega[omegaIndex].leftIndex;
					maxRightIndex_1 = j_1 + omega[omegaIndex].leftIndex;
				}			

			}

			if(maxW_0>maxW_1)
			{
				maxW = maxW_0;
				maxLeftIndex = maxLeftIndex_0;
				maxRightIndex = maxRightIndex_0;
			}
			else
			{
				maxW = maxW_1;
				maxLeftIndex = maxLeftIndex_1;
				maxRightIndex = maxRightIndex_1;
			}


			if(finaliterations!=0)			
			{
			
				RS = correlationMatrix[j][omegaSNIPIndexPlusOne];
				m = j - omegaSNIPIndex;	
				int mmin1 = m - 1;
				int mmin1multm = mmin1 * m;
				msel2 = mmin1multm / 2;					
				TS = correlationMatrix[j][i];
				int ksel2plusmsel2 = ksel2 + msel2;
				float LSplusRS = LS+RS;
				float numerator = LSplusRS / ksel2plusmsel2;
				float TSminusLS = TS - LS;
				float TSminusLSminusRS = TSminusLS - RS;
				int kmultm = k * m;
				float denominator = TSminusLSminusRS / kmultm;
				float denominatorplusOFFSET = denominator + DENOMINATOR_OFFSET; 
				tmpW =  numerator / denominatorplusOFFSET;
	
				if(tmpW>maxW)
				{
					maxW = tmpW;
					maxLeftIndex = i + omega[omegaIndex].leftIndex;
					maxRightIndex = j + omega[omegaIndex].leftIndex;
				}

			}

			maxW_0 = maxW;
			maxLeftIndex_0 = maxLeftIndex;
			maxRightIndex_0 = maxRightIndex;
	
			maxW_1 = maxW;
			maxLeftIndex_1 = maxLeftIndex;
			maxRightIndex_1 = maxRightIndex;		

#else

		
			for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
			{
				RS = correlationMatrix[j][omegaSNIPIndex+1];

				m = j - omegaSNIPIndex;
	
				msel2 = (m * (m-1)) / 2;
					
				TS = correlationMatrix[j][i];

				tmpW = computeOmega(LS, RS, TS, k, ksel2, m, msel2);
	
				if(tmpW>maxW)
				{
					maxW = tmpW;
					maxLeftIndex = i + omega[omegaIndex].leftIndex;
					maxRightIndex = j + omega[omegaIndex].leftIndex;
				}
			}

#endif

		}
		t++;
	}

	maxOmegaValueThreadsMULTI[tid]=maxW;
	maxOmegaLeftIndexThreadsMULTI[tid]=maxLeftIndex;
	maxOmegaRightIndexThreadsMULTI[tid]=maxRightIndex;
}

void computeOmegasMULTI (alignment_struct * alignment, omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, int ttid, int tthreads, int tid)
{
	computeOmegaValuesMULTI (omega, omegaIndex, correlationMatrix, ttid, tthreads, tid);
}

#endif
#else
void computeOmegasThread (alignment_struct * alignment, omega_struct * omega, int omegaIndex, threadData_t * threadData)
{	
	computeOmegaValues (omega, omegaIndex, alignment->correlationMatrix, threadData);
}

void omegasThread(threadData_t * currentThread)
{
	alignment_struct * alignment = currentThread->threadArgCO->alignment;

	omega_struct * omega = currentThread->threadArgCO->omega;

	int omegaIndex = currentThread->threadArgCO->omegaIndex;

        computeOmegasThread (alignment,omega,omegaIndex,currentThread);	
}

void setThreadArgumentsCO(threadData_t * threadData, int tid, alignment_struct * alignment, omega_struct * omega, int omegaIndex)
{
	threadData[tid].threadArgCO->alignment=alignment;
	threadData[tid].threadArgCO->omega=omega;
	threadData[tid].threadArgCO->omegaIndex=omegaIndex;
}

void getAllThreadMaxOmega(threadData_t * threadData, omega_struct * omega, int omegaIndex)
{
	int i, threads = threadData[0].threadTOTAL;

	omega[omegaIndex].maxValue = threadData[0].threadArgCO->maxValue;
	omega[omegaIndex].maxLeftIndex = threadData[0].threadArgCO->maxLeftIndex;
	omega[omegaIndex].maxRightIndex = threadData[0].threadArgCO->maxRightIndex;

	for(i=1;i<threads;i++)
	{
		if(threadData[i].threadArgCO->maxValue>omega[omegaIndex].maxValue)
		{
			omega[omegaIndex].maxValue = threadData[i].threadArgCO->maxValue;
			omega[omegaIndex].maxLeftIndex = threadData[i].threadArgCO->maxLeftIndex;
			omega[omegaIndex].maxRightIndex = threadData[i].threadArgCO->maxRightIndex;
		}
	}
}

void computeOmegaValues_THREADS (alignment_struct * alignment, omega_struct * omega, int omegaIndex, threadData_t * threadData)
{
	int i, threads = threadData[0].threadTOTAL;

	for(i=0;i<threads;i++)
		setThreadArgumentsCO(threadData, i, alignment, omega, omegaIndex);

	startThreadOperations(threadData, COMPUTEOMEGAS);

	getAllThreadMaxOmega(threadData,omega, omegaIndex);
}

void computeOmegas (alignment_struct * alignment, omega_struct * omega, int omegaIndex, void * threadData, cor_t ** correlationMatrix)
{

	threadData_t * threadDataL = (threadData_t *) threadData;

	computeOmegaValues_THREADS (alignment, omega, omegaIndex, threadDataL);
}
#endif
#else
void computeOmegas (alignment_struct * alignment, omega_struct * omega, int omegaIndex, void * threadData, cor_t ** correlationMatrix)
{
	computeOmegaValues (omega, omegaIndex, alignment->correlationMatrix, NULL);
}

void computeOmegas_gpu (alignment_struct * alignment, omega_struct * omega, int omegaIndex, void * threadData, cor_t ** correlationMatrix)
{
	computeOmegaValues_gpu13 (omega, omegaIndex, alignment->correlationMatrix, NULL);
}
#endif

void appendOmegaResultToFile (alignment_struct * alignment, omega_struct * omega, int omegaIndex, int gridIndex, FILE * fpOut, int resultType)
{
	if (resultType==RESULTS_ALL)
		if(omega[omegaIndex].valid)
			fprintf(fpOut,"%.4f\t%f\t%d\t%d\t%d\n", omega[omegaIndex].omegaRealPos, omega[omegaIndex].maxValue, alignment->positionsInd[omega[omegaIndex].maxLeftIndex], 
							      alignment->positionsInd[omega[omegaIndex].maxRightIndex], omega[omegaIndex].valid );
		else
			fprintf(fpOut,"%.4f\t%f\t%d\t%d\t%d\n", omega[omegaIndex].omegaRealPos, 0.0, 0, 0, 0);
	else
		if(omega[omegaIndex].valid)
			fprintf(fpOut,"%.4f\t%f\n", omega[omegaIndex].omegaRealPos, omega[omegaIndex].maxValue);
		else
			fprintf(fpOut,"%.4f\t%f\n", omega[omegaIndex].omegaRealPos, 0.0);
}

void maxOmegaResultReport (float * maxomegaRealPos, float * maxomegamaxValue, int * maxomegamaxLeftIndex, int * maxomegamaxRightIndex, int grid, alignment_struct * alignment, omega_struct * omega, FILE * fpInfo)
{
	*maxomegaRealPos = 0.0;
	*maxomegamaxValue = 0.0;
	*maxomegamaxLeftIndex = -1;
	*maxomegamaxRightIndex = -1;

	int i;

	for(i=0;i<grid;i++)
	{
		if(omega[i].maxValue>*maxomegamaxValue)
		{
			*maxomegamaxValue = omega[i].maxValue;
			*maxomegaRealPos = omega[i].omegaRealPos;
			*maxomegamaxLeftIndex = alignment->positionsInd[omega[i].maxLeftIndex];
			*maxomegamaxRightIndex = alignment->positionsInd[omega[i].maxRightIndex];		
		}
	}

	fprintf(stdout, "\t\tMax Omega:\t\t%f\n",*maxomegamaxValue);
	fprintf(stdout, "\t\tLocation:\t\t%.2f\n",*maxomegaRealPos);
	fprintf(stdout, "\t\tLeftmost SNP:\t\t%d\n",*maxomegamaxLeftIndex);
	fprintf(stdout, "\t\tRightmost SNP:\t\t%d\n\n\n",*maxomegamaxRightIndex);

	fprintf(fpInfo, "\t\tMax Omega:\t\t%f\n",*maxomegamaxValue);
	fprintf(fpInfo, "\t\tLocation:\t\t%.2f\n",*maxomegaRealPos);
	fprintf(fpInfo, "\t\tLeftmost SNP:\t\t%d\n",*maxomegamaxLeftIndex);
	fprintf(fpInfo, "\t\tRightmost SNP:\t\t%d\n\n\n",*maxomegamaxRightIndex);
	
}

// GPU OpenCL stuff sequential
static void create_program_with_source(cl_program *program,
                                       cl_context *context,
                                       const char *program_file)
{
    // read in program source file and create `program`
    FILE *program_handle;
    size_t program_size;
    char *program_buffer;
    int err;

    program_handle=fopen(program_file, "r");
    //assert(program_handle);
    fseek(program_handle, 0, SEEK_END);
    program_size=ftell(program_handle);
    rewind(program_handle);

    program_buffer=(char*) malloc(program_size+1);
    assert(program_buffer);
    program_buffer[program_size]='\0';

    size_t ret_code=fread(program_buffer, sizeof(char), program_size, program_handle);
    if(ret_code != program_size)
    {
        printf("error reading file %s\n", program_file);
        if(feof(program_handle))
        {
            printf("Error reading file `%s`: unexpected end of file\n",program_file);
        }
        else if(ferror(program_handle))
        {
            perror("Error reading file `%s`");
        }
        exit(1);
    }
    fclose(program_handle);

    *program=clCreateProgramWithSource(*context, 1, (const char**) &program_buffer,
                                         &program_size, &err);

    free(program_buffer);

    printCLErr(err,__LINE__,__FILE__);
}

void printCLErr(cl_int err,int line, char* file)
{
    if(err == CL_SUCCESS) {
        return;
    }
    switch (err)
    {
        case CL_SUCCESS:
            printf("CL_SUCCESS\n");
            break;
        case CL_DEVICE_NOT_FOUND:
            printf("CL_DEVICE_NOT_FOUND\n");
            break;
        case CL_DEVICE_NOT_AVAILABLE:
            printf("CL_DEVICE_NOT_AVAILABLE\n");
            break;
        case CL_COMPILER_NOT_AVAILABLE:
            printf("CL_COMPILER_NOT_AVAILABLE\n");
            break;
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
            printf("CL_MEM_OBJECT_ALLOCATION_FAILURE\n");
            break;
        case CL_OUT_OF_RESOURCES:
            printf("CL_OUT_OF_RESOURCES\n");
            break;
        case CL_OUT_OF_HOST_MEMORY:
            printf("CL_OUT_OF_HOST_MEMORY\n");
            break;
        case CL_PROFILING_INFO_NOT_AVAILABLE:
            printf("CL_PROFILING_INFO_NOT_AVAILABLE\n");
            break;
        case CL_MEM_COPY_OVERLAP:
            printf("CL_MEM_COPY_OVERLAP\n");
            break;
        case CL_IMAGE_FORMAT_MISMATCH:
            printf("CL_IMAGE_FORMAT_MISMATCH\n");
            break;
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
            printf("CL_IMAGE_FORMAT_NOT_SUPPORTED\n");
            break;
        case CL_BUILD_PROGRAM_FAILURE:
            printf("CL_BUILD_PROGRAM_FAILURE\n");
            break;
        case CL_MAP_FAILURE:
            printf("CL_MAP_FAILURE\n");
            break;
        case CL_INVALID_VALUE:
            printf("CL_INVALID_VALUE\n");
            break;
        case CL_INVALID_DEVICE_TYPE:
            printf("CL_INVALID_DEVICE_TYPE\n");
            break;
        case CL_INVALID_PLATFORM:
            printf("CL_INVALID_PLATFORM\n");
            break;
        case CL_INVALID_DEVICE:
            printf("CL_INVALID_DEVICE\n");
            break;
        case CL_INVALID_CONTEXT:
            printf("CL_INVALID_CONTEXT\n");
            break;
        case CL_INVALID_QUEUE_PROPERTIES:
            printf("CL_INVALID_QUEUE_PROPERTIES\n");
            break;
        case CL_INVALID_COMMAND_QUEUE:
            printf("CL_INVALID_COMMAND_QUEUE\n");
            break;
        case CL_INVALID_HOST_PTR:
            printf("CL_INVALID_HOST_PTR\n");
            break;
        case CL_INVALID_MEM_OBJECT:
            printf("CL_INVALID_MEM_OBJECT\n");
            break;
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
            printf("CL_INVALID_IMAGE_FORMAT_DESCRIPTOR\n");
            break;
        case CL_INVALID_IMAGE_SIZE:
            printf("CL_INVALID_IMAGE_SIZE\n");
            break;
        case CL_INVALID_SAMPLER:
            printf("CL_INVALID_SAMPLER\n");
            break;
        case CL_INVALID_BINARY:
            printf("CL_INVALID_BINARY\n");
            break;
        case CL_INVALID_BUILD_OPTIONS:
            printf("CL_INVALID_BUILD_OPTIONS\n");
            break;
        case CL_INVALID_PROGRAM:
            printf("CL_INVALID_PROGRAM\n");
            break;
        case CL_INVALID_PROGRAM_EXECUTABLE:
            printf("CL_INVALID_PROGRAM_EXECUTABLE\n");
            break;
        case CL_INVALID_KERNEL_NAME:
            printf("CL_INVALID_KERNEL_NAME\n");
            break;
        case CL_INVALID_KERNEL_DEFINITION:
            printf("CL_INVALID_KERNEL_DEFINITION\n");
            break;
        case CL_INVALID_KERNEL:
            printf("CL_INVALID_KERNEL\n");
            break;
        case CL_INVALID_ARG_INDEX:
            printf("CL_INVALID_ARG_INDEX\n");
            break;
        case CL_INVALID_ARG_VALUE:
            printf("CL_INVALID_ARG_VALUE\n");
            break;
        case CL_INVALID_ARG_SIZE:
            printf("CL_INVALID_ARG_SIZE\n");
            break;
        case CL_INVALID_KERNEL_ARGS:
            printf("CL_INVALID_KERNEL_ARGS\n");
            break;
        case CL_INVALID_WORK_DIMENSION:
            printf("CL_INVALID_WORK_DIMENSION\n");
            break;
        case CL_INVALID_WORK_GROUP_SIZE:
            printf("CL_INVALID_WORK_GROUP_SIZE\n");
            break;
        case CL_INVALID_WORK_ITEM_SIZE:
            printf("CL_INVALID_WORK_ITEM_SIZE\n");
            break;
        case CL_INVALID_GLOBAL_OFFSET:
            printf("CL_INVALID_GLOBAL_OFFSET\n");
            break;
        case CL_INVALID_EVENT_WAIT_LIST:
            printf("CL_INVALID_EVENT_WAIT_LIST\n");
            break;
        case CL_INVALID_EVENT:
            printf("CL_INVALID_EVENT\n");
            break;
        case CL_INVALID_OPERATION:
            printf("CL_INVALID_OPERATION\n");
            break;
        case CL_INVALID_GL_OBJECT:
            printf("CL_INVALID_GL_OBJECT\n");
            break;
        case CL_INVALID_BUFFER_SIZE:
            printf("CL_INVALID_BUFFER_SIZE\n");
            break;
        case CL_INVALID_MIP_LEVEL:
            printf("CL_INVALID_MIP_LEVEL\n");
            break;
        case CL_INVALID_GLOBAL_WORK_SIZE:
            printf("CL_INVALID_GLOBAL_WORK_SIZE\n");
            break;
        default:
            printf("don't know what error that was\n");
            break;
    }
    printf ("Blah error at %s (%d)\n", file, line);
}

void gpu_init(void)
{
    // ---- OpenCL stuff ---------------
    int err;

    // query number of platforms we have
    unsigned int num_platforms=0;
    err=clGetPlatformIDs(1, NULL, &num_platforms);
    printCLErr(err,__LINE__,__FILE__);

    // allocate array, an entry for each platform found
    platforms=(cl_platform_id *) malloc(sizeof(cl_platform_id) * num_platforms);
    assert(platforms);
    // place the `cl_platform_id` structures in the platforms array
    err=clGetPlatformIDs(num_platforms, platforms, NULL);
    printCLErr(err,__LINE__,__FILE__);

    // determine number of devices
    // NOTE: arbitrarily pick first platform
    unsigned int num_devices=0;
    err=clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 1, NULL, &num_devices);
    printCLErr(err,__LINE__,__FILE__);
    // allocate memory for device array
    devices=(cl_device_id*) malloc(sizeof(cl_device_id) * num_devices);
    assert(devices);
    // populate device array
    // NOTE: arbitrarily pick first platform
    err=clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, num_devices, devices, NULL);
    printCLErr(err,__LINE__,__FILE__);
   
    context=clCreateContext(NULL, 1, &devices[1], NULL, NULL, &err);
    printCLErr(err,__LINE__,__FILE__);

    create_program_with_source(&program, &context, PROGRAM_FILE);

    // add any kernel compiler options to this string
    const char* options="-cl-mad-enable";
    err=clBuildProgram(program, 1, &devices[1], options, NULL, NULL);
    // print build errors
    if(err != CL_SUCCESS)
    {
        perror("error during build");
        size_t log_size=0;
        clGetProgramBuildInfo(program, devices[1], CL_PROGRAM_BUILD_LOG, 0, NULL,
                              &log_size);
        char *program_log=(char*)malloc(log_size+1);
        assert(program_log);
        program_log[log_size]='\0';
        clGetProgramBuildInfo(program, devices[1], CL_PROGRAM_BUILD_LOG,
                              log_size+1, program_log, NULL);
        printf("%s\n", program_log);
        free(program_log);
        exit(1);
    }

    io_queue=clCreateCommandQueue(context, devices[1], CL_QUEUE_PROFILING_ENABLE, &err);
    printCLErr(err,__LINE__,__FILE__);

    compute_queue=clCreateCommandQueue(context, devices[1], CL_QUEUE_PROFILING_ENABLE,
                                       &err);
    printCLErr(err,__LINE__,__FILE__);

	// create Omega kernel
	omega_kernel=clCreateKernel(program, OMEGA_NAME, &err);
	printCLErr(err,__LINE__,__FILE__);

	// create second Omega kernel if defined
	if(strcmp(OMEGA_NAME2, "") != 0){
		omega_kernel2=clCreateKernel(program, OMEGA_NAME2, &err);
		printCLErr(err,__LINE__,__FILE__);
	}

	// Get workgroup size or preferred size
	size_t max_group_size, pref_group_size;
	err = clGetKernelWorkGroupInfo(omega_kernel, devices[1],CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &max_group_size, NULL);
	printf("Max kernel work-group size: %lu\n",max_group_size);	// 256
	err = clGetKernelWorkGroupInfo(omega_kernel, devices[1],CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &pref_group_size, NULL);
	printf("Work-group pref. multiple: %lu\n",pref_group_size);	// 64

	// Get number of work groups / compute units
	err=clGetDeviceInfo(devices[1], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint),
                        &comp_units, NULL);
    printCLErr(err,__LINE__,__FILE__);

	group_size = max_group_size; //pref_group_size; //variate!!

	work_items = group_size * 24;			//variate!!

	printf("Set work-group size: %lu, Set work-items: %lu\n", group_size, work_items);

	cl_ulong total;

    // NOTE: assume device has enough memory
    // TODO: does performance degrade if k is much less than KC?
    cl_ulong a_buffer_size=GPU_BLOCK_MC * GPU_BLOCK_KC * sizeof(inputDataType_x32);
    cl_ulong b_buffer_size=GPU_BLOCK_KC * GPU_BLOCK_NC * sizeof(inputDataType_x32);
    cl_ulong c_buffer_size=GPU_BLOCK_MC * GPU_BLOCK_NC * sizeof(inputDataType_x32);
	total = 2*a_buffer_size + 2*b_buffer_size + 2*c_buffer_size;

    cl_ulong max_alloc=0;
    err=clGetDeviceInfo(devices[1], CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(max_alloc),
                        &max_alloc, NULL);
    printCLErr(err,__LINE__,__FILE__);

    cl_ulong global_mem=0;
    err=clGetDeviceInfo(devices[1], CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(global_mem),
                        &global_mem, NULL);
    printCLErr(err,__LINE__,__FILE__);

	// Omega buffers
	cl_ulong omega_buffer_size= 2 * GPU_BLOCK_MC * GPU_BLOCK_KC * sizeof(float);
	cl_ulong LS_buffer_size= 2 * GPU_BLOCK_MC * sizeof(float);
	cl_ulong RS_buffer_size= 2 * GPU_BLOCK_MC * sizeof(float);
	cl_ulong TS_buffer_size= 2 * GPU_BLOCK_MC * GPU_BLOCK_KC * sizeof(float);
	cl_ulong k_buffer_size= 2 * GPU_BLOCK_MC * sizeof(int);
	cl_ulong m_buffer_size= 2 * GPU_BLOCK_MC * sizeof(int);
	cl_ulong LRkm_buffer_size= 4 * GPU_BLOCK_MC * sizeof(float);

	if(strcmp(OMEGA_NAME, "omega") == 0){
		total += omega_buffer_size + LS_buffer_size + RS_buffer_size + TS_buffer_size + 
			k_buffer_size + m_buffer_size;
	}
	else if(strcmp(OMEGA_NAME, "omega2") == 0){
		LS_buffer_size *= GPU_BLOCK_KC;
		RS_buffer_size *= GPU_BLOCK_KC;
		k_buffer_size *= GPU_BLOCK_KC;
		m_buffer_size *= GPU_BLOCK_KC;
		total += omega_buffer_size + LS_buffer_size + RS_buffer_size + TS_buffer_size + 
			k_buffer_size + m_buffer_size;
	}
	else if(strcmp(OMEGA_NAME, "omega3") == 0){
		total += omega_buffer_size + LRkm_buffer_size + TS_buffer_size;
	}
	else if(strcmp(OMEGA_NAME, "omega4") == 0){
		total += 2 * comp_units * sizeof(float) + LRkm_buffer_size + TS_buffer_size;
	}
	else if(strcmp(OMEGA_NAME, "omega5") == 0){
		LS_buffer_size *= GPU_BLOCK_KC;
		RS_buffer_size *= GPU_BLOCK_KC;
		k_buffer_size *= GPU_BLOCK_KC;
		m_buffer_size *= GPU_BLOCK_KC;
		total += 2 * comp_units * sizeof(float) + LS_buffer_size + RS_buffer_size + TS_buffer_size + 
			k_buffer_size + m_buffer_size;
	}
	else if(strcmp(OMEGA_NAME, "omega6") == 0 || strcmp(OMEGA_NAME, "omega7") == 0){
		omega_buffer_size = work_items * sizeof(float);
		total += 2 * omega_buffer_size + LRkm_buffer_size + TS_buffer_size;
	}
	else if(strcmp(OMEGA_NAME, "omega8") == 0 || strcmp(OMEGA_NAME, "omega9") == 0 || strcmp(OMEGA_NAME, "omega10") == 0){
		omega_buffer_size = LS_buffer_size;
		total += 2 * omega_buffer_size + LRkm_buffer_size + TS_buffer_size;
	}
	else if(strcmp(OMEGA_NAME, "omega11") == 0){
		total += 2 * omega_buffer_size + LRkm_buffer_size + TS_buffer_size;
	}
	else if(strcmp(OMEGA_NAME, "omega13") == 0){
		total += 4 * omega_buffer_size + LRkm_buffer_size + TS_buffer_size;
	}
	else{
		printf("Wrong kernel name\n");
		exit(1);
	}

    if(total > global_mem)
    {
        printf("not enough global storage!\n");
        exit(1);
    }
    if((a_buffer_size > max_alloc) ||
       (b_buffer_size > max_alloc) ||
       (c_buffer_size > max_alloc) ||
	   (omega_buffer_size > max_alloc))
    {
        printf("some buffer is too big!\n");
        exit(1);
    }

    unsigned int i;
    // NOTE: this is what's being passed in for the other GEMM implementation
    rs_c=1;
    cs_c=GPU_BLOCK_MC;
    for(i=0; i < 2; i++)
    {
        // create double buffers for input A matrix
        // printf("a buffer %u\n", i);
        a_buffers[i]=clCreateBuffer(
                context, CL_MEM_READ_ONLY,
                a_buffer_size, NULL, &err
                );
        printCLErr(err,__LINE__,__FILE__);

        // create double buffers for input B matrix
        // printf("b buffer %u\n", i);
        b_buffers[i]=clCreateBuffer(
                context, CL_MEM_READ_ONLY,
                b_buffer_size, NULL, &err
                );
        printCLErr(err,__LINE__,__FILE__);

        // TODO: eventually change this to not just write only??
        // create output buffers for C matrix
        // printf("c buffer %u\n", i);
        c_buffers[i]=clCreateBuffer(
                context, CL_MEM_READ_WRITE,
                c_buffer_size, NULL, &err
                );
        printCLErr(err,__LINE__,__FILE__);

        /* TODO:
           What's better: using an intermediate buffer between C and c_buffer, or
           many calls to C clEnqueueReadBuffer/clEnqueueWriteBuffer (using offset in
           a for loop).
           when the writes are blocking, the intermediate looks better.  try using
           event queueing...
           */
        // c_sub_matrix[i]=calloc(c_buffer_size,1);

        // create kernels
        kernels[i]=clCreateKernel(program, KERNEL_NAME, &err);
        printCLErr(err,__LINE__,__FILE__);

        // TODO: might have to move this inside the loops depending on
        // what parameters we're using (ir, mr, jr, nr, pr ?) can this be
        // substituted by global ids if A is in col major, B is in row major?
        err |= clSetKernelArg(kernels[i], 3, sizeof(cl_mem), &a_buffers[i]);
        // NOTE: don't do B here since that alternates, note later.
        err |= clSetKernelArg(kernels[i], 5, sizeof(cl_mem), &c_buffers[i]);
        err |= clSetKernelArg(kernels[i], 6, sizeof(unsigned int), &rs_c);
        err |= clSetKernelArg(kernels[i], 7, sizeof(unsigned int), &cs_c);
        printCLErr(err,__LINE__,__FILE__);
    }

	if(strcmp(OMEGA_NAME, "omega") == 0 || strcmp(OMEGA_NAME, "omega2") == 0){
		omega_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, omega_buffer_size, NULL, &err); // Should change when read!!!
		printCLErr(err,__LINE__,__FILE__);
		LS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, LS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		RS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, RS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		TS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, TS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		k_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, k_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		m_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, m_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);

		// set kernel arguements for buffers
		err |= clSetKernelArg(omega_kernel, 0, sizeof(cl_mem), &omega_buffer);
		err |= clSetKernelArg(omega_kernel, 1, sizeof(cl_mem), &LS_buffer);
		err |= clSetKernelArg(omega_kernel, 2, sizeof(cl_mem), &RS_buffer);
		err |= clSetKernelArg(omega_kernel, 3, sizeof(cl_mem), &TS_buffer);
		err |= clSetKernelArg(omega_kernel, 4, sizeof(cl_mem), &k_buffer);
		err |= clSetKernelArg(omega_kernel, 5, sizeof(cl_mem), &m_buffer);
		printCLErr(err,__LINE__,__FILE__);
	}
	// else if(strcmp(OMEGA_NAME, "omega3") == 0){
	// 	omega_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, omega_buffer_size, NULL, &err); // Should change when read!!!
	// 	printCLErr(err,__LINE__,__FILE__);
	// 	LRkm_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, LRkm_buffer_size, NULL, &err);
	// 	printCLErr(err,__LINE__,__FILE__);
	// 	TS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, TS_buffer_size, NULL, &err);
	// 	printCLErr(err,__LINE__,__FILE__);

	// 	// set kernel arguements for buffers
	// 	err |= clSetKernelArg(omega_kernel, 0, sizeof(cl_mem), &omega_buffer);
	// 	err |= clSetKernelArg(omega_kernel, 1, sizeof(cl_mem), &LRkm_buffer);
	// 	err |= clSetKernelArg(omega_kernel, 2, sizeof(cl_mem), &TS_buffer);
	// 	printCLErr(err,__LINE__,__FILE__);
	// }
	else if(strcmp(OMEGA_NAME, "omega3") == 0){
		omega_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, omega_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		LR_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, 2*LS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		TS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, TS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		km_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, 2*k_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);

		// set kernel arguements for buffers
		err |= clSetKernelArg(omega_kernel, 0, sizeof(cl_mem), &omega_buffer);
		err |= clSetKernelArg(omega_kernel, 1, sizeof(cl_mem), &LR_buffer);
		err |= clSetKernelArg(omega_kernel, 2, sizeof(cl_mem), &TS_buffer);
		err |= clSetKernelArg(omega_kernel, 3, sizeof(cl_mem), &km_buffer);
		printCLErr(err,__LINE__,__FILE__);
	}
	else if(strcmp(OMEGA_NAME, "omega4") == 0){
		omega_im=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * comp_units, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		index_im=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(unsigned int) * comp_units, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		LRkm_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, LRkm_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		TS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, TS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		// omega_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float), NULL, &err);
		// printCLErr(err,__LINE__,__FILE__);
		// index_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(int), NULL, &err);
		// printCLErr(err,__LINE__,__FILE__);

		// set kernel arguements for buffers
		err |= clSetKernelArg(omega_kernel, 0, sizeof(cl_mem), &omega_im);
		err |= clSetKernelArg(omega_kernel, 1, sizeof(cl_float) * group_size, NULL);
		err |= clSetKernelArg(omega_kernel, 2, sizeof(cl_mem), &index_im);
		err |= clSetKernelArg(omega_kernel, 3, sizeof(cl_uint) * group_size, NULL);
		err |= clSetKernelArg(omega_kernel, 4, sizeof(cl_mem), &LRkm_buffer);
		err |= clSetKernelArg(omega_kernel, 5, sizeof(cl_mem), &TS_buffer);
		// err |= clSetKernelArg(omega_kernel, 8, sizeof(cl_mem), &omega_buffer);
		// err |= clSetKernelArg(omega_kernel, 9, sizeof(cl_mem), &index_buffer);
		printCLErr(err,__LINE__,__FILE__);
	}
	else if(strcmp(OMEGA_NAME, "omega5") == 0){
		omega_im=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * comp_units, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		index_im=clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(unsigned int) * comp_units, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		LS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, LS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		RS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, RS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		TS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, TS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		k_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, k_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		m_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, m_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);

		// set kernel arguements for buffers
		err |= clSetKernelArg(omega_kernel, 0, sizeof(cl_mem), &omega_im);
		err |= clSetKernelArg(omega_kernel, 1, sizeof(cl_float) * group_size, NULL);
		err |= clSetKernelArg(omega_kernel, 2, sizeof(cl_mem), &index_im);
		err |= clSetKernelArg(omega_kernel, 3, sizeof(cl_uint) * group_size, NULL);
		err |= clSetKernelArg(omega_kernel, 4, sizeof(cl_mem), &LS_buffer);
		err |= clSetKernelArg(omega_kernel, 5, sizeof(cl_mem), &RS_buffer);
		err |= clSetKernelArg(omega_kernel, 6, sizeof(cl_mem), &TS_buffer);
		err |= clSetKernelArg(omega_kernel, 7, sizeof(cl_mem), &k_buffer);
		err |= clSetKernelArg(omega_kernel, 8, sizeof(cl_mem), &m_buffer);
		printCLErr(err,__LINE__,__FILE__);
	}
	else if(strcmp(OMEGA_NAME, "omega7") == 0 || strcmp(OMEGA_NAME, "omega8") == 0 || strcmp(OMEGA_NAME, "omega6") == 0 || strcmp(OMEGA_NAME, "omega9") == 0 || strcmp(OMEGA_NAME, "omega10") == 0){
		omega_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, omega_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		index_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, omega_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		LR_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, 2*LS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		TS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, TS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		km_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, 2*k_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);

		// set kernel arguements for buffers
		err |= clSetKernelArg(omega_kernel, 0, sizeof(cl_mem), &omega_buffer);
		err |= clSetKernelArg(omega_kernel, 1, sizeof(cl_mem), &index_buffer);
		err |= clSetKernelArg(omega_kernel, 2, sizeof(cl_mem), &LR_buffer);
		err |= clSetKernelArg(omega_kernel, 3, sizeof(cl_mem), &TS_buffer);
		err |= clSetKernelArg(omega_kernel, 4, sizeof(cl_mem), &km_buffer);
		printCLErr(err,__LINE__,__FILE__);
	}
	else if(strcmp(OMEGA_NAME, "omega11") == 0){
		omega_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, omega_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		index_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, omega_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		LR_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, 2*LS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		TS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, TS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		km_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, 2*k_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);

		// set kernel arguements for buffers
		err |= clSetKernelArg(omega_kernel, 0, sizeof(cl_mem), &omega_buffer);
		err |= clSetKernelArg(omega_kernel, 1, sizeof(cl_mem), &index_buffer);
		err |= clSetKernelArg(omega_kernel, 2, sizeof(cl_mem), &LR_buffer);
		err |= clSetKernelArg(omega_kernel, 3, sizeof(cl_float) * group_size, NULL);
		err |= clSetKernelArg(omega_kernel, 4, sizeof(cl_mem), &TS_buffer);
		err |= clSetKernelArg(omega_kernel, 5, sizeof(cl_mem), &km_buffer);
		err |= clSetKernelArg(omega_kernel, 6, sizeof(cl_int) * group_size, NULL);
		printCLErr(err,__LINE__,__FILE__);
	}
	else if(strcmp(OMEGA_NAME, "omega13") == 0){
		omega_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, omega_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		index_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, omega_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		LR_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, 4*LS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		TS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, TS_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);
		km_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, 4*k_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);

		// set kernel arguements for buffers
		err |= clSetKernelArg(omega_kernel, 0, sizeof(cl_mem), &omega_buffer);
		err |= clSetKernelArg(omega_kernel, 1, sizeof(cl_mem), &index_buffer);
		err |= clSetKernelArg(omega_kernel, 2, sizeof(cl_mem), &LR_buffer);
		err |= clSetKernelArg(omega_kernel, 3, sizeof(cl_mem), &TS_buffer);
		err |= clSetKernelArg(omega_kernel, 4, sizeof(cl_mem), &km_buffer);
		printCLErr(err,__LINE__,__FILE__);
	}
	else{
		printf("Wrong kernel name\n");
	}

	if(strcmp(OMEGA_NAME2, "omega10") == 0){
		index_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, k_buffer_size, NULL, &err);
		printCLErr(err,__LINE__,__FILE__);

		// set kernel arguements for buffers second Omega kernel
		err |= clSetKernelArg(omega_kernel2, 0, sizeof(cl_mem), &omega_buffer);
		err |= clSetKernelArg(omega_kernel2, 1, sizeof(cl_mem), &index_buffer);
		err |= clSetKernelArg(omega_kernel2, 2, sizeof(cl_mem), &LR_buffer);
		err |= clSetKernelArg(omega_kernel2, 3, sizeof(cl_mem), &TS_buffer);
		err |= clSetKernelArg(omega_kernel2, 4, sizeof(cl_mem), &km_buffer);
		printCLErr(err,__LINE__,__FILE__);
	}
	// Get device max workgroup size
	// size_t group_size;
	// err=clGetDeviceInfo(devices[1], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t),
    //                     &group_size, NULL);
    // printCLErr(err,__LINE__,__FILE__);
	// printf("max workgroup size: %lu\n",group_size);

	// Get constant memory info
	// cl_uint max_args;
	// cl_ulong max_const;
	// err=clGetDeviceInfo(devices[1], CL_DEVICE_MAX_CONSTANT_ARGS, sizeof(cl_uint),
    //                     &max_args, NULL);		// 16
    // printCLErr(err,__LINE__,__FILE__);
	// err=clGetDeviceInfo(devices[1], CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(cl_ulong),
    //                     &max_const, NULL);		// 1503238400
    // printCLErr(err,__LINE__,__FILE__);

	// printf("Arg: %u, Size: %lu\n",max_args,max_const);

	// Get local memory info
	// cl_ulong max_local;
	// err=clGetDeviceInfo(devices[1], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong),
    //                     &max_local, NULL);		// 32768
    // printCLErr(err,__LINE__,__FILE__);

	// printf("Size: %lu\n",max_local);

	// Get timer resolution info
	// cl_ulong resolution;
	// clGetDeviceInfo(devices[1], CL_DEVICE_PROFILING_TIMER_RESOLUTION, sizeof(cl_ulong), 
	// 						&resolution, NULL);		// 0
	// printCLErr(err,__LINE__,__FILE__);

	// printf("Res: %lu\n",resolution);
	// float * test = NULL;
}

void gpu_release(void)
{
    // Omega
	clReleaseKernel(omega_kernel);
	clReleaseMemObject(omega_buffer);
	clReleaseMemObject(LS_buffer);
	clReleaseMemObject(RS_buffer);
	clReleaseMemObject(TS_buffer);
	clReleaseMemObject(k_buffer);
	clReleaseMemObject(m_buffer);
	clReleaseMemObject(LRkm_buffer);
	clReleaseMemObject(LR_buffer);
	clReleaseMemObject(km_buffer);
	clReleaseMemObject(omega_im);
	clReleaseMemObject(index_im);
	clReleaseMemObject(index_buffer);

	for(int i=0; i < 2; i++)
    {
        clReleaseKernel(kernels[i]);
        clReleaseMemObject(c_buffers[i]);
        clReleaseMemObject(b_buffers[i]);
        clReleaseMemObject(a_buffers[i]);
        // free(c_sub_matrix[i]);
        //        for(int j=0; j < 4; j++)
        //            clReleaseEvent(events[(2*i)+j]);
    }
    // clReleaseCommandQueue(compute_queue);
    clReleaseCommandQueue(io_queue);
    clReleaseCommandQueue(compute_queue);
    clReleaseProgram(program);
    clReleaseContext(context);
    free(devices);
    free(platforms);
}