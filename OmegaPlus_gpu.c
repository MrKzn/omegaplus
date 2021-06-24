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
 *  !!! ADD MY OWN DETAILS? (REINOUT) !!!
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

/*   ---  GPU Omega functions  ---   */
void computeOmega_gpuF1(float * omegas, float * LR, int * km, float * T, int in_out_cnt, int inner_cnt, unsigned int total){
	// static cl_ulong p_start, p_end, p_total=0;
	int err=0;

	// //set kernel arguments
	err = clSetKernelArg(omega_kernel2, 4, sizeof(int), &inner_cnt);
	printCLErr(err,__LINE__,__FILE__);

	const size_t local = LOCAL_2;
	const size_t global = total + (LOCAL_2 - (total & (LOCAL_2 - 1)));

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
			total*sizeof(float), T,
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

	// clWaitForEvents(1, &events[1]);

    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
    //                         &p_start, NULL);
	// printCLErr(err,__LINE__,__FILE__);
    // err=clGetEventProfilingInfo(events[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
    //                         &p_end, NULL);
	// printCLErr(err,__LINE__,__FILE__);

	// p_total = p_end - p_start;

	// printf("%lu\n",p_total);

	//read back omega values in omega buffer
	err=clEnqueueReadBuffer(
			io_queue, omega_buffer, CL_FALSE, 0,
			total*sizeof(float), omegas,
			1, &events[1], NULL
			);
	printCLErr(err,__LINE__,__FILE__);

    clFinish(io_queue);
}

void computeOmega_gpuF2(float * omegas, unsigned int * indexes, float * LR, int * km, float * T, int in_out_cnt, int outer, int inner, unsigned int total){
	// static cl_ulong p_start, p_end, p_total=0;

	int err=0;
	const size_t local = group_size;
	const size_t global = work_items; 

	// set kernel arguments
	err |= clSetKernelArg(omega_kernel, 5, sizeof(cl_int), &outer);
	err |= clSetKernelArg(omega_kernel, 6, sizeof(cl_int), &inner);
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
			total*sizeof(float), T,
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

	// p_total = p_end - p_start;

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

    clFinish(io_queue);
}

void computeOmegaValues_gpuF (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	// static double mtime0, mtime1, mtimetot = 0;
	float tmpW, maxW=0.0;
	static float * omegas = NULL, * LR = NULL, * T = NULL;

	unsigned int * indexes = NULL, index=0;
	static int * km = NULL;
	int i, j, maxLeftIndex=0, maxRightIndex=0, inner_work, outer_work, work_groups,
	
	total, work_total, outer_cnt, inner_cnt, in_out_cnt, inner_i=0, Lk_i=0, Rm_i=0,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;
	
	outer_work = outer_cnt = leftMinIndex - leftMaxIndex + 1;
	inner_work = inner_cnt = rightMaxIndex - rightMinIndex + 1;
	total = outer_cnt * inner_cnt;
	
	if(total > 20000)	// Threshold
	{
		outer_work += group_size - (outer_cnt & (group_size - 1));
		inner_work += group_size - (inner_cnt & (group_size - 1));
		work_total = outer_work * inner_work;
		in_out_cnt = outer_work + inner_work;

		// work_groups = work_total / (group_size * group_size);
		// work_items = group_size * work_groups;
        work_items = work_total / group_size;

        if(work_items > max_omegas || total > max_TS || in_out_cnt > max_LRkm){
            printf("Grid point %d is skipped due to insufficient GPU memory\n",omegaIndex);
            return;
        }

		omegas = malloc(sizeof(*omegas)*work_items);
		indexes = malloc(sizeof(*indexes)*work_items);
		LR = malloc(sizeof(*LR) * in_out_cnt);
		km = malloc(sizeof(*km) * in_out_cnt);
		T = malloc(sizeof(*T) * work_total);

		if(omegas==NULL || indexes==NULL || LR==NULL || km==NULL || T==NULL)
			printf("MALLOC error\n");

		Rm_i = outer_work;

		for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
		{
			LR[Lk_i] = correlationMatrix[omegaSNIPIndex][i];				// LSs

			km[Lk_i] = omegaSNIPIndex - i + 1;							// ks

			for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
			{
				if(!Lk_i)
				{
					LR[Rm_i] = correlationMatrix[j][omegaSNIPIndex+1];	// RSs

					km[Rm_i] = j - omegaSNIPIndex;						// ms

					Rm_i++;
				}
				
				T[inner_i] = correlationMatrix[j][i];
				inner_i++;
			}
			Lk_i++;
			inner_i = Lk_i * inner_work;
		}
		
		for(i=outer_work+inner_cnt;i<in_out_cnt;i++){
			LR[i] = 1;
			km[i] = 0;
		}
		for(i=outer_cnt;i<outer_work;i++){
			LR[i] = 1;
			km[i] = 0;
		}

		// mtime0 = gettime();
		computeOmega_gpuF2(omegas, indexes, LR, km, T, in_out_cnt, outer_work, inner_work, work_total);
		// mtime1 = gettime();
		// mtimetot = mtime1 - mtime0;
		// printf("%f\n",mtimetot);

		for(i=0;i<work_items;i++)
		{
			tmpW = omegas[i];
			if(tmpW>maxW)
			{
				maxW = tmpW;
				index = i;
			}
		}

		maxLeftIndex = (leftMinIndex - (indexes[index]/inner_work)) + omega[omegaIndex].leftIndex;
		maxRightIndex = (rightMinIndex + (indexes[index]%inner_work)) + omega[omegaIndex].leftIndex;
		
		free(indexes);
	}
	else
	{
		in_out_cnt = outer_cnt + inner_cnt;

        if(work_items > max_omegas || total > max_TS || in_out_cnt > max_LRkm){
            printf("Grid point %d is skipped due to insufficient GPU memory\n",omegaIndex);
            return;
        }

		omegas = malloc(sizeof(*omegas)*total);
		LR = malloc(sizeof(*LR)*in_out_cnt);
		km = malloc(sizeof(*km)*in_out_cnt);
		T = malloc(sizeof(*T)*total);

		if(omegas==NULL || LR==NULL || km==NULL || T==NULL)
			printf("MALLOC error\n");

		Lk_i = inner_cnt;
		
		for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
		{
			LR[Lk_i] = correlationMatrix[omegaSNIPIndex][i];			// LSs

			km[Lk_i] = omegaSNIPIndex - i + 1;							// ks

			for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
			{
				if(!(Lk_i-inner_cnt))
				{
					LR[Rm_i] = correlationMatrix[j][omegaSNIPIndex+1];	// RSs

					km[Rm_i] = j - omegaSNIPIndex;						// ms

					Rm_i++;
				}
				
				T[inner_i] = correlationMatrix[j][i];

				inner_i++;
			}
			Lk_i++;
		}

		// mtime0 = gettime();
		computeOmega_gpuF1(omegas, LR, km, T, in_out_cnt, inner_cnt, total);
		// mtime1 = gettime();
		// mtimetot = mtime1 - mtime0;
		// printf("%f\n",mtimetot);

		for(i=0;i<total;i++)
		{
			tmpW = omegas[i];
			if(tmpW>maxW)
			{
				maxW = tmpW;
				index = i;
			}
		}

		maxLeftIndex = (leftMinIndex - (index/inner_cnt)) + omega[omegaIndex].leftIndex;
		maxRightIndex = (rightMinIndex + (index%inner_cnt)) + omega[omegaIndex].leftIndex;
	}
	
	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;
	
	free(omegas);
	free(LR);
	free(km);
	free(T);
}

void computeOmegas_gpu (alignment_struct * alignment, omega_struct * omega, int omegaIndex, void * threadData, cor_t ** correlationMatrix)
{
	computeOmegaValues_gpuF (omega, omegaIndex, alignment->correlationMatrix, NULL);
}


/*   ---  GPU qLD functions  ---   */
static unsigned int round_up_mult(unsigned int x, unsigned int mult) {
    return ((x + mult - 1) / mult) * mult;
}

void GPU_Pack_A(inputDataType_x32 *A,
        unsigned int lda,
        DOUBLE *A_pack,
        unsigned int m,
        unsigned int k)
{
    DOUBLE *A_pack_local;
    unsigned int m_alg;
    //#pragma omp  parallel for num_threads(*n_threads) private(A_pack_local)
    //  #pragma omp parallel
    //  #pragma omp single
    ////    #pragma omp taskloop private(A_pack_local) num_tasks((*n_threads)) label(gemm_pack)
    //  #pragma omp taskloop private(A_pack_local) grainsize(16) label(gemm_pack)
    for(unsigned int ic=0;ic<m;ic+=GPU_BLOCK_MR)
    {
        A_pack_local=&A_pack[ic*k];
        m_alg=min(GPU_BLOCK_MR,m-ic);
        for(unsigned int pc=0;pc<k;pc++)
        {
            for(unsigned int ir=0;ir<m_alg;ir++)
            {
                A_pack_local[0]=A[(ic+ir)+pc*lda]; //auto xtypaei...
                A_pack_local++;
            }
        }
    }
}

void GPU_Pack_B(inputDataType_x32 *B,
        unsigned int ldb,
        DOUBLE *B_pack,
        unsigned int k,
        unsigned int n)
{
    DOUBLE *B_pack_local;
    unsigned int n_alg;
    //#pragma omp parallel for num_threads(*n_threads) private(B_pack_local)
    //  #pragma omp parallel
    //  #pragma omp single
    ////    #pragma omp taskloop private(B_pack_local) num_tasks((*n_threads)) label(gemm_pack)
    //  #pragma omp taskloop private(B_pack_local) grainsize(16) label(gemm_pack)
    for(unsigned int jc=0;jc<n;jc+=GPU_BLOCK_NR)
    {
        B_pack_local=&B_pack[jc*k];
        n_alg=min(GPU_BLOCK_NR,n-jc);
        for(unsigned int pc=0;pc<k;pc++)
        {
            for(unsigned int jr=0;jr<n_alg;jr++)
            {
                B_pack_local[0]=B[pc+jc*ldb+jr*ldb];
                B_pack_local++;
            }
        }
    }
}

void gpu_gemm(unsigned int m,
        unsigned int n,
        unsigned int k,
        inputDataType_x32 *A,
        unsigned int lda,
        inputDataType_x32 *B,
        unsigned int ldb,
        inputDataType_x32 *C,
        unsigned int ldc,
        void *Ac_pack_v,
        void *Bc_pack_v)
{
    static cl_ulong p_start, p_end, p_total=0;
    inputDataType_x32 *Ac, *Bc;
    inputDataType_x32 *Cc;
    //unsigned int *Ar, *Br;
    //unsigned int *Cr;
    //DOUBLE beta;
    unsigned int i=0;
    int err;
    const unsigned int work_dim=2;
    size_t local[work_dim];
    local[0]=LOCAL_0;
    local[1]=LOCAL_1;
    size_t global[work_dim];

    // // FOR TESTING WRITE BUFFERS
    // int pm;     //return value used for assert
    // void *A_test=NULL, *B_test=NULL;
    // pm=posix_memalign(&(A_test), 4096,
    //         GPU_BLOCK_MC*GPU_BLOCK_KC*sizeof(inputDataType_x32));
    // assert(!pm);
    // pm=posix_memalign(&(B_test), 4096,
    //         GPU_BLOCK_MC*GPU_BLOCK_KC*sizeof(inputDataType_x32));
    // assert(!pm);

    for(unsigned int jc=0; jc<n; jc+=GPU_BLOCK_NC)
    {
        unsigned int n_alg=min(GPU_BLOCK_NC,n-jc);

        for(i=0; i < 2; i++)
        {
            err=clSetKernelArg(kernels[i], 2, sizeof(inputDataType_x32), &n_alg);
            printCLErr(err,__LINE__,__FILE__);
        }

        unsigned int pc_iter=0;
        for(unsigned int pc=0; pc<k; pc+=GPU_BLOCK_KC)
        {
            unsigned int k_alg=min(GPU_BLOCK_KC,k-pc);

            for(i=0; i < 2; i++)
            {
                err=clSetKernelArg(kernels[i], 0, sizeof(inputDataType_x32), &k_alg);
                printCLErr(err,__LINE__,__FILE__);
            }

            //beta=*betap;

            Bc=&B[pc+jc*ldb];
            GPU_Pack_B(Bc, ldb, Bc_pack_v, k_alg, n_alg);

            // TODO: double buffer B, don't just use b_buffers[0]
            // printf("\nwriting b (k: %u, n: %u)\n", k_alg, n_alg);

            // NOTE: getting CL_MEM_OBJECT_ALLOCATION_FAILURE when
            // KC*NC*sizeof(uint) is too big... but not more than max_alloc??

            // NOTE: don't have to use an event here since
            //        io_queue is processed in-order
            err=clEnqueueWriteBuffer(
                    io_queue, b_buffers[pc_iter % 2], CL_FALSE, 0,
                    k_alg*n_alg*sizeof(inputDataType_x32), Bc_pack_v,
                    0, NULL, &events[(pc_iter % 2)*4]
                    );
            printCLErr(err,__LINE__,__FILE__);

            // // ADDED
            // // Check if "Bc_pack_v" is correctly written by reading from GPU
            // err=clEnqueueReadBuffer(
            //         io_queue, b_buffers[pc_iter % 2], CL_FALSE, 0,
            //         k_alg*n_alg*sizeof(inputDataType_x32), B_test,
            //         1, &events[(pc_iter % 2)*4], NULL
            //         );
            // printCLErr(err,__LINE__,__FILE__);

            // for(int kl = 0; kl < k_alg*n_alg; kl++){
            //     if(((uint32_t*)B_test)[kl] != ((uint32_t*)Bc_pack_v)[kl])
            //         printf("FAULT\n");
            // }

            // printf("write done\n");
            // set both kernels to use this iteration's b_buffer
            for(i=0; i < 2; i++)
            {
                err=clSetKernelArg(
                        kernels[i], 4, sizeof(cl_mem), &b_buffers[pc_iter % 2]
                        );
                printCLErr(err,__LINE__,__FILE__);
            }

            unsigned int ic_iter=0;
            // write A, write C, kernel, read C (*2 for double buffer)

            for(unsigned int ic=0; ic<m; ic+=GPU_BLOCK_MC)
            {

                unsigned int m_alg=min(GPU_BLOCK_MC,m-ic);

                // m_alg is MC until it is the tail
                // global is going to be derived from m_alg/n_alg.
                // at - least as big as NR/MR (1 block).
                // ceiling division
                // TODO: figure this out for larger MR/NR
                // global[0]=((max(m_alg, LOCAL_0) + LOCAL_0 - 1) / LOCAL_0) * LOCAL_0 / GPU_BLOCK_MR;
                global[0]=round_up_mult(
                        round_up_mult(m_alg, LOCAL_0), BLOCK_SIZE_X
                        ) / BLOCK_SIZE_X * LOCAL_0;
                // (((n_alg + LOCAL_0 - 1)/LOCAL_0) * LOCAL_0) * LOCAL_0 / BLOCK_SIZE_X;
                global[1]=round_up_mult(
                        round_up_mult(n_alg, LOCAL_1), BLOCK_SIZE_Y
                        ) / BLOCK_SIZE_Y * LOCAL_1;

                global[0]=max(global[0], local[0]);
                global[1]=max(global[1], local[1]);

                // global[0]=round_up_mult(
                //   round_up_mult(n_alg, LOCAL_0) * LOCAL_0, BLOCK_SIZE_X
                // ) / BLOCK_SIZE_X;
                // // (((n_alg + LOCAL_0 - 1)/LOCAL_0) * LOCAL_0) * LOCAL_0 / BLOCK_SIZE_X;
                // global[1]=round_up_mult(
                //   round_up_mult(m_alg, LOCAL_1) * LOCAL_1, BLOCK_SIZE_Y
                // ) / BLOCK_SIZE_Y;



                //(((m_alg + LOCAL_1 - 1)/LOCAL_1) * LOCAL_1) * LOCAL_1 / BLOCK_SIZE_Y;
                // global[1]=((max(n_alg, LOCAL_1) + LOCAL_1 - 1) / LOCAL_1) * LOCAL_1 / GPU_BLOCK_NR;

                err=clSetKernelArg(
                        kernels[ic_iter % 2], 1, sizeof(inputDataType_x32), &m_alg
                        );
                printCLErr(err,__LINE__,__FILE__);

                Ac=&A[ic+pc*lda];
                GPU_Pack_A(Ac,lda,Ac_pack_v,m_alg,k_alg);

                Cc=&C[ic+jc*ldc];

                for(i=0; i < n_alg; ++i)
                {
                    // for(i=0; i < GPU_BLOCK_NC; ++i) {
                    err=clEnqueueWriteBuffer(
                            io_queue, c_buffers[ic_iter % 2], CL_FALSE,
                            i*cs_c*sizeof(inputDataType_x32),
                            m_alg*sizeof(inputDataType_x32), &Cc[i*ldc],
                            //GPU_BLOCK_MC*sizeof(unsigned int), &Cc[i*ldc],
                            0, NULL, NULL
                            );
                    printCLErr(err,__LINE__,__FILE__);
                }

                // printf("writing a\n");

                err=clEnqueueWriteBuffer(
                        io_queue, a_buffers[ic_iter % 2], CL_FALSE, 0,
                        m_alg*k_alg*sizeof(inputDataType_x32), Ac_pack_v,
                        0, NULL, &events[(ic_iter % 2)*4]
                        );
                printCLErr(err,__LINE__,__FILE__);

                // // ADDED
                // // Check if "Ac_pack_v" is correctly written by reading from GPU
                // err=clEnqueueReadBuffer(
                //         io_queue, a_buffers[ic_iter % 2], CL_FALSE, 0,
                //         k_alg*n_alg*sizeof(inputDataType_x32), A_test,
                //         1, &events[(ic_iter % 2)*4], NULL
                //         );
                // printCLErr(err,__LINE__,__FILE__);

                // for(int kl = 0; kl < k_alg*n_alg; kl++){
                //     if(((uint32_t*)A_test)[kl] != ((uint32_t*)Ac_pack_v)[kl])
                //         printf("FAULT\n");
                // }

                // assuming C starts cleared, but we want to += if iterating over k (pc)
                // for(i=0; i < n_alg; ++i) {
                //   memcpy(
                //     &(c_sub_matrix[ic_iter % 2][i*cs_c]), &Cc[i*ldc], m_alg*sizeof(unsigned int)
                //   );
                // }
                //
                // err=clEnqueueWriteBuffer(
                //   io_queue, c_buffers[ic_iter % 2], CL_FALSE, 0,
                //   GPU_BLOCK_MC*GPU_BLOCK_NC*sizeof(unsigned int),
                //                                             c_sub_matrix[ic_iter % 2],
                //   0, NULL, &events[(ic_iter % 2)*4 + 1]
                // );
                // printCLErr(err,__LINE__,__FILE__);


                // printf(
                //   "\nm: %u, n: %u, k: %u, ic: %u, jc: %u, loc: %lu,%lu, glob: %lu,%lu\n",
                //   m_alg, n_alg, k_alg, ic, jc, local[0], local[1], global[0], global[1]
                // );

                err=clEnqueueNDRangeKernel(
                        io_queue, kernels[ic_iter % 2], work_dim, NULL, global, local,	// WAS COMPUTE_QUEUE
                        1, &events[(ic_iter % 2)*4], &events[(ic_iter % 2)*4 + 2]
                        );
                printCLErr(err,__LINE__,__FILE__);

                clWaitForEvents(1, &events[(ic_iter % 2)*4 + 2]);   // With waiting for event, qLD is much faster

                err=clGetEventProfilingInfo(events[(ic_iter % 2)*4 + 2], CL_PROFILING_COMMAND_START, sizeof(cl_ulong),
                                        &p_start, NULL);
                printCLErr(err,__LINE__,__FILE__);
                err=clGetEventProfilingInfo(events[(ic_iter % 2)*4 + 2], CL_PROFILING_COMMAND_END, sizeof(cl_ulong),
                                        &p_end, NULL);
                printCLErr(err,__LINE__,__FILE__);

                p_total += p_end - p_start;

                // err=clEnqueueReadBuffer(
                //   io_queue, c_buffers[ic_iter % 2], CL_TRUE, 0,
                //   GPU_BLOCK_MC*GPU_BLOCK_NC*sizeof(unsigned int),
                //                                               c_sub_matrix[ic_iter % 2],
                //   1, &events[(ic_iter % 2)*4 + 2], &events[(ic_iter % 2)*4 + 3]
                // );
                // printCLErr(err,__LINE__,__FILE__);

                // printf("\n");
                // print_matrix(c_sub_matrix[0], GPU_BLOCK_MC, GPU_BLOCK_NC);
                // printf("\n");

                // add event timing (since kernel is done -- TODO: make this callback too)

                for(i=0; i < n_alg; ++i)
                {
                    // for(i=0; i < GPU_BLOCK_NC; ++i) {
                    err=clEnqueueReadBuffer(
                            io_queue, c_buffers[ic_iter % 2], CL_FALSE,
                            i*cs_c*sizeof(inputDataType_x32),
                            m_alg*sizeof(inputDataType_x32), &Cc[i*ldc],
                            // GPU_BLOCK_MC*sizeof(unsigned int), &Cc[i*ldc],
                            0, NULL, NULL
                            );
                    printCLErr(err,__LINE__,__FILE__);
                }

                //                getc(stdin); //used to pause program for monitoring - MPAMPIS -

                clFinish(io_queue);
                // clFinish(compute_queue);

                // TODO: do this in a callback function and make read C asynchronous
                // n_alg rows, each row ldc or cs_c
                // for(i=0; i < n_alg; ++i) {
                //   memcpy(
                //     &Cc[i*ldc], &(c_sub_matrix[ic_iter % 2][i*cs_c]),
                //     m_alg*sizeof(unsigned int)
                //   );
                // }

                ic_iter += 1;
            }
            pc_iter += 1;
        }
    }
    printf("%lu\n",p_total);
}

void mlt_gpu(unsigned int m,
        unsigned int k,
        inputDataType_x32 *A,
        inputDataType_x32 *tableA)
{
    for(unsigned int i=0;i<m;i++)
    {
        for(unsigned int j=0;j<k;j++)
        {
            ((inputDataType_x32*)A)[j*m + i]=tableA[i*k + j];
        }
    }
}

uint32_t * correlate_gpu(uint32_t* tableA,
               int tableAsize,
               int compressed_snp_size)
{
    int m=tableAsize, n=tableAsize, k=compressed_snp_size; 
    long long int i;
    long long int tableCsize = (long long int)m*(long long int)n;
    int pm;     //return value used for assert

    void *Ac_pack_v=NULL, *Bc_pack_v=NULL, *A=NULL, *C=NULL;
    pm=posix_memalign(&(Ac_pack_v), 4096,
            12*GPU_BLOCK_MC*GPU_BLOCK_KC*sizeof(inputDataType_x32));
    assert(!pm);
    pm=posix_memalign(&(Bc_pack_v), 4096,
            GPU_BLOCK_KC*GPU_BLOCK_NC*sizeof(inputDataType_x32));
    assert(!pm);

    pm=posix_memalign(&A, 4096, m*k*sizeof(inputDataType_x32) +
            m*k*sizeof(inputDataType_x32)%4096);
    assert(!pm);
    pm=posix_memalign(&C, 4096, tableCsize*sizeof(inputDataType_x32) +
            tableCsize*sizeof(inputDataType_x32)%4096);
    assert(!pm);

    for(i=0;i<tableCsize;i++)
    {
        ((inputDataType_x32*)C)[i]=0;
    }

    mlt_gpu(m, k, A, tableA);

    gpu_gemm(m,
            n,
            k,
            A,
            m,
            tableA,
            k,
            C,
            m,
            Ac_pack_v,
            Bc_pack_v);

    // int l,o;
	// unsigned int row=0,col=0,tx=0;
	
	// for(o=0;o<m/group_size;o++){
	// 	tx=o*group_size;
	// 	for(i=1;i<group_size;i++) // i is row indexing
	// 	{
	// 		row = (i+tx)*m;
	// 		for(l=i-1;l>=0;l--){ // l is col indexing
	// 			col = l+tx;
	// 			if(((uint32_t*)C)[row+col]>0)
	// 				printf("qLD row: %llu, col: %u, val: %u\n",(i+tx),col,((uint32_t*)C)[row+col]);
	// 		}
	// 	}
	// }

	free(Ac_pack_v);
	free(Bc_pack_v);
    free(A);
    // free(C);

	return ((uint32_t*)C);
}

/*   ---  GPU Correlation functions  ---   */
unsigned int precomputed16_bitcountGPU (unsigned int n)
{
	/* works only for 32-bit unsigned int*/	    
	return bits_in_16bits [n         & 0xffffu]	
	    +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}

cor_t computeCorrelationValueBIN_RSQUAREGPU(int sequences, unsigned int * accumXvec)
{	
	cor_t result;
	cor_t sequencesTotal = (cor_t) sequences;

	if (sequencesTotal==0)
		return 0.;

	cor_t numerator=0.0, 
              denominator=0.0;

	cor_t denomA = (cor_t)accumXvec[0]/sequencesTotal;
	cor_t denomB = (cor_t)accumXvec[1]/sequencesTotal;
	cor_t denomC = (cor_t)accumXvec[2]/sequencesTotal;
	cor_t denomD = (cor_t)accumXvec[3]/sequencesTotal;

	denominator = denomA * denomB * denomC * denomD;

	if (denominator==0.0)
		return 0.0;

	numerator  = ((cor_t)accumXvec[4])/sequencesTotal;
	numerator -= ((cor_t)accumXvec[1]/(cor_t)sequences) * ((cor_t)accumXvec[3]/(cor_t)sequences);
	numerator *= numerator;

	result = (sequences * numerator) /denominator;

	return result;
}

cor_t computeCorrelationValueBIN_DOMGPU(int sequences, unsigned int * accumXvec)
{
	cor_t result;

	cor_t m1  = (cor_t)accumXvec[1];
	cor_t m11 = (cor_t)accumXvec[4];
	cor_t m12 = m1 - m11;
	cor_t m2  = (cor_t)accumXvec[3];
	cor_t m21 = m2 - m11;
	cor_t m22 = (cor_t)accumXvec[0] - m21;

	if(!(accumXvec[2]-m12==m22))	
	{
		printf("Equality Failed! (%f - %f)\n",accumXvec[2]-m12,m22);
		assert(accumXvec[2]-m12==m22);
	}

	m11/=sequences;
	m22/=sequences;
	m12/=sequences;
	m21/=sequences;
	m1/=sequences;
	m2/=sequences;

	if ( ( (m1 >= 0.5) && (m2 >= 0.5) ) || 	( (m1 < 0.5) && (m2 < 0.5 ) ) )
		result =  m11*m22 - m12*m21;
	else if( ( (m1 < 0.5) && (m2 >= 0.5 ) )|| ( (m1 >= 0.5) && (m2 < 0.5) ) )
		result =  m21*m12 - m11*m22;
	else{
		fprintf(stderr, "m1: %e, m2: %e, m11: %e, m12: %e, m21: %e, m22: %e\n", m1, m2, m11, m12, m21, m22);
		assert(0);
	}

	return result;
}

cor_t computeCorrelationValueBIN_JUSTDGPU(int sequences, unsigned int * accumXvec)
{
	cor_t result;

	cor_t m1  = (cor_t)accumXvec[1];
	cor_t m11 = (cor_t)accumXvec[4];
	cor_t m12 = m1 - m11;
	cor_t m2  = (cor_t)accumXvec[3];
	cor_t m21 = m2 - m11;
	cor_t m22 = (cor_t)accumXvec[0] - m21;

	if(!(accumXvec[2]-m12==m22))	
	{
		printf("Equality Failed! (%f - %f)\n",accumXvec[2]-m12,m22);
		assert(accumXvec[2]-m12==m22);
	}

	m11/=sequences;
	m22/=sequences;
	m12/=sequences;
	m21/=sequences;
	m1/=sequences;
	m2/=sequences;

	result = m11 - m1*m2;

	return result;
}

cor_t computeCorrelationValueBIN_DOM2GPU(int sequences, unsigned int * accumXvec)
{
	cor_t result;

	cor_t m1  = (cor_t)accumXvec[1];
	cor_t m11 = (cor_t)accumXvec[4];
	cor_t m12 = m1 - m11;
	cor_t m2  = (cor_t)accumXvec[3];
	cor_t m21 = m2 - m11;
	cor_t m22 = (cor_t)accumXvec[0] - m21;

	if(!(accumXvec[2]-m12==m22))	
	{
		printf("Equality Failed! (%f - %f)\n",accumXvec[2]-m12,m22);
		assert(accumXvec[2]-m12==m22);
	}

	m11/=sequences;
	m22/=sequences;
	m12/=sequences;
	m21/=sequences;
	m1/=sequences;
	m2/=sequences;

	if ( ( (m1 >= 0.5) && (m2 >= 0.5) ) || 	( (m1 < 0.5) && (m2 < 0.5 ) ) )
		result =  m11*m22 - m12*m21;
	else if( ( (m1 < 0.5) && (m2 >= 0.5 ) )|| ( (m1 >= 0.5) && (m2 < 0.5) ) )
		result =  m21*m12 - m11*m22;
	else{
		fprintf(stderr, "m1: %e, m2: %e, m11: %e, m12: %e, m21: %e, m22: %e\n", m1, m2, m11, m12, m21, m22);
		assert(0);
	}

  	cor_t denom = (m1*m2*(1.-m1)*(1.-m2));

  	result = abs(result)/denom;

	return result;
}

cor_t computeCorrelationValueBINGPU(int sequences, unsigned int * accumXvec)
{
	if(linkage_disequilibrium == RSQUARE)
		return computeCorrelationValueBIN_RSQUAREGPU(sequences, accumXvec);

	if(linkage_disequilibrium == DOM)
		return computeCorrelationValueBIN_DOMGPU(sequences, accumXvec);

	if(linkage_disequilibrium == ABSDOM)
		return ABS(computeCorrelationValueBIN_DOMGPU(sequences, accumXvec));

	if(linkage_disequilibrium == JUSTD)
		return computeCorrelationValueBIN_JUSTDGPU(sequences, accumXvec);

	if(linkage_disequilibrium == ABSD)
		return ABS(computeCorrelationValueBIN_JUSTDGPU(sequences, accumXvec));

	if(linkage_disequilibrium == ABSDOM2)
		return ABS(computeCorrelationValueBIN_DOM2GPU(sequences, accumXvec));

	assert(linkage_disequilibrium==999);

	return 0.0;
}

void count01CombsGPU (int total, unsigned int inputL, unsigned int inputR, unsigned int * accumXvec)
{
	int limit=ENTRIES_PER_INT;
	
	if (total<limit)
		limit = total;

	int t1 = precomputed16_bitcountGPU (inputL);
	int t2 = precomputed16_bitcountGPU (inputR);

	// unsigned int combAND = inputL & inputR;

	accumXvec[1] += t1;  
	accumXvec[3] += t2;

	accumXvec[0] += limit - t1;
	accumXvec[2] += limit - t2;

	//accumXvec[4] +=  precomputed16_bitcountGPU (combAND);	
}

int count01GAPCombsGPU (unsigned int inputL, unsigned int inputR, unsigned int inputLVld, unsigned int inputRVld, unsigned int * accumXvec)
{
	
	unsigned int validComb = inputLVld & inputRVld;	
	unsigned int inputUL = inputL & validComb;
	unsigned int inputUR = inputR & validComb;	

	int validSampleSize = precomputed16_bitcountGPU (validComb);
	
	int t1 = precomputed16_bitcountGPU (inputUL); 
	int t2 = precomputed16_bitcountGPU (inputUR); 

	// unsigned int combAND = inputUL&inputUR;

	accumXvec[1] += t1;  
	accumXvec[3] += t2;

	accumXvec[0] += validSampleSize - t1;
	accumXvec[2] += validSampleSize - t2;

	//accumXvec[4] +=  precomputed16_bitcountGPU (combAND);

	return validSampleSize;	
}

cor_t computePairwiseCorrelationBINGPU (alignment_struct * alignment, int s_i, int s_j, uint32_t * qLD_res)
{	
	int i, total;
	unsigned int accumXvec[5];
	unsigned int s_i_val, s_j_val;

	total = alignment->sequences;

	for(i=0;i<5;i++)
		accumXvec[i]=0;
	
	for(i=0;i<alignment->siteSize;i++)
	{
		s_i_val = alignment->compressedArrays[0][s_i*alignment->siteSize+i];
		s_j_val = alignment->compressedArrays[0][s_j*alignment->siteSize+i];

		count01CombsGPU (total, s_i_val, s_j_val, accumXvec);
	
		total -= ENTRIES_PER_INT;
	}

    // qLD matrix
    accumXvec[4] = qLD_res[s_i*alignment->segsites+s_j];
		
	assert(accumXvec[0]+accumXvec[1]==alignment->sequences);
	assert(accumXvec[2]+accumXvec[3]==alignment->sequences);	

	return computeCorrelationValueBINGPU(alignment->sequences, accumXvec);	
}

float computePairwiseCorrelationBINGAPSGPU (alignment_struct * alignment, int s_i, int s_j, uint32_t * qLD_res)
{
	int i, sequences=0;
	unsigned int accumXvec[5];
	unsigned int s_i_val, s_j_val, s_i_vld, s_j_vld;

	for(i=0;i<5;i++)
		accumXvec[i]=0;
	
	for(i=0;i<alignment->siteSize;i++)
	{
		s_i_val = alignment->compressedArrays[0][s_i*alignment->siteSize+i];
		s_j_val = alignment->compressedArrays[0][s_j*alignment->siteSize+i];

		s_i_vld = alignment->compressedArrays[1][s_i*alignment->siteSize+i];
		s_j_vld = alignment->compressedArrays[1][s_j*alignment->siteSize+i];

		sequences += count01GAPCombsGPU (s_i_val, s_j_val, s_i_vld, s_j_vld, accumXvec);		
	}

    // qLD matrix
    accumXvec[4] = qLD_res[s_i*alignment->segsites+s_j];
	
	assert(accumXvec[0]+accumXvec[1]==sequences);
	assert(accumXvec[2]+accumXvec[3]==sequences);	

	return computeCorrelationValueBINGPU(sequences, accumXvec);	
}

void computeCorrelationsBINGPU(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, uint32_t * qLD_res)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			s_i = i + omega[omegaIndex].leftIndex;
			s_j = j + omega[omegaIndex].leftIndex;

			alignment->correlationMatrix[i][j] = computePairwiseCorrelationBINGPU(alignment, s_i, s_j, qLD_res);
		}
	}	
}

void computeCorrelationsBINGAPSGPU(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, uint32_t * qLD_res)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			s_i = i + omega[omegaIndex].leftIndex;
			s_j = j + omega[omegaIndex].leftIndex;

			alignment->correlationMatrix[i][j] = computePairwiseCorrelationBINGAPSGPU(alignment, s_i, s_j, qLD_res);
		}
	}
}

void computeCorrelationMatrixPairwiseGPU(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRowIndex, void * threadData, cor_t ** myCorrelationMatrix, char * lookuptable, uint32_t * qLD_res)
{

	if (firstRowIndex==-1)
		return;

	switch(alignment->states)
	{
		case 2: computeCorrelationsBINGPU(alignment,omega,omegaIndex, firstRowIndex, qLD_res); 
			break;
		case 3: computeCorrelationsBINGAPSGPU(alignment,omega,omegaIndex, firstRowIndex, qLD_res); 
			break;
		case 4: printf("GPU LD calculation only works with binary data\n"); //computeCorrelationsDNA(alignment, omega, omegaIndex, firstRowIndex); 
			break;
		case 5: printf("GPU LD calculation only works with binary data\n"); //computeCorrelationsDNAGAPS(alignment, omega, omegaIndex, firstRowIndex); 
			break;
		default: assert(0);
	}	
}

/*   ---  GPU general functions  ---   */
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
    unsigned int i;

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

    // Let user choose best GPU
    char * vendor = NULL, * name = NULL;
    size_t valueSize;
    
    for(i = 0; i < num_devices; i++)
    {
        err |= clGetDeviceInfo(devices[i], CL_DEVICE_NAME, 0, NULL, &valueSize);
        name = (char*)malloc(valueSize);
        err |= clGetDeviceInfo(devices[i], CL_DEVICE_NAME, valueSize, name, NULL);
        err |= clGetDeviceInfo(devices[i], CL_DEVICE_VENDOR, 0, NULL, &valueSize);
        vendor = (char*)malloc(valueSize);
        err |= clGetDeviceInfo(devices[i], CL_DEVICE_VENDOR, valueSize, vendor, NULL);
        printCLErr(err,__LINE__,__FILE__);
        printf("GPU %u: %s - %s\n",i, vendor, name);
    }
    free(name);
    free(vendor);

    int gpu, result;
    do{
        printf("Which GPU do you want to use (enter # and press enter): ");
        result = scanf("%d", &gpu);
        if(result == 0)
            while (fgetc(stdin) != '\n'); // Read until a newline is found
    } while (result == EOF || result == 0 || gpu >= num_devices);


    context=clCreateContext(NULL, 1, &devices[gpu], NULL, NULL, &err);
    printCLErr(err,__LINE__,__FILE__);

    create_program_with_source(&program, &context, PROGRAM_FILE);

    // add any kernel compiler options to this string
    const char* options="-cl-mad-enable";
    err=clBuildProgram(program, 1, &devices[gpu], options, NULL, NULL);
    // print build errors
    if(err != CL_SUCCESS)
    {
        perror("error during build");
        size_t log_size=0;
        clGetProgramBuildInfo(program, devices[gpu], CL_PROGRAM_BUILD_LOG, 0, NULL,
                              &log_size);
        char *program_log=(char*)malloc(log_size+1);
        assert(program_log);
        program_log[log_size]='\0';
        clGetProgramBuildInfo(program, devices[gpu], CL_PROGRAM_BUILD_LOG,
                              log_size+1, program_log, NULL);
        printf("%s\n", program_log);
        free(program_log);
        exit(1);
    }

    io_queue=clCreateCommandQueue(context, devices[gpu], CL_QUEUE_PROFILING_ENABLE, &err);
    printCLErr(err,__LINE__,__FILE__);

    compute_queue=clCreateCommandQueue(context, devices[gpu], CL_QUEUE_PROFILING_ENABLE,
                                       &err);
    printCLErr(err,__LINE__,__FILE__);

	// create Omega kernel
	omega_kernel=clCreateKernel(program, OMEGA_NAME, &err);
	printCLErr(err,__LINE__,__FILE__);

	// create second Omega kernel
	omega_kernel2=clCreateKernel(program, OMEGA_NAME2, &err);
	printCLErr(err,__LINE__,__FILE__);

	// Get workgroup size or preferred size
	size_t max_group_size, pref_group_size;
	err = clGetKernelWorkGroupInfo(omega_kernel, devices[gpu],CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &max_group_size, NULL);
	// printf("Max kernel work-group size: %lu\n",max_group_size);	// 256
	err = clGetKernelWorkGroupInfo(omega_kernel, devices[gpu],CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &pref_group_size, NULL);
	// printf("Work-group pref. multiple: %lu\n",pref_group_size);	// 64

	// Get number of work groups / compute units
	err=clGetDeviceInfo(devices[gpu], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint),
                        &comp_units, NULL);
    printCLErr(err,__LINE__,__FILE__);

    printf("Compute units: %u\n", comp_units);

	group_size = pref_group_size;

	printf("Set work-group size: %lu\n", group_size);

	cl_ulong total;

    // NOTE: assume device has enough memory
    // TODO: does performance degrade if k is much less than KC?
    cl_ulong a_buffer_size=GPU_BLOCK_MC * GPU_BLOCK_KC * sizeof(inputDataType_x32);
    cl_ulong b_buffer_size=GPU_BLOCK_KC * GPU_BLOCK_NC * sizeof(inputDataType_x32);
    cl_ulong c_buffer_size=GPU_BLOCK_MC * GPU_BLOCK_NC * sizeof(inputDataType_x32);
	total = 2*a_buffer_size + 2*b_buffer_size + 2*c_buffer_size;

    cl_ulong max_alloc=0;
    err=clGetDeviceInfo(devices[gpu], CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(max_alloc),
                        &max_alloc, NULL);
    printCLErr(err,__LINE__,__FILE__);

    cl_ulong global_mem=0;
    err=clGetDeviceInfo(devices[gpu], CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(global_mem),
                        &global_mem, NULL);
    printCLErr(err,__LINE__,__FILE__);

    // cl_ulong omega_buffer_size 	= 512 * 40000 * sizeof(float);		// # work items is unknown here, see omega2
	// cl_ulong LRkm_buffer_size 	= 2 * GPU_BLOCK_MC * sizeof(float);
	// cl_ulong TS_buffer_size 	= omega_buffer_size;

    // Divide remaining memory correctly over needed buffers
    cl_ulong remain = global_mem - total;

    double omega_portion = 34.0128, LRkm_portion = 5314.5, TS_portion = 1.0629;

	cl_ulong omega_buffer_size 	= ((cl_ulong)(remain / omega_portion) / 256) * 256;		// # work items is unknown here, see omega2
	cl_ulong LRkm_buffer_size 	= ((cl_ulong)(remain / LRkm_portion) / 256) * 256;
	cl_ulong TS_buffer_size 	= ((cl_ulong)(remain / TS_portion) / 256) * 256;

	total += 2 * omega_buffer_size + 2 * LRkm_buffer_size + TS_buffer_size;

    if(total > global_mem)
    {
        printf("not enough global storage!\n");
        exit(1);
    }
    if((a_buffer_size > max_alloc) ||
       (b_buffer_size > max_alloc) ||
       (c_buffer_size > max_alloc) ||
	   (omega_buffer_size > max_alloc) ||
	   (LRkm_buffer_size > max_alloc) ||
	   TS_buffer_size > max_alloc)
    {
        printf("some buffer is too big!\n");
        exit(1);
    }

    max_omegas = omega_buffer_size / sizeof(float);
    max_LRkm = LRkm_buffer_size / sizeof(float);
    max_TS = TS_buffer_size / sizeof(float);

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

	omega_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, omega_buffer_size, NULL, &err);
	printCLErr(err,__LINE__,__FILE__);
	index_buffer=clCreateBuffer(context, CL_MEM_WRITE_ONLY, omega_buffer_size, NULL, &err);
	printCLErr(err,__LINE__,__FILE__);
	LR_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, LRkm_buffer_size, NULL, &err);
	printCLErr(err,__LINE__,__FILE__);
	TS_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, TS_buffer_size, NULL, &err);
	printCLErr(err,__LINE__,__FILE__);
	km_buffer=clCreateBuffer(context, CL_MEM_READ_ONLY, LRkm_buffer_size, NULL, &err);
	printCLErr(err,__LINE__,__FILE__);

	// set kernel arguements for buffers
	err |= clSetKernelArg(omega_kernel, 0, sizeof(cl_mem), &omega_buffer);
	err |= clSetKernelArg(omega_kernel, 1, sizeof(cl_mem), &index_buffer);
	err |= clSetKernelArg(omega_kernel, 2, sizeof(cl_mem), &LR_buffer);
	err |= clSetKernelArg(omega_kernel, 3, sizeof(cl_mem), &TS_buffer);
	err |= clSetKernelArg(omega_kernel, 4, sizeof(cl_mem), &km_buffer);
	printCLErr(err,__LINE__,__FILE__);
	
	// set second kernel arguements for buffers
	err |= clSetKernelArg(omega_kernel2, 0, sizeof(cl_mem), &omega_buffer);
	err |= clSetKernelArg(omega_kernel2, 1, sizeof(cl_mem), &LR_buffer);
	err |= clSetKernelArg(omega_kernel2, 2, sizeof(cl_mem), &TS_buffer);
	err |= clSetKernelArg(omega_kernel2, 3, sizeof(cl_mem), &km_buffer);
	printCLErr(err,__LINE__,__FILE__);
	

	// Get device max workgroup size
	// size_t group_size;
	// err=clGetDeviceInfo(devices[gpu], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t),
    //                     &group_size, NULL);
    // printCLErr(err,__LINE__,__FILE__);
	// printf("max workgroup size: %lu\n",group_size);
	// Get constant memory info
	// cl_uint max_args;
	// cl_ulong max_const;
	// err=clGetDeviceInfo(devices[gpu], CL_DEVICE_MAX_CONSTANT_ARGS, sizeof(cl_uint),
    //                     &max_args, NULL);		// 16
    // printCLErr(err,__LINE__,__FILE__);
	// err=clGetDeviceInfo(devices[gpu], CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(cl_ulong),
    //                     &max_const, NULL);		// 1503238400
    // printCLErr(err,__LINE__,__FILE__);
	// printf("Arg: %u, Size: %lu\n",max_args,max_const);
	// Get local memory info
	// cl_ulong max_local;
	// err=clGetDeviceInfo(devices[gpu], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong),
    //                     &max_local, NULL);		// 32768
    // printCLErr(err,__LINE__,__FILE__);
	// printf("Size: %lu\n",max_local);
	// Get timer resolution info
	// cl_ulong resolution;
	// clGetDeviceInfo(devices[gpu], CL_DEVICE_PROFILING_TIMER_RESOLUTION, sizeof(cl_ulong), 
	// 						&resolution, NULL);		// 0
	// printCLErr(err,__LINE__,__FILE__);
	// printf("Res: %lu\n",resolution);
	// float * test = NULL;
}

void gpu_release(void)
{
    // Omega
	clReleaseKernel(omega_kernel);
    clReleaseKernel(omega_kernel2);
	clReleaseMemObject(omega_buffer);
	clReleaseMemObject(TS_buffer);
	clReleaseMemObject(LR_buffer);
	clReleaseMemObject(km_buffer);
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
    clReleaseCommandQueue(io_queue);
    clReleaseCommandQueue(compute_queue);
    clReleaseProgram(program);
    clReleaseContext(context);
    free(devices);
    free(platforms);
}