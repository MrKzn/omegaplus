
#include "OmegaPlus.h"

int compute_genLoad(int startIndex, int finishIndex, int prev_totalLoad, int prev_startIndex, int prev_finishIndex)
{

	if((prev_startIndex==startIndex)&&(prev_finishIndex==finishIndex))
	{
		return prev_totalLoad;
	}


	int counter2;
	int a = finishIndex+1-startIndex;
	int b = prev_finishIndex+1-startIndex;	

	counter2 = (a*a+a)*0.5;

	if(startIndex<=prev_finishIndex)
	{
		counter2 = counter2 - ((b*b+b)*0.5);
	}

	return counter2 + prev_totalLoad;

}

inline int get_SNP_groupID(int SNP_index)
{
	return SNP_index / SNP_GROUP_SIZE;
}

inline float get_Mem_GroupSize()
{
	return SNP_GROUP_SIZE_POW2 * 4.0 / 1024.0 / 1024.0;
}

void update_workgroup_map_ptr(float *** workgroup_map_ptr, int start, int finish, int first_group_index)
{
	int i,j;
	
	for(i=start;i<=finish;i++)
	{
		for(j=start;j<=i;j++)
		{	
			if(workgroup_map_ptr[i-first_group_index][j-first_group_index]==NULL)
			{
				workgroup_map_ptr[i-first_group_index][j-first_group_index] = malloc(sizeof(float)*SNP_GROUP_SIZE_POW2);
			}
		}
	}

}

void update_workgroup_map_partial_ptr(float *** workgroup_map_ptr, int start, int finish, int prev_start, int prev_finish, int first_group_index)
{

	if(prev_start==start)
		if(prev_finish==finish)
			return;

	int i,j;
	if(prev_finish>=start)
	{
		for(i=prev_finish+1;i<=finish;i++)
		{
			for(j=start;j<=i;j++)
			{
				if(workgroup_map_ptr[i-first_group_index][j-first_group_index]==NULL)
				{
					workgroup_map_ptr[i-first_group_index][j-first_group_index] = malloc(sizeof(float)*SNP_GROUP_SIZE_POW2);
				}
			}
		}
		
	}
	else
	{
		for(i=start;i<=finish;i++)
		{
			for(j=start;j<=i;j++)
			{
				if(workgroup_map_ptr[i-first_group_index][j-first_group_index]==NULL)
				{
					workgroup_map_ptr[i-first_group_index][j-first_group_index] = malloc(sizeof(float)*SNP_GROUP_SIZE_POW2);
				}
			}
		}
	}
}

void count01Combs_light (int total, unsigned int inputL, unsigned int inputR, unsigned int * accumXvec)
{
//	int limit=ENTRIES_PER_INT;
	
//	if (total<limit)
//		limit = total;

	int t1 = precomputed16_bitcount (inputL);
	int t2 = precomputed16_bitcount (inputR);

	unsigned int combAND = inputL & inputR;

	accumXvec[1] += t1;  
	accumXvec[3] += t2;

//	accumXvec[0] += limit - t1;
//	accumXvec[2] += limit - t2;

	accumXvec[4] +=  precomputed16_bitcount (combAND);	
}

cor_t computePairwiseCorrelationBIN_with_border_check (alignment_struct * alignment, int s_i, int s_j, int leftSNPindex, int rightSNPindex, uint32_t * qLD_res)
{

	if(s_i<leftSNPindex)
		return 0.0;
	
	if(s_j<leftSNPindex)
		return 0.0;

	if(s_i>rightSNPindex)
		return 0.0;

	if(s_j>rightSNPindex)
		return 0.0;


	int i, total;
	unsigned int accumXvec[5];
	unsigned int s_i_val, s_j_val;

//	total = alignment->sequences;

	for(i=0;i<5;i++)
		accumXvec[i]=0;
	
	for(i=0;i<alignment->siteSize;i++)
	{
		s_i_val = alignment->compressedArrays[0][s_i*alignment->siteSize+i];
		s_j_val = alignment->compressedArrays[0][s_j*alignment->siteSize+i];

		count01Combs_light (total, s_i_val, s_j_val, accumXvec);
	
//		total -= ENTRIES_PER_INT;
	}

	// accumXvec[4] = qLD_res[s_i*alignment->segsites+s_j];

	accumXvec[0]=alignment->sequences-accumXvec[1];
	accumXvec[2]=alignment->sequences-accumXvec[3];

	if(accumXvec[0]==0 || accumXvec[1]==0 || accumXvec[2]==0 || accumXvec[3]==0)
		return 0.0;
		
//	assert(accumXvec[0]+accumXvec[1]==alignment->sequences);
//	assert(accumXvec[2]+accumXvec[3]==alignment->sequences);	

	return computeCorrelationValueBIN(alignment->sequences, accumXvec);
		
}

cor_t computePairwiseCorrelationBINGAPS_with_border_check (alignment_struct * alignment, int s_i, int s_j, int leftSNPindex, int rightSNPindex)
{

	if(s_i<leftSNPindex)
		return 0.0;
	
	if(s_j<leftSNPindex)
		return 0.0;

	if(s_i>rightSNPindex)
		return 0.0;

	if(s_j>rightSNPindex)
		return 0.0;



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

		sequences += count01GAPCombs (s_i_val, s_j_val, s_i_vld, s_j_vld, accumXvec);		
	}
	
	assert(accumXvec[0]+accumXvec[1]==sequences);
	assert(accumXvec[2]+accumXvec[3]==sequences);	

	return computeCorrelationValueBIN(sequences, accumXvec);	
}


cor_t computePairwiseCorrelationDNA_with_border_check (alignment_struct * alignment, int s_i, int s_j, int leftSNPindex, int rightSNPindex)
{

	if(s_i<leftSNPindex)
		return 0.0;
	
	if(s_j<leftSNPindex)
		return 0.0;

	if(s_i>rightSNPindex)
		return 0.0;

	if(s_j>rightSNPindex)
		return 0.0;


	int i,j,k, total=0;
	unsigned int accumXvec[5],
                     validVec[8];

	cor_t pairwiseCorrelationMatrix [4][4]; // A->0 C->1 G->2 T->3

	unsigned int s_i_val, s_j_val;
   	
	for(i=0;i<8;i++)
		validVec[i]=0;
	
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			for(k=0;k<5;k++)
				accumXvec[k]=0;

			total = alignment->sequences;

			for(k=0;k<alignment->siteSize;k++)
			{
				s_i_val = alignment->compressedArrays[i][s_i*alignment->siteSize+k];
				s_j_val = alignment->compressedArrays[j][s_j*alignment->siteSize+k];
				
				count01Combs (total, s_i_val, s_j_val, accumXvec);

				total -= ENTRIES_PER_INT;
			}
			
			if(accumXvec[1]==0)
				break;

			validVec[i] = 1;			
			
			if(accumXvec[3]==0)
				continue;

			validVec[j+4] = 1;

			pairwiseCorrelationMatrix[i][j]=computeCorrelationValueBIN(alignment->sequences,accumXvec);
		}
	}
	return computeCorrelationValueDNA (alignment->sequences, pairwiseCorrelationMatrix, validVec);	
}


cor_t computePairwiseCorrelationDNAGAPS_with_border_check (alignment_struct * alignment, int s_i, int s_j, int leftSNPindex, int rightSNPindex)
{

	if(s_i<leftSNPindex)
		return 0.0;
	
	if(s_j<leftSNPindex)
		return 0.0;

	if(s_i>rightSNPindex)
		return 0.0;

	if(s_j>rightSNPindex)
		return 0.0;


	int i,j,k,sequences;
	
	unsigned int    accumXvec[5],
                     validVec[8];

	float pairwiseCorrelationMatrix [4][4]; // A->0 C->1 G->2 T->3

	unsigned int s_i_val, s_j_val, s_i_vld, s_j_vld;        

        for(i=0;i<8;i++)
		validVec[i]=0;	
	
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			sequences=0;

			for(k=0;k<5;k++)
				accumXvec[k]=0;

			for(k=0;k<alignment->siteSize;k++)
			{
				s_i_val = alignment->compressedArrays[i][s_i*alignment->siteSize+k];
				s_j_val = alignment->compressedArrays[j][s_j*alignment->siteSize+k];

				s_i_vld = alignment->compressedArrays[4][s_i*alignment->siteSize+k];
				s_j_vld = alignment->compressedArrays[4][s_j*alignment->siteSize+k];

				sequences += count01GAPCombs (s_i_val, s_j_val, s_i_vld, s_j_vld, accumXvec);
			}
				
                        if(accumXvec[1]==0)
				break;

			validVec[i] = 1;
			
			if(accumXvec[3]==0)
				continue;

			validVec[j+4] = 1;

			pairwiseCorrelationMatrix[i][j]=computeCorrelationValueBIN(sequences,accumXvec);
		}
	}
	return computeCorrelationValueDNA (sequences, pairwiseCorrelationMatrix, validVec);

		
}

cor_t computePairwiseCorrelationGENERIC_with_border_check (alignment_struct * alignment, int s_i, int s_j, int leftSNPindex, int rightSNPindex, uint32_t* qLD_res)
{

	switch(alignment->states)
	{	
		case 2: return computePairwiseCorrelationBIN_with_border_check(alignment, s_i, s_j, leftSNPindex, rightSNPindex, qLD_res); 
			break;
		case 3: return computePairwiseCorrelationBINGAPS_with_border_check(alignment, s_i, s_j, leftSNPindex, rightSNPindex); 
			break;
		case 4: return computePairwiseCorrelationDNA_with_border_check(alignment, s_i, s_j, leftSNPindex, rightSNPindex);  
			break;
		case 5: return computePairwiseCorrelationDNAGAPS_with_border_check(alignment, s_i, s_j, leftSNPindex, rightSNPindex);  
			break;
		default: assert(0);
	}
	
}

void apply_dp_on_diag_triangles_ptr (float * gen_output_vector_group_ptr, int t_x_SNP_GROUP_SIZE, int t_y_SNP_GROUP_SIZE, alignment_struct * alignment, int leftSNPindex, int rightSNPindex, uint32_t * qLD_res)
{
	int m=0;
	int l=0;
	int m_SNP_GROUP_SIZE;

	//Pairwise correlation calculation
	for(m=0;m<SNP_GROUP_SIZE;m++) // m is row indexing 
		gen_output_vector_group_ptr[m*SNP_GROUP_SIZE+m] = 0.0;
	
	for(m=1;m<SNP_GROUP_SIZE;m++) // m is row indexing
	{
		m_SNP_GROUP_SIZE = m * SNP_GROUP_SIZE;
		for(l=m-1;l>=0;l--) // l is col indexing
			gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+l] = computePairwiseCorrelationGENERIC_with_border_check(alignment, t_x_SNP_GROUP_SIZE+m, t_y_SNP_GROUP_SIZE+l, leftSNPindex, rightSNPindex, qLD_res);
	}


	//DP on lower triangular part of tile
	for(m=2;m<SNP_GROUP_SIZE;m++) // m is row indexing
	{
		m_SNP_GROUP_SIZE = m * SNP_GROUP_SIZE;
		for(l=m-2;l>=0;l--) // l is col indexing
			gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+l] = gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+l] + 
									  gen_output_vector_group_ptr[m_SNP_GROUP_SIZE-SNP_GROUP_SIZE+l] + 										  gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+(l+1)] -
									  gen_output_vector_group_ptr[m_SNP_GROUP_SIZE-SNP_GROUP_SIZE+(l+1)];
	}	
}

void apply_dp_on_tiles_diag_one_down_ptr (float * gen_output_vector_group_ptr, float * gen_output_vector_group_above_ptr, float * gen_output_vector_group_right_ptr, int t_x_SNP_GROUP_SIZE, int t_y_SNP_GROUP_SIZE, alignment_struct * alignment, int leftSNPindex, int rightSNPindex, uint32_t * qLD_res)
{
	int m=0;
	int l=0;
	int m_SNP_GROUP_SIZE;

	//Pairwise correlation calculation
	for(m=0;m<SNP_GROUP_SIZE;m++) // m is row indexing
	{
		m_SNP_GROUP_SIZE = m * SNP_GROUP_SIZE;
		for(l=SNP_GROUP_SIZE-1;l>=0;l--) // l is col indexing
			gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+l] = computePairwiseCorrelationGENERIC_with_border_check(alignment, t_y_SNP_GROUP_SIZE+m, t_x_SNP_GROUP_SIZE+l, leftSNPindex, rightSNPindex, qLD_res);
	}

	//DP on tile
	//First row
	for(l=SNP_GROUP_SIZE-2;l>=0;l--)
		gen_output_vector_group_ptr[l] = gen_output_vector_group_ptr[l] + 
						 gen_output_vector_group_above_ptr[(SNP_GROUP_SIZE-1)*SNP_GROUP_SIZE+l] +
						 gen_output_vector_group_ptr[l+1] - 								
						 gen_output_vector_group_above_ptr[(SNP_GROUP_SIZE-1)*SNP_GROUP_SIZE+l+1];

	//Last column
	l=SNP_GROUP_SIZE-1;
	for(m=1;m<SNP_GROUP_SIZE;m++)
		gen_output_vector_group_ptr[m*SNP_GROUP_SIZE+l] = gen_output_vector_group_ptr[m*SNP_GROUP_SIZE+l] + 
								  gen_output_vector_group_ptr[(m-1)*SNP_GROUP_SIZE+l] +
						 		  gen_output_vector_group_right_ptr[m*SNP_GROUP_SIZE+(0)] - 
								  gen_output_vector_group_right_ptr[(m-1)*SNP_GROUP_SIZE+0];

	//Rest of tile
	for(m=1;m<SNP_GROUP_SIZE;m++) // m is row indexing
	{
		m_SNP_GROUP_SIZE = m * SNP_GROUP_SIZE;
		for(l=SNP_GROUP_SIZE-2;l>=0;l--) // l is col indexing
			gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+l] =  gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+l] + 										
									   gen_output_vector_group_ptr[m_SNP_GROUP_SIZE-SNP_GROUP_SIZE+l] + 										   gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+(l+1)] - 										   gen_output_vector_group_ptr[m_SNP_GROUP_SIZE-SNP_GROUP_SIZE+(l+1)];
	}	
}

void apply_dp_on_tiles_all_diags_ptr (float * gen_output_vector_group_ptr, float * gen_output_vector_group_above_ptr, float * gen_output_vector_group_right_ptr, float * gen_output_vector_group_diag_ptr, int t_x_SNP_GROUP_SIZE, int t_y_SNP_GROUP_SIZE, alignment_struct * alignment, int leftSNPindex, int rightSNPindex, uint32_t * qLD_res)
{
	int m=0; 
	int l=0;
	int m_SNP_GROUP_SIZE;

	//Pairwise correlation calculation
	for(m=0;m<SNP_GROUP_SIZE;m++) // m is row indexing
	{
		m_SNP_GROUP_SIZE = m * SNP_GROUP_SIZE;
		for(l=SNP_GROUP_SIZE-1;l>=0;l--) // l is col indexing
			gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+l] = computePairwiseCorrelationGENERIC_with_border_check(alignment, t_y_SNP_GROUP_SIZE+m, t_x_SNP_GROUP_SIZE+l, leftSNPindex, rightSNPindex, qLD_res);
	}

	//DP on tile
	//First row last column
	l=SNP_GROUP_SIZE-1;
	gen_output_vector_group_ptr[l] = gen_output_vector_group_ptr[l] + 
					 gen_output_vector_group_above_ptr[(SNP_GROUP_SIZE-1)*SNP_GROUP_SIZE+l] +
					 gen_output_vector_group_right_ptr[0] - 
					 gen_output_vector_group_diag_ptr[(SNP_GROUP_SIZE-1)*SNP_GROUP_SIZE+0] ;

	//Rest first row
	for(l=SNP_GROUP_SIZE-2;l>=0;l--)
		gen_output_vector_group_ptr[l] = gen_output_vector_group_ptr[l] + 
						     gen_output_vector_group_above_ptr[(SNP_GROUP_SIZE-1)*SNP_GROUP_SIZE+l] +
								    gen_output_vector_group_ptr[(l+1)] - 						
								     gen_output_vector_group_above_ptr[(SNP_GROUP_SIZE-1)*SNP_GROUP_SIZE+l+1];

	//Last column
	l=SNP_GROUP_SIZE-1;
	for(m=1;m<SNP_GROUP_SIZE;m++)
		gen_output_vector_group_ptr[m*SNP_GROUP_SIZE+l] = gen_output_vector_group_ptr[m*SNP_GROUP_SIZE+l] +
										      gen_output_vector_group_ptr[(m-1)*SNP_GROUP_SIZE+l] + 
										      gen_output_vector_group_right_ptr[m*SNP_GROUP_SIZE+(0)] - 
										      gen_output_vector_group_right_ptr[(m-1)*SNP_GROUP_SIZE+0] ;

	//Rest of tile
	for(m=1;m<SNP_GROUP_SIZE;m++) // m is row indexing
	{
		m_SNP_GROUP_SIZE = m * SNP_GROUP_SIZE;
		for(l=SNP_GROUP_SIZE-2;l>=0;l--) // l is col indexing
			gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+l] = gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+l] + 												      gen_output_vector_group_ptr[m_SNP_GROUP_SIZE-SNP_GROUP_SIZE+l] + 												      gen_output_vector_group_ptr[m_SNP_GROUP_SIZE+(l+1)] - 											      		gen_output_vector_group_ptr[m_SNP_GROUP_SIZE-SNP_GROUP_SIZE+(l+1)];
	}
}


void dp_on_tiles_overlap_ptr (int first_DP_tile, int last_DP_tile, float *** workgroup_map_ptr, float *** overlap_workgroup_map_ptr, int overlap, int first_group_index, alignment_struct * alignment, int leftSNPindex, int rightSNPindex, uint32_t * qLD_res)	
{

	int i, j, t_x;
	int diag_total = last_DP_tile - first_DP_tile + 1; 

	int diag_i_first_tile_y = first_DP_tile, 
	    diag_i_last_tile_x = last_DP_tile;


	#pragma omp parallel for private (t_x)
	for(t_x=first_DP_tile;t_x<=diag_i_last_tile_x;t_x++) // t_x is indexing x dimension
	{
		int t_y = t_x-first_DP_tile+diag_i_first_tile_y;
		int t_x_SNP_GROUP_SIZE=t_x*SNP_GROUP_SIZE;
		int t_y_SNP_GROUP_SIZE=t_y*SNP_GROUP_SIZE;
		int skip=0;
		if(t_y<first_DP_tile+overlap-1)
		{
			if(overlap_workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile]!=NULL)
			{	
				skip=1;
			}
		}	
		
		if(workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile]!=NULL)
		{
			if(skip)
			{
				free(workgroup_map_ptr[t_y-first_group_index][t_x-first_group_index]);
				workgroup_map_ptr[t_y-first_group_index][t_x-first_group_index] = overlap_workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile];
			}
			else
			{
				apply_dp_on_diag_triangles_ptr (workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile], t_x_SNP_GROUP_SIZE, t_y_SNP_GROUP_SIZE, alignment, leftSNPindex, rightSNPindex, qLD_res);
			}
		}	
	}

	diag_i_first_tile_y++;
	diag_i_last_tile_x--;

	// One diagonal below the main diagonal: special tile processing	
	#pragma omp parallel for private (t_x)		
	for(t_x=first_DP_tile;t_x<=diag_i_last_tile_x;t_x++) // t_x is indexing x dimension
	{
		int t_y = t_x-first_DP_tile+diag_i_first_tile_y;
		int t_x_SNP_GROUP_SIZE=t_x*SNP_GROUP_SIZE;
		int t_y_SNP_GROUP_SIZE=t_y*SNP_GROUP_SIZE;
		int skip=0;

		if(t_y<first_DP_tile+overlap-1)
		{
			if(overlap_workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile]!=NULL)
			{	
				skip=1;
			}
		}

		
		if(workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile]!=NULL) 
		{
			if(skip)
			{

				free(workgroup_map_ptr[t_y-first_group_index][t_x-first_group_index]);
				workgroup_map_ptr[t_y-first_group_index][t_x-first_group_index] = overlap_workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile];
			}
			else
			{
				apply_dp_on_tiles_diag_one_down_ptr (workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile],
									workgroup_map_ptr[t_y-1-first_DP_tile][t_x-first_DP_tile],
									workgroup_map_ptr[t_y-first_DP_tile][t_x+1-first_DP_tile],
									t_x_SNP_GROUP_SIZE, t_y_SNP_GROUP_SIZE,  alignment, leftSNPindex, rightSNPindex, qLD_res);
			}
		}
	}

	diag_i_first_tile_y++;
	diag_i_last_tile_x--;

	int early_stop;
	int * early_stop_vec=calloc(omp_get_max_threads(),sizeof(int));

	// All remaining diagonals: standard tile processing
	for(i=2;i<diag_total;i++) // i is indexing the diagonals
	{	
		early_stop=0;
		#pragma omp parallel for private (t_x)		
		for(t_x=first_DP_tile;t_x<=diag_i_last_tile_x;t_x++) // t_x is indexing x dimension
		{
			int t_y = (t_x-first_DP_tile)+diag_i_first_tile_y;
			int t_x_SNP_GROUP_SIZE=t_x*SNP_GROUP_SIZE;
			int t_y_SNP_GROUP_SIZE=t_y*SNP_GROUP_SIZE;		
			int skip=0;

			if(t_y<first_DP_tile+overlap-1)
			{
				if(overlap_workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile]!=NULL)
				{	
					skip=1;
				}
			}

			if(workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile]!=NULL)
			{	

				if(skip)
				{
					free(workgroup_map_ptr[t_y-first_group_index][t_x-first_group_index]);
					workgroup_map_ptr[t_y-first_group_index][t_x-first_group_index] = overlap_workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile];
				}
				else
				{
					apply_dp_on_tiles_all_diags_ptr (workgroup_map_ptr[t_y-first_DP_tile][t_x-first_DP_tile],
									workgroup_map_ptr[t_y-1-first_DP_tile][t_x-first_DP_tile],
									workgroup_map_ptr[t_y-first_DP_tile][t_x+1-first_DP_tile],
									workgroup_map_ptr[t_y-1-first_DP_tile][t_x+1-first_DP_tile],
									t_x_SNP_GROUP_SIZE, t_y_SNP_GROUP_SIZE,  alignment, leftSNPindex, rightSNPindex, qLD_res);
				}
			}
			else
			{
				early_stop_vec[omp_get_thread_num()]++;
			}
		
		}

		for(j=0;j<omp_get_max_threads();j++)
			early_stop+=early_stop_vec[j];

		if(early_stop==diag_i_last_tile_x-first_DP_tile+1)
			i=diag_total;

		diag_i_first_tile_y++;
		diag_i_last_tile_x--;					
	}

	free(early_stop_vec);
}


void computeOmegaValues_generic (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData, float *** workgroup_map_ptr, int first_group_index)
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

	int omegaSNIPindex_real = omega[omegaIndex].omegaPos;

#ifdef _UNROLL
	int vw = 2;
	int iter, iterations = (rightMaxIndex-rightMinIndex+1) / vw;
	int finaliterations = (rightMaxIndex-rightMinIndex+1) % vw;
	float maxW0=0.0, maxW1=0.0;
	int maxLeftIndex0=0, maxRightIndex0=0, maxLeftIndex1=0, maxRightIndex1=0;
	int indexj01 = (omegaSNIPindex_real+1) / SNP_GROUP_SIZE; 
	int indexj01p = indexj01 - first_group_index;
	int indexj02pp = (omegaSNIPindex_real+1)%SNP_GROUP_SIZE;
#endif
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{
		LS=workgroup_map_ptr[get_SNP_groupID(omegaSNIPindex_real)-first_group_index][get_SNP_groupID(i+omega[omegaIndex].leftIndex)-first_group_index][(omegaSNIPindex_real%SNP_GROUP_SIZE)*SNP_GROUP_SIZE+(i+omega[omegaIndex].leftIndex)%SNP_GROUP_SIZE];

		k = omegaSNIPIndex - i + 1;
		
		ksel2 = (k * (k-1)) * 0.5;

		if(borderTol > 0)
		{
			rightMinIndex = rightMinIndexORIG;

			rightMaxIndex = rightMaxIndexORIG;

			int leftSNPs = omegaSNIPIndex - i + 1;
			int equalRightPosition = omegaSNIPIndex + leftSNPs;

			rightMinIndex = max(rightMinIndex, equalRightPosition - borderTol);
			rightMaxIndex = min(rightMaxIndex, equalRightPosition + borderTol);
		    
		}
#ifdef _UNROLL
		int indexj03 = (i+omega[omegaIndex].leftIndex) / SNP_GROUP_SIZE;
		int indexj13 = indexj03;//(i+omega[omegaIndex].leftIndex) / SNP_GROUP_SIZE;

		int indexj04 = (i+omega[omegaIndex].leftIndex)%SNP_GROUP_SIZE;
		int indexj14 = indexj04;//(i+omega[omegaIndex].leftIndex)%SNP_GROUP_SIZE;
	
		j = rightMinIndex;

		for(iter=0;iter<iterations;iter++)
		{
			int j0 = j;
			int j1 = j0 + 1;

			j=j+vw;

			int indexj00 = (j0+omega[omegaIndex].leftIndex) / SNP_GROUP_SIZE;
			int indexj10 = (j1+omega[omegaIndex].leftIndex) / SNP_GROUP_SIZE;

			int indexj00p = indexj00 - first_group_index;
			int indexj10p = indexj10 - first_group_index;


			int indexj02 = (j0+omega[omegaIndex].leftIndex)%SNP_GROUP_SIZE;
			int indexj12 = (j1+omega[omegaIndex].leftIndex)%SNP_GROUP_SIZE;

			int indexj02p = indexj02*SNP_GROUP_SIZE;
			int indexj12p = indexj12*SNP_GROUP_SIZE;
 
			int indexj02ppp = indexj02p +indexj02pp;			 
			int indexj12ppp = indexj12p +indexj02pp;

			float RS0 = workgroup_map_ptr[indexj00p][indexj01p][indexj02ppp];
			float RS1 = workgroup_map_ptr[indexj10p][indexj01p][indexj12ppp];

			int m0 = j0 - omegaSNIPIndex;
			int m1 = j1 - omegaSNIPIndex;

			int msel20 = (m0 * (m0-1)) * 0.5;
			int msel21 = (m1 * (m1-1)) * 0.5;
			
			int indexj03p = indexj03 - first_group_index;
			int indexj13p = indexj13 - first_group_index;
			
			int indexj04p = indexj04 + indexj02p;
			int indexj14p = indexj14 + indexj12p;

			float TS0 = workgroup_map_ptr[indexj00p][indexj03p][indexj04p];
			float TS1 = workgroup_map_ptr[indexj10p][indexj13p][indexj14p];


			float numerator0 = (LS + RS0) / (ksel2 + msel20);
			float numerator1 = (LS + RS1) / (ksel2 + msel21);


			float denominator0 = (TS0 - LS - RS0) / (k*m0) + DENOMINATOR_OFFSET;
			float denominator1 = (TS1 - LS - RS1) / (k*m1) + DENOMINATOR_OFFSET;


			float tmpW0 = numerator0 / denominator0;
			float tmpW1 = numerator1 / denominator1;

			if(tmpW0>maxW0)
			{
				maxW0 = tmpW0;
				maxLeftIndex0 = i + omega[omegaIndex].leftIndex;
				maxRightIndex0 = j0 + omega[omegaIndex].leftIndex;
			}

			if(tmpW1>maxW1)
			{
				maxW1 = tmpW1;
				maxLeftIndex1 = i + omega[omegaIndex].leftIndex;
				maxRightIndex1 = j1 + omega[omegaIndex].leftIndex;
			}
			
		}

		if(maxW0>maxW1)
		{
			maxW = maxW0;
			maxLeftIndex = maxLeftIndex0;
			maxRightIndex = maxRightIndex0;
		}
		else
		{
			maxW = maxW1;
			maxLeftIndex = maxLeftIndex1;
			maxRightIndex = maxRightIndex1;
		}



		if(finaliterations!=0)
			for(iter=0;iter<finaliterations;iter++)
			{
				RS = workgroup_map_ptr[ get_SNP_groupID(j+omega[omegaIndex].leftIndex)-first_group_index] [get_SNP_groupID(omegaSNIPindex_real+1)-first_group_index ][((j+omega[omegaIndex].leftIndex)%SNP_GROUP_SIZE)*SNP_GROUP_SIZE+(omegaSNIPindex_real+1)%SNP_GROUP_SIZE];

				m = j - omegaSNIPIndex;

				msel2 = (m * (m-1)) * 0.5;

				TS = workgroup_map_ptr[ get_SNP_groupID(j+omega[omegaIndex].leftIndex)-first_group_index] [get_SNP_groupID(i+omega[omegaIndex].leftIndex)-first_group_index ][((j+omega[omegaIndex].leftIndex)%SNP_GROUP_SIZE)*SNP_GROUP_SIZE+(i+omega[omegaIndex].leftIndex)%SNP_GROUP_SIZE];

				tmpW = computeOmega(LS, RS, TS, k, ksel2, m, msel2);

				if(tmpW>maxW)
				{
					maxW = tmpW;
					maxLeftIndex = i + omega[omegaIndex].leftIndex;
					maxRightIndex = j + omega[omegaIndex].leftIndex;
				}

			j++;		
			}


		maxW0 = maxW;
		maxLeftIndex0 = maxLeftIndex;
		maxRightIndex0 = maxRightIndex;
		
		maxW1 = maxW;
		maxLeftIndex1 = maxLeftIndex;
		maxRightIndex1 = maxRightIndex;
#else		
		
		for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
		{

			RS = workgroup_map_ptr[ get_SNP_groupID(j+omega[omegaIndex].leftIndex)-first_group_index] [get_SNP_groupID(omegaSNIPindex_real+1)-first_group_index ][((j+omega[omegaIndex].leftIndex)%SNP_GROUP_SIZE)*SNP_GROUP_SIZE+(omegaSNIPindex_real+1)%SNP_GROUP_SIZE];

			m = j - omegaSNIPIndex;

			msel2 = (m * (m-1)) * 0.5;

			TS = workgroup_map_ptr[ get_SNP_groupID(j+omega[omegaIndex].leftIndex)-first_group_index] [get_SNP_groupID(i+omega[omegaIndex].leftIndex)-first_group_index ][((j+omega[omegaIndex].leftIndex)%SNP_GROUP_SIZE)*SNP_GROUP_SIZE+(i+omega[omegaIndex].leftIndex)%SNP_GROUP_SIZE];

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

	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;	

}


void computeOmegas_generic (omega_struct * omega, int omegaIndex, float *** workgroup_map_ptr,  int first_group_index)
{
	computeOmegaValues_generic (omega, omegaIndex, NULL, NULL, workgroup_map_ptr,  first_group_index);
}

// qLD ADDED
void Pack_A(uint32_t *A,
            unsigned int lda,
            uint32_t *A_pack,
            unsigned int m,
            unsigned int k)
{
	uint32_t *A_pack_local;
	for(unsigned int ic=0;ic<m;ic+=BLOCK_MR)
    {
		A_pack_local=&A_pack[ic*k];
		unsigned int m_alg=fmin(BLOCK_MR,m-ic);
		for(unsigned int pc=0;pc<k;pc++)
        {
			for(unsigned int ir=0;ir<m_alg;ir++)
            {
				A_pack_local[0]=A[(ic+ir)+pc*lda];
				A_pack_local++;
			}
		}
	}
}

void Pack_B(uint32_t *B,
            unsigned int ldb,
            uint32_t *B_pack,
            unsigned int k,
            unsigned int n)
{
    uint32_t *B_pack_local;
	for(unsigned int jc=0;jc<n;jc+=BLOCK_NR)
    {
	    B_pack_local=&B_pack[jc*k];
        unsigned int n_alg=fmin(BLOCK_NR,n-jc);
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

void dgemm_ref(int k,
               int mr_alg,
               int nr_alg,
               uint32_t alpha,
               uint32_t* restrict a,
               uint32_t* restrict b,
               uint32_t* restrict c,
               int rs_c,
               int cs_c)
{
	uint32_t ab[mr_alg*nr_alg];
	uint32_t bj, ai;
	unsigned int l,i,j;
	for(i=0;i < (unsigned int)(mr_alg * nr_alg); ++i) //set 0s
        *(ab+i)=0;

	for(l=0;l < (unsigned int)k; ++l)
    {
		uint32_t *abij=ab;
        for(j=0;j < (unsigned int)nr_alg; ++j)
        {
            bj = *(b + j);
            for(i=0;i < (unsigned int)mr_alg; ++i)
            {
                ai = *(a + i);
                if(sizeof(uint32_t) <= 4)
                    abij[0] += __builtin_popcount(ai & bj);
                else if(sizeof(uint32_t) < 8)
                    abij[0] += __builtin_popcountl(ai & bj);
                else if(sizeof(uint32_t) >= 8)
                    abij[0] += __builtin_popcountll(ai & bj);
                abij++;
            }
        }
        a+=mr_alg;
        b+=nr_alg;
    }

	for(i=0;i < (unsigned int)(mr_alg * nr_alg); ++i) //Scale by alpha
        ab[i]*=alpha;

	uint32_t *Cr_l;
	uint32_t *ab_l=ab;

	for(j=0;j < (unsigned int)nr_alg; ++j)
    {
        Cr_l=&c[j*cs_c];
        for(i=0;i < (unsigned int)mr_alg; ++i)
        {
            *Cr_l += *ab_l;
            ab_l++;
            Cr_l += rs_c;
        }
    }
}

void gemm(unsigned int m,
          unsigned int n,
          unsigned int k,
		  uint32_t alphap,
          uint32_t *A,
          unsigned int lda,
          uint32_t *B,
          unsigned int ldb,
          uint32_t *C,
          unsigned int ldc,
          void * Ac_pack_v,
          void * Bc_pack_v)
{
	uint32_t *Ac, *Bc;
	uint32_t *Cc;
    uint32_t *Ar, *Br;
	uint32_t *Cr;
    uint32_t *Ac_pack=(uint32_t *)Ac_pack_v;
	uint32_t *Bc_pack=(uint32_t *)Bc_pack_v;

    for (unsigned int jc=0; jc<n; jc+=BLOCK_NC)
	{
		unsigned int n_alg=fmin(BLOCK_NC,n-jc);
		for (unsigned int pc=0; pc<k; pc+=BLOCK_KC)
		{
			unsigned int k_alg=fmin(BLOCK_KC,k-pc);

			Bc=&B[pc+jc*ldb];
            Pack_B(Bc, ldb, Bc_pack, k_alg, n_alg);  //PACK B
			for (unsigned int ic=0; ic<m; ic+=BLOCK_MC)
			{
				unsigned int m_alg=fmin(BLOCK_MC,m-ic);
				Ac=&A[ic+pc*lda];
				uint32_t *Ac_pack_local=Ac_pack; // Ac pack pointer per Loop 3 thread
                Pack_A(Ac,lda,(uint32_t*)Ac_pack_local,m_alg,k_alg); //PACK A
				Cc=&C[ic+jc*ldc];
                for(unsigned jr=0;jr<n_alg;jr+=BLOCK_NR)
				{
					unsigned int nr_alg=fmin(BLOCK_NR,n_alg-jr);
					for(unsigned int ir=0;ir<m_alg;ir+=BLOCK_MR)
					{
						unsigned int mr_alg=fmin(BLOCK_MR,m_alg-ir);
						Ar=&Ac_pack_local[ir*k_alg];
						Br=&Bc_pack[jr*k_alg];
						Cr=&Cc[ir+jr*ldc];

                        dgemm_ref(k_alg,mr_alg,nr_alg,alphap,Ar,Br,Cr,1,ldc);
					}
				}
			}
		}
	}
}

void mlt(unsigned int m,
         unsigned int k,
         uint32_t* A,
         uint32_t* tableA)
{
	for(unsigned int i=0;i<m;i++)
	{
		for(unsigned int j=0;j<k;j++)
		{
			((uint32_t*)A)[j*m + i] = tableA[i*k + j];
		}
	}
}

uint32_t * correlate(uint32_t* tableA,
               int tableAsize,
               int compressed_snp_size,
			   int group_size)
{
	int m=tableAsize, n=tableAsize, k=compressed_snp_size;
    long long int i;
    long long int tableCsize = (long long int)m*(long long int)n;
    int pm;     //return value used for assert

	void *Ac_pack_v=NULL, *Bc_pack_v=NULL, *A=NULL, *C=NULL;

    pm = posix_memalign(&(Ac_pack_v), 4096,
                        12*BLOCK_MC*BLOCK_KC*sizeof(uint32_t));
    assert(!pm);
    pm = posix_memalign(&(Bc_pack_v), 4096,
                        BLOCK_KC*BLOCK_NC*sizeof(uint32_t));
    assert(!pm);

    pm = posix_memalign(&A, 4096, m*k*sizeof(uint32_t) +
                         m*k*sizeof(uint32_t)%4096);
    assert(!pm);

    pm = posix_memalign(&C, 4096, tableCsize*sizeof(uint32_t) +
                        tableCsize*sizeof(uint32_t)%4096);
    assert(!pm);

	uint32_t  alphap=1.0;

	for(i=0;i<tableCsize;i++)
	{
		((uint32_t*)C)[i] = 0;
	}

    mlt(m, k, A, tableA);

	gemm(m,
			n,
			k,
			alphap,
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

// OPENCL STUFF GPU
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
                            1, &events[(ic_iter % 2)*4 + 2], NULL
                            );
                    printCLErr(err,__LINE__,__FILE__);
                }

                //                getc(stdin); //used to pause program for monitoring - MPAMPIS -

                clFinish(io_queue);
                clFinish(compute_queue);

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
               int compressed_snp_size,
			   int group_size)
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