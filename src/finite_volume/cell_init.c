#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/finite_volume.h"


#define MAX(a,b) ((a) > (b) ? (a) : (b))


#define cv_init_mem(v, n)												\
	do {																\
		cv.v = malloc((n) * sizeof(double));							\
		if(cv.v == NULL)												\
			{															\
				fprintf(stderr, "Not enough memory in cell var init!\n"); \
				goto return_NULL;										\
			}															\
	} while (0)															\

#define cp_init_mem(v, n)												\
	do {																\
		cv.v = malloc((n) * sizeof(void *));							\
		if(cv.v == NULL)												\
			{															\
				fprintf(stderr, "Not enough memory in cell var init!\n"); \
				goto return_NULL;										\
			}															\
		init_mem(cv.v, n, mv.cell_pt);									\
	} while (0)															\

#define cp_init_mem_int(v, n)											\
	do {																\
		cv.v = malloc((n) * sizeof(void *));							\
		if(cv.v == NULL)												\
			{															\
				fprintf(stderr, "Not enough memory in cell var init!\n"); \
				goto return_NULL;										\
			}															\
		init_mem_int(cv.v, n, mv.cell_pt);								\
	} while (0)															\

struct cell_var cell_mem_init(const struct mesh_var mv)
{
	const int dim = (int)config[0];
	const int order = (int)config[9];
	const int num_cell_ghost = mv.num_ghost + (int)config[3];
	const int num_cell = (int)config[3];
	
	struct cell_var cv;

	cp_init_mem_int(cell_cell, num_cell_ghost);
	cv_init_mem(vol, num_cell_ghost);	
	cp_init_mem(n_x, num_cell_ghost);
	
	cp_init_mem(F_u, num_cell);
	cv_init_mem(U_u, num_cell_ghost);
	cv_init_mem(U0_u, num_cell);

	cp_init_mem(F_rho, num_cell);
	cv_init_mem(U_rho, num_cell_ghost);
	cv_init_mem(U0_rho, num_cell);

	cp_init_mem(F_e, num_cell);
	cv_init_mem(U_e, num_cell_ghost);
	cv_init_mem(U0_e, num_cell);

	if (order > 1)
		{
			cv_init_mem(X_c, num_cell_ghost);
			cv_init_mem(gradx_rho, num_cell_ghost);
			cv_init_mem(gradx_e, num_cell_ghost);
			cv_init_mem(gradx_u, num_cell_ghost);			
		}
	
	if (dim > 1)
		{
			cp_init_mem(n_y, num_cell_ghost);
			cp_init_mem(F_v, num_cell);
			cv_init_mem(U_v, num_cell_ghost);
			cv_init_mem(U0_v, num_cell);
			if (order > 1)
				{
					cv_init_mem(Y_c, num_cell_ghost);
					cv_init_mem(grady_rho, num_cell_ghost);
					cv_init_mem(grady_e, num_cell_ghost);
					cv_init_mem(grady_u, num_cell_ghost);
					cv_init_mem(grady_v, num_cell_ghost);
					cv_init_mem(gradx_v, num_cell_ghost);
				}
		}
	if (dim > 2)
		{					
			cp_init_mem(n_z, num_cell_ghost);
			cp_init_mem(F_w, num_cell);
			cv_init_mem(U_w, num_cell_ghost);
			cv_init_mem(U0_w, num_cell);		
			if (order > 1)
				{
					cv_init_mem(Z_c, num_cell_ghost);
					cv_init_mem(gradz_rho, num_cell_ghost);
					cv_init_mem(gradz_e, num_cell_ghost);
					cv_init_mem(gradz_u, num_cell_ghost);
					cv_init_mem(gradz_v, num_cell_ghost);
					cv_init_mem(gradz_w, num_cell_ghost);
					cv_init_mem(grady_w, num_cell_ghost);
					cv_init_mem(gradx_w, num_cell_ghost);
				}
		}

	if ((int)config[2] == 2)
		{	
			cp_init_mem(F_phi, num_cell);
			cv_init_mem(U_phi, num_cell_ghost);
			cv_init_mem(U0_phi, num_cell);				
			if (order > 1)
				{
					cv_init_mem(gradx_phi, num_cell_ghost);
					if (dim > 2)
						cv_init_mem(grady_phi, num_cell_ghost);
					if (dim > 3)
						cv_init_mem(gradz_phi, num_cell_ghost);	
				}
		}
	cp_init_mem(F_gamma, num_cell);
	cv_init_mem(U_gamma, num_cell_ghost);
	cv_init_mem(U0_gamma, num_cell);
	if (!isinf(config[60]))
		{					
			if (order > 1)
				{
					cv_init_mem(gradx_gamma, num_cell_ghost);
					if (dim > 2)
						cv_init_mem(grady_gamma, num_cell_ghost);
					if (dim > 3)
						cv_init_mem(gradz_gamma, num_cell_ghost);	
				}
		}

	return cv;
	
 return_NULL:
	exit(5);
}


//Calculate volume.
void vol_comp(struct cell_var * cv, const struct mesh_var mv)
{
	const int dim = (int)config[0];
	const int num_cell = mv.num_ghost + (int)config[3];
	int **cp = mv.cell_pt;
	
	int p_p, p_n;

	if (dim == 1)
		for(int k = 0; k < num_cell; k++)
			{
				p_p = cp[k][1];
				p_n = cp[k][2];					
				cv->vol[k] = fabs(mv.X[p_p] - mv.X[p_n]);						
			}
	else if (dim == 2)
		for(int k = 0; k < num_cell; k++)
			{			
				cv->vol[k] = 0.0;
				for(int j = 0; j < cp[k][0]; j++)
					{
						if(j == cp[k][0]-1) 
							{
								p_p = cp[k][1];
								p_n = cp[k][j+1];
							}				  
						else
							{
								p_p = cp[k][j+2];
								p_n = cp[k][j+1];
							} 
						cv->vol[k] = cv->vol[k] + 0.5 * (mv.X[p_n]*mv.Y[p_p] - mv.Y[p_n]*mv.X[p_p]);
					}
			}
}



//Determine the normal direction and relationship between cells.
void cell_rel(struct cell_var * cv, const struct mesh_var mv)
{
	const int dim = (int)config[0];
	const int num_cell = mv.num_ghost + (int)config[3];
	
	int **cp = mv.cell_pt;
	int p_p, p_n, p2_p, p2_n;

	int cell_rec, n_border;
	int i, l, ts;
	double length;

	if (dim == 1)
		for(int k = 0; k < num_cell; k++)
			{ 						
				for(int j = 0; j < 2; j++)
					{
						if(j == 1) 
							{
								p_p = cp[k][1];
								p_n = cp[k][2];
							}				  
						else
							{
								p_p = cp[k][2];
								p_n = cp[k][1];
							}						
						length = fabs(mv.X[p_p] - mv.X[p_n]);
						cv->n_x[k][j] = (mv.X[p_p] - mv.X[p_n]) / length;
				
						cell_rec = 0;
						ts = 1;
						while (ts <= MAX(num_cell-k-1, k))
							{
								// seek in two side
								i = k + ts;
								if (ts > 0)
									ts = -ts;
								else
									ts = -ts + 1;
								if (i < 0 || i >= num_cell)
									continue;
								
								for(l = 0; l < 2; l++)
									{
										if(l == 1) 
											{
												p2_p = cp[i][1];
												p2_n = cp[i][2];
											}				  
										else
											{
												p2_p = cp[i][2];
												p2_n = cp[i][1];
											}		
										if(p_p == p2_p || p_p == p2_n)
											{
												cv->cell_cell[k][j] = i;
												cell_rec = 1;
												break;
											}
									}
								if (cell_rec)
									break;
							}				
						if (cell_rec)
							continue;
						
						for(i = 0; i < 2; i++)
							{				
								p2_p = mv.border_pt[i];
								if(p_p == p2_p)
									{
										cv->cell_cell[k][j] = mv.border_cond[i];
										cell_rec = 1;
										break;
									}
							}
								
						if(!cell_rec && k < (int)config[3])
							{
								fprintf(stderr, "Ther are some wrong cell relationships!\n");
								exit(2);
							}
					}
			}
	else if (dim == 2)
		for(int k = 0; k < num_cell; k++)
			{ 						
				for(int j = 0; j < cp[k][0]; j++)
					{
						if(j == cp[k][0]-1) 
							{
								p_p = cp[k][1];
								p_n = cp[k][j+1];
							}				  
						else
							{
								p_p = cp[k][j+2];
								p_n = cp[k][j+1];
							}
						length = sqrt((mv.Y[p_p] - mv.Y[p_n])*(mv.Y[p_p] - mv.Y[p_n])+(mv.X[p_n] - mv.X[p_p])*(mv.X[p_n] - mv.X[p_p]));
						cv->n_x[k][j] = (mv.Y[p_p] - mv.Y[p_n]) / length;
						cv->n_y[k][j] = (mv.X[p_n] - mv.X[p_p]) / length;
					
						cell_rec = 0;							   		
						ts = 1;
						while (ts <= MAX(num_cell-k-1, k))
							{
								// seek in two side
								i = k + ts;
								if (ts > 0)
									ts = -ts;
								else
									ts = -ts + 1;
								if (i < 0 || i >= num_cell)
									continue;
								
								for(l = 0; l < cp[i][0]; l++)
									{
										if(l == cp[i][0]-1) 
											{
												p2_p = cp[i][1];
												p2_n = cp[i][l+1];
											}				  
										else
											{
												p2_p = cp[i][l+2];
												p2_n = cp[i][l+1];
											}
										if((p_p == p2_n) && (p2_p == p_n))
											{
												cv->cell_cell[k][j] = i;
												cell_rec = 1;;
												break;
											}
									}
								if (cell_rec)
									break;
							}
				
						if (cell_rec)
							continue;								
					
						for(l = 1, n_border = -1; l <= mv.num_border[0]; l++)
							{
								n_border += mv.num_border[l] + 1; 
								for(i = n_border-mv.num_border[l]; i < n_border; i++)
									{				
										p2_p = mv.border_pt[i+1];
										p2_n = mv.border_pt[i];
										if((p_p == p2_p && p_n == p2_n) || (p_p == p2_n && p_n == p2_p))
											{
												cv->cell_cell[k][j] = mv.border_cond[i];
												cell_rec = 1;
												break;
											}
									}							
								if (cell_rec)
									break;
							}
			
						if(!cell_rec && k < (int)config[3])
							{
								fprintf(stderr, "Ther are some wrong cell relationships!\n");
								exit(2);
							}
					}							
			}
}


void cell_centroid(struct cell_var * cv, const struct mesh_var mv)
{
	const int dim = (int)config[0];
	const int num_cell = mv.num_ghost + (int)config[3];
	const double *X = mv.X, *Y = mv.Y;
	int **cp = mv.cell_pt;
	
	double S, S_tri;

	if (dim == 1)		
		for(int k = 0; k < num_cell; ++k)
			cv->X_c[k] = (X[cp[k][1]] + X[cp[k][2]])/2.0;							
	else if (dim == 2)
		for(int k = 0; k < num_cell; ++k)
			{
				S = 0.0;
				cv->X_c[k] = 0.0;
				cv->Y_c[k] = 0.0;

				for(int j = 2; j < cp[k][0]; j++)
					{
						S_tri = X[cp[k][1]]*Y[cp[k][j]] + X[cp[k][j+1]]*Y[cp[k][1]] + X[cp[k][j]]*Y[cp[k][j+1]] - X[cp[k][j+1]]*Y[cp[k][j]] - X[cp[k][1]]*Y[cp[k][j+1]] - X[cp[k][j]]*Y[cp[k][1]];
						cv->X_c[k] += (X[cp[k][1]] + X[cp[k][j]] + X[cp[k][j+1]]) * S_tri;
						cv->Y_c[k] += (Y[cp[k][1]] + Y[cp[k][j]] + Y[cp[k][j+1]]) * S_tri;
						S += S_tri;
					}			 
				cv->X_c[k] /= S*3.0;
				cv->Y_c[k] /= S*3.0;
			}
}
