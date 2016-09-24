#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"



#define copy(var)  cv->var[i] = cv->var[pc[i]]


void period_ghost(struct cell_var * cv, struct mesh_var mv, double t)
{
	const int dim = (int)config[0];
	const int order = (int)config[9];
	const int num_cell = (int)config[3];
	const int num_cell_ghost = mv.num_ghost + num_cell;
	const int *pc = mv.peri_cell;	
	
	for(int i = num_cell; i < num_cell_ghost; i++)
		{
			copy(U_rho);
			copy(U_gamma);
			copy(U_e);
			copy(U_u);
			if (order > 1)
				{
					copy(gradx_rho);
					copy(gradx_e);
					copy(gradx_u);;
				}
			if (dim > 1)
				{
					copy(U_v);
					if (order > 1)
						{
							copy(grady_rho);
							copy(grady_e);
							copy(grady_u);
							copy(grady_v);
							copy(gradx_v);
						}
				}
			if (dim > 2)
				{
					copy(U_w);
					if (order > 1)
						{
							copy(gradz_rho);
							copy(gradz_e);
							copy(gradz_u);
							copy(gradz_v);
							copy(gradz_w);
							copy(grady_w);
							copy(gradx_w);
						}
				}
			if ((int)config[2] == 2)
				{
					copy(U_phi);
					if (order > 1)
						{
							copy(gradx_phi);
							if (dim > 2)
								copy(grady_phi);
							if (dim > 3)
								copy(gradz_phi);
						}
				}
		}
}



void period_cell_modi(struct mesh_var * mv)
{
	const int num_cell = mv->num_ghost + (int)config[3];
	int *pc = mv->peri_cell;
	
	int per_num[num_cell], per_n = 0;
	int i, j;
	for (i = 0; i < num_cell; i++)
		{
			if (pc[i] >= 0)
				per_n++;
			per_num[i] = per_n;
		}

	for (i = 0; i < num_cell; i++)
		if (pc[i] >= 0)
			pc[i] -= per_num[pc[i]];

	int *cc_tmp, pc_tmp; 
	for (i = 1; i < num_cell; i++)
		{
			j = i;
			for(j = i; pc[j-1] >= 0 && j >= 1; j--)
				{
					if(pc[j] < 0)
						{
							pc_tmp = pc[j-1];
							pc[j-1] = pc[j];
							pc[j] = pc_tmp;
							cc_tmp = mv->cell_pt[j];
							mv->cell_pt[j] = mv->cell_pt[j-1];
							mv->cell_pt[j-1] = cc_tmp;
						}
					else
						break;
				}	 				
		}
	mv->bc = period_ghost;
}
