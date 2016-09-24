#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"


static int line_mesh(struct mesh_var * mv)
{
	int k;
	const int num_cell = (int)config[3];

	if(isinf(config[13]))
		{
			fprintf(stderr, "The initial data is not mentioned in the 1-D mesh!\n");
			exit(2);
		}
	if(isinf(config[10]) || config[10] < config[4])
		{
			fprintf(stderr, "Without a proper spatial grid size!\n");
			exit(2);
		}
	
	const int n_x = (int)config[13];

	mv->num_pt = n_x+1;
	mv->X = malloc(mv->num_pt * sizeof(double));
	if(mv->X == NULL)
		{
			printf("Not enough memory in 1-D mesh constructed!\n");
			goto return_0;
		}
	for(k = 0; k < mv->num_pt; ++k)
		mv->X[k] = k * config[10];		
	
	mv->cell_pt = malloc(num_cell * sizeof(void *));
	if(mv->cell_pt == NULL)
		{
			fprintf(stderr, "Not enough memory in 1-D mesh constructed!\n");
			goto return_0;
		}
	for(k = 0; k < num_cell; ++k)
		{
			mv->cell_pt[k] = malloc(3 * sizeof(int));
			if(mv->cell_pt[k] == NULL)
				{
					for(int i = 0; i < k; ++i)
						{
							free(mv->cell_pt[i]);
							mv->cell_pt[i] = NULL;
						}
					fprintf(stderr, "Not enough memory in CELL_POINT[%d]!\n", k);
					goto return_0;
				}
			
			mv->cell_pt[k][0] = 2;
			mv->cell_pt[k][1] = k;
			mv->cell_pt[k][2] = mv->cell_pt[k][1] + 1;
		}

	mv->num_border[0] = 1;	
	mv->num_border[1] = 2;	
	mv->border_pt = malloc(2 * sizeof(int));	
	mv->border_cond = malloc(2 * sizeof(int));
	if(mv->border_pt == NULL || mv->border_cond == NULL)
		{
			printf("Not enough memory in square mesh constructed!\n");
			goto return_0;
		}
	
	mv->border_pt[0] = 0;
	mv->border_pt[1] = num_cell;

	return 1;

 return_0:
	free(mv->X);
	mv->X = NULL;	
	free(mv->border_pt);
	mv->border_pt = NULL;
	free(mv->border_cond);
	mv->border_cond = NULL;	
	if (mv->cell_pt != NULL)
		{
			for(int i = 0; i < num_cell; ++i)
				{
					free(mv->cell_pt[i]);
					mv->cell_pt[i] = NULL;
				}
			free(mv->cell_pt);
			mv->cell_pt = NULL;
		}
	return 0;	
}


void free_1D_mesh(struct mesh_var * mv)
{
	if (line_mesh(mv) == 0)
		exit(5);

	mv->border_cond[0] = -3;
	mv->border_cond[1] = -3;
}
