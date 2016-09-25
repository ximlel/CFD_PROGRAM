#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/finite_volume.h"


void Euler_scheme(struct flu_var *FV, const struct mesh_var mv, const char *scheme)
{
	clock_t start_clock;
	double cpu_time = 0.0;	

	const int dim = (int)config[0];
	const int order = (int)config[9];
	const int el = isinf(config[8]) ? 0 : (int)config[8];
	const int num_cell = (int)config[3];
	int ** cp = mv.cell_pt;
	
	struct cell_var cv = cell_mem_init(mv, FV);

	cons_qty_init(&cv, *FV);

	vol_comp(&cv, mv);

	cell_rel(&cv, mv);

	if (order > 1)
		cell_centroid(&cv, mv);

	printf("Grid has been constructed.\n");
	
	
	double tau; // the length of the time step
	double t_all = 0.0;
	struct i_f_var ifv, ifv_R;

	int k, j, ivi, stop_step = 0;
	for(int i = 0; i < (int)config[5] && stop_step == 0; ++i)
		{
			start_clock = clock();
			
			if (order > 1)
				{
					if (el != 0 && i > 0)
						cell_centroid(&cv, mv);
					if (isinf(config[31]))
						fluid_var_update(FV, cv);
					if (mv.bc != NULL)
						mv.bc(&cv, mv, FV, t_all);		
					slope_limiter(&cv, mv, *FV);
				}
			
			if (mv.bc != NULL)
				mv.bc(&cv, mv, FV, t_all);
			
			if (dim == 2)
				tau = tau_calc(cv, mv);

			t_all += tau;			
			if(tau < 0.00001)
				{
					printf("\nThe length of the time step is so small at step %d, t_all=%lf, tau=%lf.\n",i,t_all,tau);
					stop_step = 1;
				}
			
			if(t_all > config[1])
				{
					printf("\nThe time is enough at step %d.\n",i);
					tau = tau - (t_all - config[1]);
					t_all = config[1];
					stop_step = 1;
				} // Time

			for(k = 0; k < num_cell; k++)
				{
					if (stop_step == 1)
						break;
					for(j = 0; j < cp[k][0]; j++)
						{
							ivi = interface_var_init(cv, mv, &ifv, &ifv_R, k, j);
							if (ivi == 0)
								;
							else if(ivi == -1)
								{
									stop_step = 1;
									break;									
								}
							else if (order == 1)
								{													
									if (strcmp(scheme,"Roe") == 0)
										Roe_scheme(&ifv, &ifv_R);
									else if (strcmp(scheme,"HLL") == 0)
										HLL_scheme(&ifv, &ifv_R);
									else if(strcmp(scheme,"Riemann_exact") == 0)
										Riemann_exact_scheme(&ifv, &ifv_R);
									else										
										{
											printf("No Riemann solver!\n");
											exit(2);
										}
								}
							else if (order == 2)
								{													
									if(strcmp(scheme,"GRP") == 0)
										GRP_scheme(&ifv, &ifv_R, tau);
									else										
										{
											printf("No Riemann solver!\n");
											exit(2);
										}
								}


							cv.F_rho[k][j] = ifv.F_rho;
							cv.F_e[k][j]   = ifv.F_e;
							cv.F_u[k][j]   = ifv.F_u;
							if (dim > 1)
								cv.F_v[k][j] = ifv.F_v;
							if (dim > 2)
								cv.F_w[k][j] = ifv.F_w;
							if ((int)config[2] == 2)
								cv.F_phi[k][j] = ifv.F_phi;
							if(!isinf(config[60]))
								cv.F_gamma[k][j] = ifv.F_gamma;
						}					
				}

			cons_qty_update(&cv, mv, tau);			

			DispPro(t_all*100.0/config[1], i);
						
			cpu_time += (clock() - start_clock) / (double)CLOCKS_PER_SEC;

			if (stop_step == 1)
				break;				
		}

	fluid_var_update(FV, cv);

	printf("\nThe cost of CPU time for the Eulerian method is %g seconds.\n", cpu_time);
}
