
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"
#include "../include/finite_volume.h"



#define CONS_QTY_COPY(ifvU, cvU, c)				\
	do{											\
		ifvU##_rho = cvU##_rho[c];				\
		ifvU##_gamma = cvU##_gamma[c];			\
		ifvU##_e = cvU##_e[c];					\
		ifvU##_u = cvU##_u[c];					\
		if (dim > 1)							\
			ifvU##_v = cvU##_v[c];				\
		if (dim > 2)							\
			ifvU##_w = cvU##_w[c];				\
		if ((int)config[2] == 2)				\
			ifvU##_phi = cvU##_phi[c];			\
	} while (0)									\

int fluid_var_update(struct flu_var *FV, struct cell_var cv)
{
	const int dim = (int)config[0];		
	const int num_cell = (int)config[3];
	struct i_f_var ifv;
	
	for(int k = 0; k < num_cell; k++)
		{
			CONS_QTY_COPY(ifv.U, cv.U, k);
			
			if(cons2prim(&ifv) == 0)
				{
					fprintf(stderr, "Wrong in cons var to prim var!\n");
					return 0;
				}

			FV->RHO[k] = ifv.RHO;
			FV->P[k]   = ifv.P;
			FV->U[k]   = ifv.U;
			if (dim > 1)
				FV->V[k] = ifv.V;
			if (dim > 2)
				FV->W[k] = ifv.W;
			if ((int)config[2] == 2)
				FV->PHI[k] = ifv.PHI;
			FV->gamma[k] = ifv.gamma;
		}
	return 1;
}


static int order2_i_f_var_init(const struct cell_var cv, struct i_f_var * ifv, const int k)
{
	const int dim = (int)config[0];	
	const double n_x = ifv->n_x, n_y = ifv->n_y, n_z = ifv->n_z;
	const double delta_x = ifv->delta_x, delta_y = ifv->delta_y, delta_z = ifv->delta_z;
												
	ifv->d_rho  = cv.gradx_rho[k]*n_x;
	ifv->d_e    = cv.gradx_e[k]  *n_x;
	ifv->d_u    = cv.gradx_u[k]  *n_x;
	if ((int)config[2] == 2)									
		ifv->d_phi  = cv.gradx_phi[k]*n_x;
	if (dim > 1)
		{
			ifv->d_rho += cv.grady_rho[k]*n_y;
			ifv->d_e   += cv.grady_e[k]  *n_y;
			ifv->d_u   += cv.grady_u[k]  *n_y;
			ifv->d_v    = cv.gradx_v[k]  *n_x     + cv.grady_v[k]*n_y;
			if ((int)config[2] == 2)
				ifv->d_phi += cv.grady_phi[k]*n_y;
		}
	if (dim > 2)
		{
			ifv->d_rho += cv.gradz_rho[k]*n_z;
			ifv->d_e   += cv.gradz_e[k]  *n_z;
			ifv->d_u   += cv.gradz_u[k]  *n_z;
			ifv->d_v   += cv.gradz_v[k]  *n_z;
			ifv->d_w    = cv.gradx_w[k]  *n_x     + cv.grady_w[k]*n_y     + cv.gradz_w[k]*n_z;
			if ((int)config[2] == 2)
				ifv->d_phi += cv.gradz_phi[k]*n_z;
		}

	if(cons2prim(ifv) == 0)
		{
			fprintf(stderr, "Error happens on primitive variable!\n");
			return -1;
		}
	if (isinf(config[31]))
		{					
			ifv->d_p = ifv->d_e;

			ifv->RHO += cv.gradx_rho[k]*delta_x;
			ifv->P   += cv.gradx_e[k]  *delta_x;
			ifv->U   += cv.gradx_u[k]  *delta_x;
			if ((int)config[2] == 2)									
				ifv->PHI += cv.gradx_phi[k]*delta_x;
			if (dim > 1)
				{
					ifv->RHO += cv.grady_rho[k]*delta_y;
					ifv->P   += cv.grady_e[k]  *delta_y;
					ifv->U   += cv.grady_u[k]  *delta_y;
					ifv->V   += cv.gradx_v[k]  *delta_x + cv.grady_v[k]*delta_y;
					if ((int)config[2] == 2)				
						ifv->PHI += cv.grady_phi[k]*delta_y;				
				}
			if (dim > 2)
				{
					ifv->RHO += cv.gradz_rho[k]*delta_z;
					ifv->P   += cv.gradz_e[k]  *delta_z;
					ifv->U   += cv.gradz_u[k]  *delta_z;
					ifv->V   += cv.gradz_v[k]  *delta_z;
					ifv->W   += cv.gradx_w[k]  *delta_x + cv.grady_w[k]*delta_y + cv.gradz_w[k]*delta_z;
					if ((int)config[2] == 2)
						ifv->PHI += cv.gradz_phi[k]*delta_z;				
				}
		}
	else
		{
			ifv->d_u = (ifv->d_u - ifv->U*ifv->d_rho)/ifv->RHO;	
			ifv->d_p = (ifv->d_e - 0.5*ifv->d_rho*ifv->U*ifv->U - ifv->RHO*ifv->U*ifv->d_u) * (ifv->gamma-1.0);	
			if (dim > 1)
				{
					ifv->d_v  = (ifv->d_v - ifv->V*ifv->d_rho)/ifv->RHO;
					ifv->d_p += (- 0.5*ifv->d_rho*ifv->V*ifv->V - ifv->RHO*ifv->V*ifv->d_v) * (ifv->gamma-1.0);
				}
			if (dim > 2)
				{			
					ifv->d_w  = (ifv->d_w - ifv->W*ifv->d_rho)/ifv->RHO;
					ifv->d_p += (- 0.5*ifv->d_rho*ifv->W*ifv->W - ifv->RHO*ifv->W*ifv->d_w) * (ifv->gamma-1.0);
				}
			if ((int)config[2] == 2)
				ifv->d_phi = (ifv->d_phi - ifv->PHI*ifv->d_rho)/ifv->RHO;
			if (!isinf(config[60]))
				ifv->d_gamma = (ifv->d_gamma - ifv->gamma*ifv->d_rho)/ifv->RHO;

	ifv->d_rho= 0.0;
	ifv->d_e  = 0.0;
	ifv->d_u  = 0.0;
	if (dim > 1)						
		ifv->d_v = 0.0;
	if (dim > 2)
		ifv->d_w  = 0.0;				
	if ((int)config[2] == 2)
		ifv->d_phi = 0.0;
/*
			ifv->U_rho += cv.gradx_rho[k]*delta_x;
			ifv->U_e   += cv.gradx_e[k]  *delta_x;
			ifv->U_u   += cv.gradx_u[k]  *delta_x;
			if ((int)config[2] == 2)									
				ifv->U_phi += cv.gradx_phi[k]*delta_x;
			if (dim > 1)
				{
					ifv->U_rho += cv.grady_rho[k]*delta_y;
					ifv->U_e   += cv.grady_e[k]  *delta_y;
					ifv->U_u   += cv.grady_u[k]  *delta_y;
					ifv->U_v   += cv.gradx_v[k]  *delta_x + cv.grady_v[k]*delta_y;
					if ((int)config[2] == 2)				
						ifv->U_phi += cv.grady_phi[k]*delta_y;				
				}
			if (dim > 2)
				{
					ifv->U_rho += cv.gradz_rho[k]*delta_z;
					ifv->U_e   += cv.gradz_e[k]  *delta_z;
					ifv->U_u   += cv.gradz_u[k]  *delta_z;
					ifv->U_v   += cv.gradz_v[k]  *delta_z;
					ifv->U_w   += cv.gradx_w[k]  *delta_x + cv.grady_w[k]*delta_y + cv.gradz_w[k]*delta_z;
					if ((int)config[2] == 2)
						ifv->U_phi += cv.gradz_phi[k]*delta_z;				
				}
*/
			if(cons2prim(ifv) == 0)
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return -1;
				}
		}
	
	return 1;
}


static int order2_i_f_var0(struct i_f_var * ifv)
{
	const int dim = (int)config[0];
		
	ifv->d_rho= 0.0;
	ifv->d_e  = 0.0;
	ifv->d_u  = 0.0;
	if (dim > 1)						
		ifv->d_v = 0.0;
	if (dim > 2)
		ifv->d_w  = 0.0;				
	if ((int)config[2] == 2)
		ifv->d_phi = 0.0;

	if(cons2prim(ifv) == 0)
		{
			fprintf(stderr, "Error happens on primitive variable!\n");
			return -1;
		}
	
	return 1;
}

	  
int interface_var_init
(const struct cell_var cv, const struct mesh_var mv,
 struct i_f_var * ifv, struct i_f_var * ifv_R,
 const int k, const int j)
{
	const int dim = (int)config[0];
	const int order = (int)config[9];
	int **cc = cv.cell_cell;
	int **cp = mv.cell_pt;	

	int p_p, p_n;
	if (dim == 2)
		{	
			if(j == cp[k][0]-1) 
				{
					p_p=cp[k][1];
					p_n=cp[k][j+1];
				}				  
			else
				{
					p_p=cp[k][j+2];
					p_n=cp[k][j+1];
				}
			ifv->length = sqrt((mv.X[p_p] - mv.X[p_n])*(mv.X[p_p] - mv.X[p_n]) + (mv.Y[p_p] - mv.Y[p_n])*(mv.Y[p_p] - mv.Y[p_n]));
		}
	ifv->n_x = cv.n_x[k][j];
	if (dim > 1)	
		ifv->n_y = cv.n_y[k][j];
	if (dim > 2)
		ifv->n_z = cv.n_z[k][j];

	CONS_QTY_COPY(ifv->U, cv.U, k);
	
	if (order == 2)
		{
			if (dim == 1)
				ifv->delta_x = mv.X[cp[k][j+1]] - cv.X_c[k];
			else if (dim == 2)
				{					
					ifv->delta_x = 0.5*(mv.X[p_p] + mv.X[p_n]) - cv.X_c[k];
					ifv->delta_y = 0.5*(mv.Y[p_p] + mv.Y[p_n]) - cv.Y_c[k];
				}
			
			if(order2_i_f_var_init(cv, ifv, k) == -1)			
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return -1;
				}

		}
	
	int cR; //cell_right
	
	ifv_R->n_x = ifv->n_x;
	if (dim > 1)			
		ifv_R->n_y = ifv->n_y;
	if (dim > 2)
		ifv_R->n_z = ifv->n_z;
	
	if (cc[k][j] >= 0)
		{
			cR = cc[k][j];
			CONS_QTY_COPY(ifv_R->U, cv.U, cR);

			if (order == 2)
				{
					if (dim == 1)
						ifv_R->delta_x = mv.X[cp[k][j+1]] - cv.X_c[cR];
					else if (dim == 2)
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
							ifv_R->delta_x = 0.5*(mv.X[p_p] + mv.X[p_n]) - cv.X_c[cR];
							ifv_R->delta_y = 0.5*(mv.Y[p_p] + mv.Y[p_n]) - cv.Y_c[cR];
						}
					if(order2_i_f_var_init(cv, ifv_R, cR) == -1)
						{
							fprintf(stderr, "Error happens on primitive variable!\n");
							return -1;
						}
				}
		}
	else if (cc[k][j] == -1)//initial boundary condition.		
		{
			CONS_QTY_COPY(ifv_R->U, cv.U0, k);

			if (order == 2)
				if(order2_i_f_var0(ifv_R) == -1)
					{
						fprintf(stderr, "Error happens on primitive variable!\n");
						return -1;
					}
		}
	else if (cc[k][j] == -2)//reflecting boundary condition.
		{
			if(cons2prim(ifv) == 0)
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return -1;
				}
			ifv->F_rho = 0.0;
			ifv->F_u = ifv->P*ifv->n_x;
			ifv->F_v = ifv->P*ifv->n_y;
			ifv->F_e = 0.0;
			ifv->F_gamma = 0.0;
			return 0;
		}
	else if (cc[k][j] == -3)//prescribed boundary condition.
		{
			CONS_QTY_COPY(ifv_R->U, cv.U, k);

			if (order == 2)
				if(order2_i_f_var0(ifv_R) == -1)
					{
						fprintf(stderr, "Error happens on primitive variable!\n");
						return -1;
					}

		}		
	else
		{
			printf("No suitable boundary!\n");
			return -1;
		}

	if (order == 1)
		{
			if(cons2prim(ifv) == 0)
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return -1;
				}
			if(cons2prim(ifv_R) == 0)
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return -1;
				}
		}

	return 1;
}

double tau_calc(const struct cell_var cv, const struct mesh_var mv)
{
	const double CFL = config[7];
	if (CFL < 0.0)
		return -CFL;
	const int dim = (int)config[0];		
	const int num_cell = (int)config[3];
	int ** cp = mv.cell_pt;
	
	double tau = config[1];
	struct i_f_var ifv, ifv_R;
	double cum, lambda_max;
	int ivi;
	
	double qn, qn_R;
	double c, c_R;	
	
	for(int k = 0; k < num_cell; ++k)
		{
			cum = 0.0;
			
			for(int j = 0; j < cp[k][0]; ++j)
				{
					ivi = interface_var_init(cv, mv, &ifv, &ifv_R, k, j);
					if (ivi == 0)
						;
					else if(ivi == -1)
						return -1.0;
					else if (dim == 2)
						{
							qn = ifv.U*ifv.n_x + ifv.V*ifv.n_y; 
							qn_R = ifv_R.U*ifv_R.n_x + ifv_R.V*ifv_R.n_y;
							c = sqrt(ifv.gamma * ifv.P / ifv.RHO);
							c_R = sqrt(ifv_R.gamma * ifv_R.P / ifv_R.RHO);
							lambda_max = fmax(c+fabs(qn), c_R+fabs(qn_R));
							cum +=  lambda_max * ifv.length;
						}
				}		
			tau = fmin(tau, 2.0*cv.vol[k]/cum * CFL);
		}	//To decide tau.
	return tau;
}
