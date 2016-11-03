#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/finite_volume.h"



int fluid_var_update(struct flu_var *FV, struct cell_var cv)
{
	const int dim = (int)config[0];		
	const int num_cell = (int)config[3];
	struct i_f_var ifv;
	
	for(int k = 0; k < num_cell; k++)
		{
			cons_qty_copy_cv2ifv(&ifv, cv, k);			
			if(cons2prim(&ifv) == 0)
				{
					fprintf(stderr, "Wrong in copying cons_var to prim_var!\n");
					return 0;
				}
			prim_var_copy_ifv2FV(ifv, FV, k);
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
			ifv->d_v    = cv.gradx_v[k]  *n_x + cv.grady_v[k]*n_y;
			if ((int)config[2] == 2)
				ifv->d_phi += cv.grady_phi[k]*n_y;
		}
	if (dim > 2)
		{
			ifv->d_rho += cv.gradz_rho[k]*n_z;
			ifv->d_e   += cv.gradz_e[k]  *n_z;
			ifv->d_u   += cv.gradz_u[k]  *n_z;
			ifv->d_v   += cv.gradz_v[k]  *n_z;
			ifv->d_w    = cv.gradx_w[k]  *n_x + cv.grady_w[k]*n_y + cv.gradz_w[k]*n_z;
			if ((int)config[2] == 2)
				ifv->d_phi += cv.gradz_phi[k]*n_z;
		}

	if (cons2prim(ifv) == 0)
		{
			fprintf(stderr, "Error happens on primitive variable!\n");
			return 0;
		}
	if ((int)config[31] == 0)
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
	else if ((int)config[31] == 1)
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
			if(cons2prim(ifv) == 0)
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return 0;
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
			return 0;
		}
	
	return 1;
}

	  
int interface_var_init
(const struct cell_var cv, const struct mesh_var mv,
 struct i_f_var * ifv, struct i_f_var * ifv_R,
 const int k, const int j, const int i)
{
	const int dim = (int)config[0];
	const int order = (int)config[9];
	int **cc = cv.cell_cell;
	int **cp = mv.cell_pt;	

	int p_p, p_n;
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

	if (dim == 1)
		ifv->n_x = 1.0;
	if (dim > 1)
		{
			ifv->n_x = cv.n_x[k][j];
			ifv->n_y = cv.n_y[k][j];
		}
	if (dim > 2)
		ifv->n_z = cv.n_z[k][j];
		if (dim == 2)
		ifv->length = sqrt((mv.X[p_p] - mv.X[p_n])*(mv.X[p_p] - mv.X[p_n]) + (mv.Y[p_p] - mv.Y[p_n])*(mv.Y[p_p] - mv.Y[p_n]));

	cons_qty_copy_cv2ifv(ifv, cv, k);
	
	if (order == 2)
		{
			if (dim == 1)
				ifv->delta_x = mv.X[p_n] - cv.X_c[k];
			else if (dim == 2)
				{					
					ifv->delta_x = 0.5*(mv.X[p_p] + mv.X[p_n]) - cv.X_c[k];
					ifv->delta_y = 0.5*(mv.Y[p_p] + mv.Y[p_n]) - cv.Y_c[k];
				}			
			if(order2_i_f_var_init(cv, ifv, k) == 0)			
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return 0;
				}
		}
		
	ifv_R->n_x = ifv->n_x;
	if (dim > 1)
		{
			ifv_R->n_y = ifv->n_y;
			ifv_R->length = ifv->length;
		}
	if (dim > 2)
		ifv_R->n_z = ifv->n_z;
	
	int cR; //cell_right	
	if (cc[k][j] >= 0)
		{
			cR = cc[k][j];
			cons_qty_copy_cv2ifv(ifv_R, cv, cR);			

			if (order == 2)
				{
					if (dim == 1)
						ifv_R->delta_x = mv.X[p_n] - cv.X_c[cR];
					else if (dim == 2)
						{
							ifv_R->delta_x = 0.5*(mv.X[p_p] + mv.X[p_n]) - cv.X_c[cR];
							ifv_R->delta_y = 0.5*(mv.Y[p_p] + mv.Y[p_n]) - cv.Y_c[cR];
						}
					if(order2_i_f_var_init(cv, ifv_R, cR) == 0)
						{
							fprintf(stderr, "Error happens on primitive variable!\n");
							return 0;
						}
				}
		}
	else if (cc[k][j] == -1)//initial boundary condition.		
		{
			if (i > 0)
				return -1;
			cons_qty_copy_cv2ifv(ifv_R, cv, k);

			if (order == 2)
				if(order2_i_f_var0(ifv_R) == 0)
					{
						fprintf(stderr, "Error happens on primitive variable!\n");
						return 0;
					}
		}
	else if (cc[k][j] == -2)//reflecting boundary condition.
		{
			if(cons2prim(ifv) == 0)
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return 0;
				}
			ifv->F_rho = 0.0;
			ifv->F_u = ifv->P*ifv->n_x;
			if (dim > 1)
				ifv->F_v = ifv->P*ifv->n_y;
			if (dim > 2)
				ifv->F_w = ifv->P*ifv->n_z;
			ifv->F_e = 0.0;
			if (!isinf(config[60]))
				ifv->F_gamma = 0.0;
			if ((int)config[2] == 2)
				ifv->F_phi = 0.0;

			ifv->F_rho_star = ifv->F_rho;
			ifv->F_u_star = ifv->P*ifv->n_x;
			if (dim > 1)
				ifv->F_v_star = ifv->F_v;
			if (dim > 2)
				ifv->F_w_star = ifv->F_w;
			ifv->F_e_star = ifv->F_e;
			if (!isinf(config[60]))
				ifv->F_gamma_star = ifv->F_gamma;
			if ((int)config[2] == 2)
				ifv->F_phi_star = ifv->F_phi;
			ifv->F_rho_minus_c = ifv->F_rho;
			ifv->F_u_minus_c = ifv->P*ifv->n_x;
			if (dim > 1)
				ifv->F_v_minus_c = ifv->F_v;
			if (dim > 2)
				ifv->F_w_minus_c = ifv->F_w;
			ifv->F_e_minus_c = ifv->F_e;
			if (!isinf(config[60]))
				ifv->F_gamma_minus_c = ifv->F_gamma;
			if ((int)config[2] == 2)
				ifv->F_phi_minus_c = ifv->F_phi;
			ifv->F_rho_add_c = ifv->F_rho;
			ifv->F_u_add_c = ifv->P*ifv->n_x;
			if (dim > 1)
				ifv->F_v_add_c = ifv->F_v;
			if (dim > 2)
				ifv->F_w_add_c = ifv->F_w;
			ifv->F_e_add_c = ifv->F_e;
			if (!isinf(config[60]))
				ifv->F_gamma_add_c = ifv->F_gamma;
			if ((int)config[2] == 2)
				ifv->F_phi_add_c = ifv->F_phi;

			ifv->F_delta_e = 0.0;
			ifv->u_star = 0.0;
			ifv->u_minus_c = 0.0;
			ifv->u_add_c = 0.0;

			return -2;
		}
	else if (cc[k][j] == -3)//prescribed boundary condition.
		{
			cons_qty_copy_cv2ifv(ifv_R, cv, k);			

			if (order == 2)
				if(order2_i_f_var0(ifv_R) == 0)
					{
						fprintf(stderr, "Error happens on primitive variable!\n");
						return 0;
					}
		}		
	else
		{
			printf("No suitable boundary!\n");
			return 0;
		}

	if (order == 1)
		{
			if(cons2prim(ifv) == 0)
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return 0;
				}
			if(cons2prim(ifv_R) == 0)
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return 0;
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
					ivi = interface_var_init(cv, mv, &ifv, &ifv_R, k, j, 0);
					if (ivi < 0)
						;
					else if(ivi == 0)
						return -1.0;
					else if (dim == 1)
						{
							c = sqrt(ifv.gamma * ifv.P / ifv.RHO);
							c_R = sqrt(ifv_R.gamma * ifv_R.P / ifv_R.RHO);
							lambda_max = fmax(c+fabs(ifv.U), c_R+fabs(ifv_R.U));
							cum += lambda_max;	
						}
					else if (dim == 2)
						{
							qn = ifv.U*ifv.n_x + ifv.V*ifv.n_y; 
							qn_R = ifv_R.U*ifv_R.n_x + ifv_R.V*ifv_R.n_y;
							c = sqrt(ifv.gamma * ifv.P / ifv.RHO);
							c_R = sqrt(ifv_R.gamma * ifv_R.P / ifv_R.RHO);
							lambda_max = fmax(c+fabs(qn), c_R+fabs(qn_R));
							cum += 0.5*lambda_max * ifv.length;
						}
				}
			tau = fmin(tau, cv.vol[k]/cum * CFL);
		}	//To decide tau.
	return tau;
}
