#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/Riemann_solver.h"
#include "../include/var_struc.h"

void Roe_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R)
{
	const double delta = 0.2;
	
	double F[4];
	double lambda_max;

	ROE_2D_solver(F, ifv->gamma, ifv->P, ifv->RHO, ifv->U, ifv->V, ifv->n_x, ifv->n_y, ifv_R->P, ifv_R->RHO, ifv_R->U, ifv_R->V, &lambda_max, delta);
	ifv->F_rho = F[0];
	ifv->F_u   = F[1];
	ifv->F_v   = F[2];
	ifv->F_e   = F[3];
}

void HLL_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R)
{	
	double F[4];
	double lambda_max;

	HLL_solver(F, ifv->gamma, ifv->P, ifv->RHO, ifv->U, ifv->V, ifv->n_x, ifv->n_y, ifv_R->P, ifv_R->RHO, ifv_R->U, ifv_R->V, &lambda_max);
	ifv->F_rho = F[0];
	ifv->F_u   = F[1];
	ifv->F_v   = F[2];
	ifv->F_e   = F[3];
}

void Riemann_exact_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R)
{
	const double eps = config[4];	

	const double n_x = ifv->n_x, n_y = ifv->n_y;
	const double u   = ifv->U  *n_x + ifv->V  *n_y; 
	const double u_R = ifv_R->U*n_x + ifv_R->V*n_y;

	double dire[3], mid[3];
		
	linear_GRP_solver_Edir(dire, mid, ifv->RHO, ifv_R->RHO, 0.0, 0.0, u, u_R, 0.0, 0.0, ifv->P, ifv_R->P, 0.0, 0.0, ifv->gamma, eps);

	double rho_mid, p_mid, u_mid, v_mid, mid_qt;
	rho_mid = mid[0];
	p_mid = mid[2];

	if(mid[1]>0)
		mid_qt = -ifv->U  *n_y + ifv->V  *n_x;
	else
		mid_qt = -ifv_R->U*n_y + ifv_R->V*n_x;
	u_mid = mid[1]*n_x - mid_qt*n_y;
	v_mid = mid[1]*n_y + mid_qt*n_x;
											
	ifv->F_rho = rho_mid*(u_mid*n_x + v_mid*n_y);
	ifv->F_u   = ifv->F_rho*u_mid + p_mid*n_x;
	ifv->F_v   = ifv->F_rho*v_mid + p_mid*n_y;
	ifv->F_e   = (ifv->gamma/(ifv->gamma-1.0))*p_mid/rho_mid + 0.5*(u_mid*u_mid + v_mid*v_mid);
	ifv->F_e   = ifv->F_rho*ifv->F_e;
	if (!isinf(config[60]))
		ifv->F_gamma = ifv->F_rho*ifv->gamma;
}


void GRP_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R, double tau)
{
	const double eps = config[4];	
	const double n_x =  ifv->n_x, n_y = ifv->n_y;

	const double u     =  ifv->U    *n_x + ifv->V    *n_y; 
	const double u_R   =  ifv_R->U  *n_x + ifv_R->V  *n_y;
	const double d_u   =  ifv->d_u  *n_x + ifv->d_v  *n_y; 
	const double d_u_R =  ifv_R->d_u*n_x + ifv_R->d_v*n_y;
	const double v     = -ifv->U    *n_y + ifv->V    *n_x;
	const double v_R   = -ifv_R->U  *n_y + ifv_R->V  *n_x;
	const double d_v   = -ifv->d_u  *n_y + ifv->d_v  *n_x;
	const double d_v_R = -ifv_R->d_u*n_y + ifv_R->d_v*n_x;
	
	double wave_speed[2], dire[4], mid[4];
	double gamma = ifv->gamma;

	linear_GRP_2D_solver_Edir(wave_speed, dire, mid, 0.0, ifv->RHO, ifv->RHO, 0, 0, u, u_R, 0, 0, v, v_R, 0, 0, ifv->P, ifv_R->P, 0, 0, gamma, eps);
//	linear_GRP_2D_solver_Edir(wave_speed, dire, mid, 0.0, ifv->RHO, ifv->RHO, ifv->d_rho, ifv_R->d_rho, u, u_R, d_u, d_u_R, v, v_R, d_v, d_v_R, ifv->P, ifv_R->P, ifv->d_p, ifv_R->d_p, gamma, eps);

	double rho_mid, p_mid, u_mid, v_mid, mid_qt;
	rho_mid = mid[0] + 0.5*tau*dire[0];
	u_mid   = (mid[1] + 0.5*tau*mid[0]*dire[1]/rho_mid)*n_x - (mid[2] + 0.5*tau*mid[0]*dire[2]/rho_mid)*n_y;
	v_mid   = (mid[1] + 0.5*tau*mid[0]*dire[1]/rho_mid)*n_y + (mid[2] + 0.5*tau*mid[0]*dire[2]/rho_mid)*n_x;
	p_mid   = mid[3] + (gamma-1.0)*0.5*(tau*(dire[0]*0.5*(mid[1]*mid[1]+mid[2]*mid[2])+mid[0]*(mid[1]*dire[1]+mid[2]*dire[2])) + mid[0]*mid[1]*mid[1]+mid[0]*mid[2]*mid[2]-rho_mid*u_mid*u_mid-rho_mid*v_mid*v_mid);
										   
	ifv->F_rho = rho_mid*(u_mid*n_x + v_mid*n_y);
	ifv->F_u   = ifv->F_rho*u_mid + p_mid*n_x;
	ifv->F_v   = ifv->F_rho*v_mid + p_mid*n_y;
	ifv->F_e   = (gamma/(gamma-1.0))*p_mid/rho_mid + 0.5*(u_mid*u_mid + v_mid*v_mid);
	ifv->F_e   = ifv->F_rho*ifv->F_e;
}