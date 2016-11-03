#include <stdio.h>
#include <math.h>

#include "../include/var_struc.h"



void cons_qty_copy_cv2ifv(struct i_f_var *ifv, struct cell_var cv, int c)
{
	const int dim = (int)config[0];
	
	ifv->U_rho = cv.U_rho[c];
	ifv->U_e = cv.U_e[c];
	ifv->U_u = cv.U_u[c];
	if (dim > 1)
		ifv->U_v = cv.U_v[c];
	if (dim > 2)
		ifv->U_w = cv.U_w[c];
	if ((int)config[2] == 2)
		ifv->U_phi = cv.U_phi[c];
	ifv->U_gamma = cv.U_gamma[c];

	ifv->delta_U_e = cv.delta_U_e[c];
}

void prim_var_copy_ifv2FV(struct i_f_var ifv, struct flu_var *FV, int c)
{
	const int dim = (int)config[0];

	FV->RHO[c] = ifv.RHO;
	FV->P[c]   = ifv.P;
	FV->U[c]   = ifv.U;
	if (dim > 1)
		FV->V[c] = ifv.V;
	if (dim > 2)
		FV->W[c] = ifv.W;
	if ((int)config[2] == 2)
		FV->PHI[c] = ifv.PHI;
	FV->gamma[c] = ifv.gamma;
}

void flux_copy_ifv2cv(struct i_f_var ifv, struct cell_var *cv, int k, int j)
{
	const int dim = (int)config[0];
	const int order = (int)config[9];
	
	cv->F_rho[k][j] = ifv.F_rho;
	cv->F_e[k][j]   = ifv.F_e;
	cv->F_u[k][j]   = ifv.F_u;
	if (dim > 1)
		cv->F_v[k][j] = ifv.F_v;
	if (dim > 2)
		cv->F_w[k][j] = ifv.F_w;
	if ((int)config[2] == 2)
		cv->F_phi[k][j] = ifv.F_phi;
	if (!isinf(config[60]))
		cv->F_gamma[k][j] = ifv.F_gamma;

	cv->U_p[k][j]   = ifv.U;
	if (dim > 1)
		cv->V_p[k][j] = ifv.V;
	cv->F_p_x[k][j] = ifv.P;
	cv->RHO_p[k][j] = ifv.RHO;	
	if ((int)config[2] == 2)
		cv->PHI_p[k][j] = ifv.PHI;
	cv->gamma_p[k][j] = ifv.gamma;


	cv->F_rho_star[k][j] = ifv.F_rho_star;
	cv->F_e_star[k][j]   = ifv.F_e_star;
	cv->F_u_star[k][j]   = ifv.F_u_star;
	if (dim > 1)
		cv->F_v_star[k][j] = ifv.F_v_star;
	if (dim > 2)
		cv->F_w_star[k][j] = ifv.F_w_star;
	if ((int)config[2] == 2)
		cv->F_phi_star[k][j] = ifv.F_phi_star;
	if (!isinf(config[60]))
		cv->F_gamma_star[k][j] = ifv.F_gamma_star;
	cv->F_rho_add_c[k][j] = ifv.F_rho_add_c;
	cv->F_e_add_c[k][j]   = ifv.F_e_add_c;
	cv->F_u_add_c[k][j]   = ifv.F_u_add_c;
	if (dim > 1)
		cv->F_v_add_c[k][j] = ifv.F_v_add_c;
	if (dim > 2)
		cv->F_w_add_c[k][j] = ifv.F_w_add_c;
	if ((int)config[2] == 2)
		cv->F_phi_add_c[k][j] = ifv.F_phi_add_c;
	if (!isinf(config[60]))
		cv->F_gamma_add_c[k][j] = ifv.F_gamma_add_c;
	cv->F_rho_minus_c[k][j] = ifv.F_rho_minus_c;
	cv->F_e_minus_c[k][j]   = ifv.F_e_minus_c;
	cv->F_u_minus_c[k][j]   = ifv.F_u_minus_c;
	if (dim > 1)
		cv->F_v_minus_c[k][j] = ifv.F_v_minus_c;
	if (dim > 2)
		cv->F_w_minus_c[k][j] = ifv.F_w_minus_c;
	if ((int)config[2] == 2)
		cv->F_phi_minus_c[k][j] = ifv.F_phi_minus_c;
	if (!isinf(config[60]))
		cv->F_gamma_minus_c[k][j] = ifv.F_gamma_minus_c;

	cv->u_star[k][j]  = ifv.u_star;
	cv->u_minus_c[k][j] = ifv.u_minus_c;
	cv->u_add_c[k][j] = ifv.u_add_c;

	cv->F_delta_e[k][j] = ifv.F_delta_e;
}
