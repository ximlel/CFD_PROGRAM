#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"


//Initialize conserved quantities.
void cons_qty_init(struct cell_var * cv, const struct flu_var FV)
{
	const int dim = (int)config[0];
	const int num_cell = (int)config[3];
	for(int k = 0; k < num_cell; k++)
		{
			cv->U_rho[k]   = FV.RHO[k];
			cv->U_gamma[k] = FV.RHO[k] * FV.gamma[k];			
			cv->U_e[k]     = FV.P[k]/(FV.gamma[k]-1.0) + 0.5*FV.RHO[k]*FV.U[k]*FV.U[k];
			cv->U_u[k]     = FV.RHO[k] * FV.U[k];			
			if (dim > 1)
				{									
					cv->U_v[k]  = FV.RHO[k] * FV.V[k];
					cv->U_e[k] += 0.5*FV.RHO[k]*FV.V[k]*FV.V[k];
				}
			if (dim > 2)
				{									
					cv->U_w[k]  = FV.RHO[k] * FV.W[k];
					cv->U_e[k] += 0.5*FV.RHO[k]*FV.W[k]*FV.W[k];
				}
			if ((int)config[2] == 2)					
				cv->U_phi[k] = FV.RHO[k] * FV.PHI[k];
		}
}


int cons2prim(struct i_f_var * ifv)
{
	const int dim = (int)config[0];
	const double eps = (int)config[4];
	
	ifv->RHO   = ifv->U_rho;	
	ifv->gamma = ifv->U_gamma / ifv->U_rho;
	ifv->U     = ifv->U_u/ifv->U_rho;
	ifv->P     = (ifv->U_e - 0.5*(ifv->U_u*ifv->U_u)/ifv->U_rho) * (ifv->gamma-1.0);
	if (dim > 1)
		{
			ifv->V  = ifv->U_v/ifv->U_rho;
			ifv->P -= (0.5*(ifv->U_v*ifv->U_v)/ifv->U_rho) * (ifv->gamma-1.0);
			if (isnan(ifv->V) || isinf(ifv->V))									
				return 0;
		}
	if (dim > 2)
		{			
			ifv->W  = ifv->U_w/ifv->U_rho;
			ifv->P -= (0.5*(ifv->U_w*ifv->U_w)/ifv->U_rho) * (ifv->gamma-1.0);
			if (isnan(ifv->W) || isinf(ifv->W))
				return 0;
		}
	if ((int)config[2] == 2)
		{
			ifv->PHI = ifv->U_phi/ifv->U_rho;
			if (isnan(ifv->PHI) || ifv->PHI < -0.1 || ifv->PHI > 1.0 + 0.1)//-100*eps
				return 0;
		}
	if (isnan(ifv->RHO + ifv->U + ifv->P) || isinf(ifv->RHO + ifv->U + ifv->P) || ifv->RHO < -100*eps || ifv->P < -100*eps)
		return 0;

	return 1;
}


int cons_qty_update
(struct cell_var * cv, const struct mesh_var mv,
 const struct flu_var FV, const double tau)
{
	const int dim = (int)config[0];
	const int num_cell = (int)config[3];
	int ** cp = mv.cell_pt;
	
	int p_p, p_n;
	double length, gamma, flux_v_fix, length2;
	int i;
	for(int k = 0; k < num_cell; ++k)
		{
			flux_v_fix = 0.0;				
			if(isinf(config[60]))
				gamma = cv->U_gamma[k]/cv->U_rho[k];
			for(int j = 0; j < cp[k][0]; j++)
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
					if (dim == 1)
						length = cv->n_x[k][j];
					else if (dim == 2)
						length = sqrt((mv.X[p_p] - mv.X[p_n])*(mv.X[p_p]-mv.X[p_n]) + (mv.Y[p_p] - mv.Y[p_n])*(mv.Y[p_p]-mv.Y[p_n]));
					
					cv->U_rho[k] += - tau*cv->F_rho[k][j] * length / cv->vol[k];
					cv->U_e[k]   += - tau*cv->F_e[k][j]   * length / cv->vol[k];	
					cv->U_u[k]   += - tau*cv->F_u[k][j]   * length / cv->vol[k];
					if (dim > 1)
						cv->U_v[k] += - tau*cv->F_v[k][j] * length / cv->vol[k];
					if (dim > 2)
						cv->U_w[k] += - tau*cv->F_w[k][j] * length / cv->vol[k];
					if ((int)config[2] == 2)
						cv->U_phi[k] += - tau*cv->F_phi[k][j] * length / cv->vol[k];
					if(!isinf(config[60]))
						cv->U_gamma[k] += - tau*cv->F_gamma[k][j] * length / cv->vol[k];
					if ((int)config[61] == 1)
						{													
							//							flux_v_fix += tau*length/cv->vol[k]*FV.RHO[k]*cv->RHO_p[k][j]*(cv->U_p[k][j]*cv->n_x[k][j]+cv->V_p[k][j]*cv->n_y[k][j])*((cv->U_p[k][j]-FV.U[k])*(cv->U_p[k][j]-FV.U[k])+(cv->V_p[k][j]-FV.V[k])*(cv->V_p[k][j]-FV.V[k]));
							for (i = j+1; i < cp[k][0]; i++)
								{
									if(i == cp[k][0]-1) 
										{
											p_p=cp[k][1];
											p_n=cp[k][i+1];
										}				  
									else
										{
											p_p=cp[k][i+2];
											p_n=cp[k][i+1];
										}
									if (dim == 1)
										length2 = cv->n_x[k][i];
									else if (dim == 2)
										length2 = sqrt((mv.X[p_p] - mv.X[p_n])*(mv.X[p_p]-mv.X[p_n]) + (mv.Y[p_p] - mv.Y[p_n])*(mv.Y[p_p]-mv.Y[p_n]));
						
									//									flux_v_fix -= tau*length/cv->vol[k]*cv->RHO_p[k][j]*(cv->U_p[k][j]*cv->n_x[k][j]+cv->V_p[k][j]*cv->n_y[k][j])*tau*length2/cv->vol[k]*cv->RHO_p[k][i]*(cv->U_p[k][i]*cv->n_x[k][i]+cv->V_p[k][i]*cv->n_y[k][i])*((cv->U_p[k][j]-cv->U_p[k][i])*(cv->U_p[k][j]-cv->U_p[k][i])+(cv->V_p[k][j]-cv->V_p[k][i])*(cv->V_p[k][j]-cv->V_p[k][i]));
								}

							flux_v_fix += tau*length/cv->vol[k]*cv->RHO_p[k][j]*(cv->U_p[k][j]*cv->n_x[k][j]+cv->V_p[k][j]*cv->n_y[k][j])*((cv->U_p[k][j]-FV.U[k])*(cv->U_p[k][j]-FV.U[k])+(cv->V_p[k][j]-FV.V[k])*(cv->V_p[k][j]-FV.V[k]))/2.0;

						}
				}
			if(isinf(config[60]))
				cv->U_gamma[k] = gamma*cv->U_rho[k];
			
			cv->U_e[k] += flux_v_fix;///cv->U_rho[k]/2.0;
		}
	return 0;
}


struct sub_cell_var {
	double length[8], vol, cos_a, sin_a;
	double F_rho[8], F_e[8], F_gamma[8], F_phi[8], F_u[8], F_v[8], F_w[8], F_delta_e[8];
	double U_rho, U_e, U_gamma, U_phi, U_u, U_v, U_w, delta_U_e;
	double   RHO,   P,   gamma,   PHI,   U,   V,   W; 
};


static void sub_cell_init
(const struct cell_var cv, const struct mesh_var mv,
 const double tau, struct sub_cell_var * scv, const int k)
{
	const int dim = (int)config[0];
	int ** cp = mv.cell_pt;
	
	int p_p, p, p_n, l_n, l_p; //(l_n)_j(line)

	for(int j = 0; j < cp[k][0]*2+2; j++)
		{
			scv[j].U_rho = cv.U_rho[k];
			scv[j].U_e   = cv.U_e[k];
			scv[j].U_u   = cv.U_u[k];
			scv[j].U_v   = cv.U_v[k];
			scv[j].delta_U_e = cv.delta_U_e[k];
		}
	
	for(int j = 0; j < cp[k][0]; j++)
		{
			if(j == cp[k][0]-1) 
				{
					p_p=cp[k][1];
					p  =cp[k][j+1];
				}
			else
				{
					p_p=cp[k][j+2];
					p  =cp[k][j+1];
				}
			
			if (dim == 1)
				scv[cp[k][0]*2+1].length[j] = cv.n_x[k][j];
			else if (dim == 2)
				scv[cp[k][0]*2+1].length[j] = sqrt((mv.X[p_p] - mv.X[p])*(mv.X[p_p]-mv.X[p]) + (mv.Y[p_p] - mv.Y[p])*(mv.Y[p_p]-mv.Y[p]));
			scv[cp[k][0]*2+1].vol      = cv.vol[k];
			scv[cp[k][0]*2+1].F_rho[j] = cv.F_rho[k][j];
			scv[cp[k][0]*2+1].F_e[j]   = cv.F_e[k][j];
			scv[cp[k][0]*2+1].F_u[j]   = cv.F_u[k][j];
			scv[cp[k][0]*2+1].F_v[j]   = cv.F_v[k][j];
			scv[cp[k][0]*2+1].F_delta_e[j] = cv.F_delta_e[k][j];

			scv[cp[k][0]*2].F_rho[j] = cv.F_rho_star[k][j];
			scv[cp[k][0]*2].F_e[j]   = cv.F_e_star[k][j];
			scv[cp[k][0]*2].F_u[j]   = cv.F_u_star[k][j];
			scv[cp[k][0]*2].F_v[j]   = cv.F_v_star[k][j];
		}

	for(int j = 0; j < cp[k][0]; j++)
		{
			if(j == cp[k][0]-1) 
				{
					p_p=cp[k][1];
					p  =cp[k][j+1];
					p_n=cp[k][j];
					l_n=j-1;
				}
			else if(j)
				{
					p_p=cp[k][j+2];
					p  =cp[k][j+1];
					p_n=cp[k][j];
					l_n=j-1;
				}
			else
				{
					p_p=cp[k][j+2];
					p  =cp[k][j+1];
					p_n=cp[k][cp[k][0]];
					l_n=cp[k][0]-1;
				}
			scv[j].cos_a  = (mv.X[p_p] - mv.X[p])*(mv.X[p_n]-mv.X[p])+(mv.Y[p_p] - mv.Y[p])*(mv.Y[p_n]-mv.Y[p]);
			scv[j].cos_a /= scv[cp[k][0]*2+1].length[j]*scv[cp[k][0]*2+1].length[l_n];
			scv[j].sin_a  = sqrt(1 - scv[j].cos_a*scv[j].cos_a);
			scv[j].length[0] = -cv.u_star[k][j]  *tau/scv[j].sin_a;
			scv[j].length[1] = -cv.u_star[k][l_n]*tau/scv[j].sin_a;
			scv[j].length[2] = scv[j].length[0];
			scv[j].length[3] = scv[j].length[1];
			scv[j].vol = scv[j].sin_a*scv[j].length[1]*scv[j].length[0];
			
			scv[j].F_rho[0] = cv.F_rho[k][l_n];
			scv[j].F_e[0]   = cv.F_e[k][l_n];
			scv[j].F_u[0]   = cv.F_u[k][l_n];
			scv[j].F_v[0]   = cv.F_v[k][l_n];
			scv[j].F_rho[1] = cv.F_rho[k][j];
			scv[j].F_e[1]   = cv.F_e[k][j];
			scv[j].F_u[1]   = cv.F_u[k][j];
			scv[j].F_v[1]   = cv.F_v[k][j];
			scv[j].F_rho[2] = -cv.F_rho_star[k][l_n];
			scv[j].F_e[2]   = -cv.F_e_star[k][l_n];
			scv[j].F_u[2]   = -cv.F_u_star[k][l_n];
			scv[j].F_v[2]   = -cv.F_v_star[k][l_n];
			scv[j].F_rho[3] = -cv.F_rho_star[k][j];
			scv[j].F_e[3]   = -cv.F_e_star[k][j];
			scv[j].F_u[3]   = -cv.F_u_star[k][j];
			scv[j].F_v[3]   = -cv.F_v_star[k][j];
		}

	for(int j = 0; j < cp[k][0]; j++)
		{
			if(j == cp[k][0]-1) 
				{
					l_p=0;
					l_n=j-1;
				}
			else if(j)
				{
					l_p=j+1;
					l_n=j-1;
				}
			else
				{
					l_p=j+1;
					l_n=cp[k][0]-1;
				}
			scv[j+cp[k][0]].length[0] = scv[j].length[2];
			scv[j+cp[k][0]].length[2] = scv[l_p].length[3];
			scv[j+cp[k][0]].length[1] = scv[cp[k][0]*2+1].length[j] - scv[j].length[1] - scv[l_p].length[0];
			scv[j+cp[k][0]].length[3] = scv[j+cp[k][0]].length[1] - scv[j+cp[k][0]].length[0]*scv[j].cos_a - scv[j+cp[k][0]].length[2]*scv[l_p].cos_a;	
			scv[j+cp[k][0]].vol = (scv[j+cp[k][0]].length[1] + scv[j+cp[k][0]].length[3]) * (-cv.u_star[k][j]*tau) / 2.0;
			
			scv[j].F_rho[0] = cv.F_rho_star[k][l_n];
			scv[j].F_e[0]   = cv.F_e_star[k][l_n];
			scv[j].F_u[0]   = cv.F_u_star[k][l_n];
			scv[j].F_v[0]   = cv.F_v_star[k][l_n];
			scv[j].F_rho[1] = cv.F_rho[k][j];
			scv[j].F_e[1]   = cv.F_e[k][j];
			scv[j].F_u[1]   = cv.F_u[k][j];
			scv[j].F_v[1]   = cv.F_v[k][j];
			scv[j].F_rho[2] = cv.F_rho_star[k][l_p];
			scv[j].F_e[2]   = cv.F_e_star[k][l_p];
			scv[j].F_u[2]   = cv.F_u_star[k][l_p];
			scv[j].F_v[2]   = cv.F_v_star[k][l_p];
			scv[j].F_rho[3] = -cv.F_rho_star[k][j];
			scv[j].F_e[3]   = -cv.F_e_star[k][j];
			scv[j].F_u[3]   = -cv.F_u_star[k][j];
			scv[j].F_v[3]   = -cv.F_v_star[k][j];
		}
	
	for(int j = 0; j < cp[k][0]; j++)
		scv[cp[k][0]*2].length[j] = scv[j+cp[k][0]].length[3];
	
	scv[cp[k][0]*2].vol = cv.vol[k];
	for(int j = 0; j < cp[k][0]*2; j++)
		scv[cp[k][0]*2].vol -= scv[j].vol;
}

static int sub_cell_update
(const struct cell_var * cv, const struct mesh_var mv,
 struct sub_cell_var * scv, const int k, double gamma)
{
	const int dim = (int)config[0];
	const double eps = (int)config[4];
	int ** cp = mv.cell_pt;

	for(int j = 0; j < cp[k][0]*2+2; j++)
		{
			scv[j].RHO   = scv[j].U_rho;	
			scv[j].gamma = gamma;
			scv[j].U     = scv[j].U_u/scv[j].U_rho;
			scv[j].P     = (scv[j].U_e - 0.5*(scv[j].U_u*scv[j].U_u)/scv[j].U_rho) * (scv[j].gamma-1.0);
			if (dim > 1)
				{
					scv[j].V  = scv[j].U_v/scv[j].U_rho;
					scv[j].P -= (0.5*(scv[j].U_v*scv[j].U_v)/scv[j].U_rho) * (scv[j].gamma-1.0);
					if (isnan(scv[j].V) || isinf(scv[j].V))									
						return 0;
				}
			if (dim > 2)
				{			
					scv[j].W  = scv[j].U_w/scv[j].U_rho;
					scv[j].P -= (0.5*(scv[j].U_w*scv[j].U_w)/scv[j].U_rho) * (scv[j].gamma-1.0);
					if (isnan(scv[j].W) || isinf(scv[j].W))
						return 0;
				}
			if ((int)config[2] == 2)
				{
					scv[j].PHI = scv[j].U_phi/scv[j].U_rho;
					if (isnan(scv[j].PHI) || scv[j].PHI < -0.1 || scv[j].PHI > 1.0 + 0.1)//-100*eps
						return 0;
				}
			if (isnan(scv[j].RHO + scv[j].U + scv[j].P) || isinf(scv[j].RHO + scv[j].U + scv[j].P) || scv[j].RHO < -100*eps || scv[j].P < -100*eps)
				return 0;
		}

	double P_ave = 0.0, delta_U_e;
	for(int j = 0; j < cp[k][0]*2+1; j++)
		P_ave += scv[j].vol*scv[j].P;
	P_ave /= cv->vol[k];
	delta_U_e = scv[cp[k][0]*2+1].U_e - (scv[cp[k][0]*2+1].P/(scv[cp[k][0]*2+1].gamma-1.0) + (0.5*(scv[cp[k][0]*2+1].U_u*scv[cp[k][0]*2+1].U_u + scv[cp[k][0]*2+1].U_v*scv[cp[k][0]*2+1].U_v)/scv[cp[k][0]*2+1].U_rho));
										 
	cv->U_rho[k] = scv[cp[k][0]*2+1].U_rho;
	cv->U_e[k]   = scv[cp[k][0]*2+1].U_e - delta_U_e;
	cv->U_u[k]   = scv[cp[k][0]*2+1].U_u;
	cv->U_v[k]   = scv[cp[k][0]*2+1].U_v;
	cv->delta_U_e[k] += delta_U_e;

	return 1;	
}

int cons_qty_update_corr_ave_P
(struct cell_var * cv, const struct mesh_var mv,
 const struct flu_var FV, const double tau)
{
	const int dim = (int)config[0];
	const int num_cell = (int)config[3];
	int ** cp = mv.cell_pt;
	
	double gamma;

	struct sub_cell_var scv[18];
	int i, j, num_border;
	for(int k = 0; k < num_cell; ++k)
		{
			if(isinf(config[60]))
				gamma = cv->U_gamma[k]/cv->U_rho[k];
			
			sub_cell_init(*cv, mv, tau, scv, k);
			for(i = 0; i < cp[k][0]*2+2; i++)
				{					
					if (i == cp[k][0]*2 || i == cp[k][0]*2+1)
						num_border = cp[k][0];					
					else
						num_border = 4;
					for(j = 0; j < num_border; j++)
						{					
							scv[i].U_rho += -tau*scv[i].F_rho[j] * scv[i].length[j] / scv[i].vol;
							scv[i].U_e   += -tau*scv[i].F_e[j]   * scv[i].length[j] / scv[i].vol; 
							scv[i].U_u   += -tau*scv[i].F_u[j]   * scv[i].length[j] / scv[i].vol;
							if (dim > 1)
								scv[i].U_v += -tau*scv[i].F_v[j] * scv[i].length[j] / scv[i].vol;
							if ((int)config[2] == 2)
								scv[i].U_phi += -tau*scv[i].F_phi[j] * scv[i].length[j] / scv[i].vol;
							if (i == cp[k][0]*2+1)
								scv[i].delta_U_e += -tau*scv[i].F_delta_e[j] * scv[i].length[j] / scv[i].vol;
						}
				}			
			if(sub_cell_update(cv, mv, scv, k, gamma) == 0)
				{
					fprintf(stderr, "Error happens on sub cell updating primitive variable!\n");
					return 0;
				}
			
			if(isinf(config[60]))
				cv->U_gamma[k] = gamma*cv->U_rho[k];
		}
	return 1;
}
