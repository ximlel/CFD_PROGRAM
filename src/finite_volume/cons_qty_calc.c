#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"


//Initialize conserved quantities.
void cons_qty_init(struct cell_var * cv, const struct flu_var FV)
{
	const int dim = (int)config[0];
	for(int k = 0; k < (int)config[3]; k++)
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

			cv->U0_rho[k]   = cv->U_rho[k];
			cv->U0_gamma[k] = cv->U_gamma[k];
			cv->U0_e[k]     = cv->U_e[k];
			cv->U0_u[k]     = cv->U_u[k];			
			if (dim > 1)									
				cv->U0_v[k] = cv->U_v[k];
			if (dim > 2)									
				cv->U0_w[k] = cv->U_w[k];
			if ((int)config[2] == 2)					
				cv->U0_phi[k] = cv->U_phi[k];
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
			if (isnan(ifv->PHI) || ifv->PHI < -10*eps || ifv->PHI > 1.0 + 10*eps)
				return 0;
		}

	if (isnan(ifv->RHO + ifv->U + ifv->P) || isinf(ifv->RHO + ifv->U + ifv->P) || ifv->RHO < -10*eps || ifv->P < -10*eps)
		return 0;

	return 1;
}


void cons_qty_update(struct cell_var * cv, const struct mesh_var mv, const double tau)
{
	const int dim = (int)config[0];
	const int num_cell = (int)config[3];
	int ** cp = mv.cell_pt;
	
	int p_p, p_n;
	double length;
	for(int k = 0; k < num_cell; ++k)
		{
			for(int j = 0; j < cp[k][0]; ++j)
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
				}
			if(isinf(config[60]))
				cv->U_gamma[k] = cv->U0_gamma[k]/cv->U0_rho[k]*cv->U_rho[k];
		}
}