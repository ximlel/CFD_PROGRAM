#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/Riemann_solver.h"




void linear_GRP_solver_Edir
(double * direvative, double * mid,
 double rho_L, double rho_R, double s_rho_L, double s_rho_R,
 double   u_L, double   u_R, double   s_u_L, double   s_u_R,
 double   p_L, double   p_R, double   s_p_L, double   s_p_R,
 double gamma, double eps)
{
  double dist;
  double c_L, c_R;
  int CRW[2];
  double u_star, p_star, rho_star_L, rho_star_R, c_star_L, c_star_R;

  double PI, H1, H2, H3;
  double a_L, b_L, d_L, a_R, b_R, d_R;
  double L_u, L_p, L_rho;

  double u_t_mat, p_t_mat;

  double shk_spd, zeta = (gamma-1.0)/(gamma+1.0), zts = zeta*zeta;

  double g_rho, g_u, g_p, f;

  double speed_L, speed_R;

  //printf("rho_L=%lf\tu_L=%lf\tp_L=%lf\n", rho_L, u_L, p_L);
  //printf("rho_R=%lf\tu_R=%lf\tp_R=%lf\n", rho_R, u_R, p_R);


  c_L = sqrt(gamma * p_L / rho_L);
  c_R = sqrt(gamma * p_R / rho_R);

  //dist = (rho_L-rho_R)*(rho_L-rho_R) + (u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R);
  dist = (u_L-u_R)*(u_L-u_R) + (p_L-p_R)*(p_L-p_R);
  dist = sqrt(dist);
//=========acoustic case==========
  if(dist < eps)
  {
    //printf("AC\n");
    //------trivial case------
    if(u_L-c_L > 0.0) //the t-axe is on the left side of all the three waves
    {
      direvative[0] = -s_rho_L*u_L - rho_L*s_u_L;
      direvative[1] = (direvative[0]*u_L + s_rho_L*u_L*u_L + 2.0*rho_L*u_L*s_u_L + s_p_L) / -rho_L;
      direvative[2] = -(gamma-1.0) * (0.5*direvative[0]*u_L*u_L + rho_L*u_L*direvative[1]);
      direvative[2] = direvative[2] - s_u_L * (gamma*p_L + 0.5*(gamma-1.0)*rho_L*u_L*u_L);
      direvative[2] = direvative[2] - u_L * (gamma * s_p_L + (gamma-1.0)*(0.5*s_rho_L*u_L*u_L + rho_L*u_L*s_u_L));

      mid[0] = rho_L;
      mid[1] =   u_L;
      mid[2] =   p_L;
    }
    else if(u_R+c_R < 0.0) //the t-axe is on the right side of all the three waves
    {
      direvative[0] = -s_rho_R*u_R - rho_R*s_u_R;
      direvative[1] = (direvative[0]*u_R + s_rho_R*u_R*u_R + 2.0*rho_R*u_R*s_u_R + s_p_R) / -rho_R;
      direvative[2] = -(gamma-1.0) * (0.5*direvative[0]*u_R*u_R + rho_R*u_R*direvative[1]);
      direvative[2] = direvative[2] - s_u_R * (gamma*p_R + 0.5*(gamma-1.0)*rho_R*u_R*u_R);
      direvative[2] = direvative[2] - u_R * (gamma * s_p_R + (gamma-1.0)*(0.5*s_rho_R*u_R*u_R + rho_R*u_R*s_u_R));

      mid[0] = rho_R;
      mid[1] =   u_R;
      mid[2] =   p_R;
    }
    //------non-trivial case------
    else
    {
      //printf("%.9lf\t%.9lf\t%.9lf\n", lambda, A, gamma);
      //printf("%.9lf\t%.9lf\t%.9lf\t%.9lf\n", t_rho_L, rho_L, rho_R, t_rho_R);
      //printf("%.9lf\t%.9lf\t%.9lf\t%.9lf\n", t_u_L, u_L, u_R, t_u_R);
      //printf("%.9lf\t%.9lf\t%.9lf\t%.9lf\n", t_p_L, p_L, p_R, t_p_R);

      //printf("AC\n");
      rho_star_L = rho_L;
      rho_star_R = rho_R;
      c_star_L =   c_L;
      c_star_R =   c_R;
      u_star = 0.5*(u_R+u_L);
      p_star = 0.5*(p_R+p_L);

      if(u_star > 0.0)
      {
	mid[0] = rho_star_L;
	mid[1] =   u_star;
	mid[2] =   p_star;

	/*
	direvative[1] = -0.5*((u_star+c_star_L)*(s_u_L+s_p_L/rho_star_L/c_star_L) + (u_star-c_star_R)*(s_u_R-s_p_R/rho_star_R/c_star_R));
	direvative[2] = -0.5*rho_star_L*c_star_L*((u_star+c_star_L)*(s_u_L+s_p_L/rho_star_L/c_star_L) - (u_star-c_star_R)*(s_u_R-s_p_R/rho_star_R/c_star_R));

	direvative[0] = (u_star*(s_p_L - s_rho_L*c_star_L*c_star_L) + direvative[2])/c_star_L/c_star_L;
	*/
	PI = (u_star+c_star_R)*rho_star_L*c_star_L*c_star_L / (u_star-c_star_L)/rho_star_R/c_star_R/c_star_R;
	direvative[1] = (s_p_L/rho_L+c_L*s_u_L)*PI/(1.0-PI) + (s_p_R/rho_R-c_R*s_u_R)/(PI-1.0);
	direvative[2] = ((u_star+c_star_R)/rho_star_R/c_star_R/c_star_R) - ((u_star-c_star_L)/rho_star_L/c_star_L/c_star_L);
	direvative[2] = (s_p_R/rho_R-c_R*s_u_R-s_p_L/rho_L-c_L*s_u_L) / direvative[2];
	direvative[2] = direvative[2] * (1.0 - (u_star*u_star/c_star_L/c_star_L)) + rho_star_L*u_star*direvative[1];
	direvative[0] = (u_star*(s_p_L - s_rho_L*c_star_L*c_star_L) + direvative[2])/c_star_L/c_star_L;
      }
      else
      {
	mid[0] = rho_star_R;
	mid[1] =   u_star;
	mid[2] =   p_star;

	/*
	direvative[1] = -0.5*((u_star+c_star_L)*(s_u_L+s_p_L/rho_star_L/c_star_L) + (u_star-c_star_R)*(s_u_R-s_p_R/rho_star_R/c_star_R));
	direvative[2] = -0.5*rho_star_L*c_star_L*((u_star+c_star_L)*(s_u_L+s_p_L/rho_star_L/c_star_L) - (u_star-c_star_R)*(s_u_R-s_p_R/rho_star_R/c_star_R));

	direvative[0] = (u_star*(s_p_R - s_rho_R*c_star_R*c_star_R) + direvative[2])/c_star_R/c_star_R;
	*/
	PI = (u_star+c_star_R)*rho_star_L*c_star_L*c_star_L / (u_star-c_star_L)/rho_star_R/c_star_R/c_star_R;
	direvative[1] = (s_p_L/rho_L+c_L*s_u_L)*PI/(1.0-PI) + (s_p_R/rho_R-c_R*s_u_R)/(PI-1.0);
	direvative[2] = ((u_star+c_star_R)/rho_star_R/c_star_R/c_star_R) - ((u_star-c_star_L)/rho_star_L/c_star_L/c_star_L);
	direvative[2] = (s_p_R/rho_R-c_R*s_u_R-s_p_L/rho_L-c_L*s_u_L) / direvative[2];
	direvative[2] = direvative[2] * (1.0 - (u_star*u_star/c_star_R/c_star_R)) + rho_star_R*u_star*direvative[1];
	direvative[0] = (u_star*(s_p_R - s_rho_R*c_star_R*c_star_R) + direvative[2])/c_star_R/c_star_R;
      }
      //printf("%lf\t%lf\t%lf\n\n", direvative[0], direvative[1], direvative[2]);
    }

    //printf("rho=%lf\tu=%lf\tt=%lf\n", mid[0], mid[1], mid[2]);
    return;
  }

//=========non-acoustic case==========
  //printf("u_L=%lf\tp_L=%lf\n", u_L, p_L);
  //printf("u_R=%lf\tp_R=%lf\n", u_R, p_R);
  Riemann_solver_exact(&u_star, &p_star, gamma, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, 50);
  //CRW[0] = 1;
  //CRW[1] = 1;
  if(p_star > p_L)
    rho_star_L = rho_L*(p_star+zeta*p_L)/(p_L+zeta*p_star);
  else
    rho_star_L = rho_L*pow(p_star/p_L,1.0/gamma);
  if(p_star > p_R)
    rho_star_R = rho_R*(p_star+zeta*p_R)/(p_R+zeta*p_star);
  else
    rho_star_R = rho_R*pow(p_star/p_R,1.0/gamma);
  c_star_L = sqrt(gamma * p_star / rho_star_L);
  c_star_R = sqrt(gamma * p_star / rho_star_R);



  //printf("%d\t%d\n", CRW[0], CRW[1]);
  //printf("u*=%lf\tp*=%lf\n", u_star, p_star);
  //printf("rho*L=%lf\trho*R=%lf\n\n", rho_star_L, rho_star_R);
  /*
  if(CRW[0])
    printf("speed_L:%lf\t%lf\n", u_L-c_L, u_star-c_star_L);
  else
    printf("speed_L:%lf\t%lf\t%lf\n", (rho_star_L*u_star - rho_L*u_L) / (rho_star_L - rho_L), u_L-c_L, u_star-c_star_L);
  if(CRW[1])
    printf("speed_R:%lf\t%lf\n\n", u_star+c_star_R, u_R+c_R);
  else
    printf("speed_R:%lf\n\n", (rho_star_R*u_star - rho_R*u_R) / (rho_star_R - rho_R));
  //*/
//----------solving the LINEAR GRP----------
  if(CRW[0])
    speed_L = u_L - c_L;
  else
    speed_L = (rho_star_L*u_star - rho_L*u_L) / (rho_star_L - rho_L);
  if(CRW[1])
    speed_R = u_R + c_R;
  else
    speed_R = (rho_star_R*u_star - rho_R*u_R) / (rho_star_R - rho_R);

  //------trivial case------
  if(speed_L > 0.0) //the t-axe is on the left side of all the three waves
  {
    direvative[0] = -s_rho_L*u_L - rho_L*s_u_L;
    direvative[1] = (direvative[0]*u_L + s_rho_L*u_L*u_L + 2.0*rho_L*u_L*s_u_L + s_p_L) / -rho_L;
    direvative[2] = (s_u_L*p_L + u_L*s_p_L)*gamma/(1.0-gamma) - 0.5*s_rho_L*u_L*u_L*u_L - 1.5*rho_L*u_L*u_L*s_u_L;
    direvative[2] = direvative[2] - 0.5*direvative[0]*u_L*u_L - rho_L*u_L*direvative[1];
    direvative[2] = direvative[2] * (gamma-1.0);

    mid[0] = rho_L;
    mid[1] =   u_L;
    mid[2] =   p_L;
  }
  else if(speed_R < 0.0) //the t-axe is on the right side of all the three waves
  {
    direvative[0] = -s_rho_R*u_R - rho_R*s_u_R;
    direvative[1] = (direvative[0]*u_R + s_rho_R*u_R*u_R + 2.0*rho_R*u_R*s_u_R + s_p_R) / -rho_R;
    direvative[2] = -(gamma-1.0) * (0.5*direvative[0]*u_R*u_R + rho_R*u_R*direvative[1]);
    direvative[2] = direvative[2] - s_u_R * (gamma*p_R + 0.5*(gamma-1.0)*rho_R*u_R*u_R);
    direvative[2] = direvative[2] - u_R * (gamma * s_p_R + (gamma-1.0)*(0.5*s_rho_R*u_R*u_R + rho_R*u_R*s_u_R));

    mid[0] = rho_R;
    mid[1] =   u_R;
    mid[2] =   p_R;
  }
  //----non-trivial case----
  else
  {
    if((CRW[0]) && ((u_star-c_star_L) > 0.0)) // the t-axe is in a 1-CRW
    {
      shk_spd = (rho_star_L*u_star - rho_L*u_L)/(rho_star_L - rho_L);

      mid[1] = zeta*(u_L+2.0*c_L/(gamma-1.0));
      mid[2] = mid[1]*mid[1]*rho_L/gamma/pow(p_L, 1.0/gamma);
      mid[2] = pow(mid[2], gamma/(gamma-1.0));
      mid[0] = gamma*mid[2]/mid[1]/mid[1];

      direvative[1] = 0.5*(pow(mid[1]/c_L, 0.5/zeta)*(1.0+zeta) + pow(mid[1]/c_L, (1.0+zeta)/zeta)*zeta)/(0.5+zeta);
      direvative[1] = direvative[1] * (s_p_L - s_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
      //printf("ddd%lf\n", c_L*pow(mid[1]/c_L, 0.5/zeta)*(s_u_L + (gamma*s_p_L/c_L - c_L*s_rho_L)/(gamma-1.0)/rho_L));
      direvative[1] = direvative[1] - c_L*pow(mid[1]/c_L, 0.5/zeta)*(s_u_L + (gamma*s_p_L/c_L - c_L*s_rho_L)/(gamma-1.0)/rho_L);

      direvative[2] = mid[0]*mid[1]*direvative[1];

      direvative[0] = mid[0]*mid[1]*pow(mid[1]/c_L, (1.0+zeta)/zeta)*(s_p_L - s_rho_L*c_L*c_L)/rho_L;
      direvative[0] = (direvative[0] + direvative[2]) / mid[1]/mid[1];
    }
    else if((CRW[1]) && ((u_star+c_star_R) < 0.0)) // the t-axe is in a 3-CRW
    {
      shk_spd = (rho_star_R*u_star - rho_R*u_R)/(rho_star_R - rho_R);

      mid[1] = zeta*(u_R-2.0*c_R/(gamma-1.0));
      //c_sonic = -mid[1]
      mid[2] = mid[1]*mid[1]*rho_R/gamma/pow(p_R, 1.0/gamma);
      mid[2] = pow(mid[2], gamma/(gamma-1.0));
      mid[0] = gamma*mid[2]/mid[1]/mid[1];

      direvative[1] = 0.5*(pow(-mid[1]/c_R, 0.5/zeta)*(1.0+zeta) + pow(-mid[1]/c_R, (1.0+zeta)/zeta)*zeta)/(0.5+zeta);
      direvative[1] = direvative[1] * (s_p_R - s_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
      //printf("ddd%lf\n", c_R*pow(-mid[1]/c_R, 0.5/zeta)*(s_u_R - (gamma*s_p_R/c_R - c_R*s_rho_R)/(gamma-1.0)/rho_R));
      direvative[1] = direvative[1] + c_R*pow(-mid[1]/c_R, 0.5/zeta)*(s_u_R - (gamma*s_p_R/c_R - c_R*s_rho_R)/(gamma-1.0)/rho_R);

      direvative[2] = mid[0]*mid[1]*direvative[1];

      direvative[0] = mid[0]*mid[1]*pow(-mid[1]/c_R, (1.0+zeta)/zeta)*(s_p_R - s_rho_R*c_R*c_R)/rho_R;
      direvative[0] = (direvative[0] + direvative[2]) / mid[1]/mid[1];
    }
    //--non-sonic case--
    else
    {
      //printf("CRW\n");
    //determine a_L, b_L and d_L
      if(CRW[0]) //the 1-wave is a CRW
      {
	a_L = 1.0;
        b_L = 1.0 / rho_star_L / c_star_L;
	d_L = 0.5*(pow(c_star_L/c_L, 0.5/zeta)*(1.0+zeta) + pow(c_star_L/c_L, (1.0+zeta)/zeta)*zeta)/(0.5+zeta);
	d_L = d_L * (s_p_L - s_rho_L*c_L*c_L)/(gamma-1.0)/rho_L;
	d_L = d_L - c_L*pow(c_star_L/c_L, 0.5/zeta)*(s_u_L + (gamma*s_p_L/c_L - c_L*s_rho_L)/(gamma-1.0)/rho_L);
      }
      else //the 1-wave is a shock
      {
	H1 = 0.5*sqrt((1.0-zeta)/(rho_L*(p_star+zeta*p_L))) * (p_star + (1.0+2.0*zeta)*p_L)/(p_star+zeta*p_L);
	H2 = -0.5*sqrt((1.0-zeta)/(rho_L*(p_star+zeta*p_L))) * ((2.0+zeta)*p_star + zeta*p_L)/(p_star+zeta*p_L);
	H3 = -0.5*sqrt((1.0-zeta)/(rho_L*(p_star+zeta*p_L))) * (p_star-p_L) / rho_L;
	shk_spd = (rho_star_L*u_star - rho_L*u_L)/(rho_star_L - rho_L);

	a_L = 1.0 - rho_star_L*(shk_spd-u_star)*H1;
	b_L = (u_star - shk_spd)/rho_star_L/c_star_L/c_star_L + H1;

	L_rho = (u_L-shk_spd) * H3;
	L_u = shk_spd - u_L + rho_L*c_L*c_L*H2 + rho_L*H3;
	L_p = (u_L-shk_spd)*H2 - 1.0/rho_L;

	d_L = L_rho*s_rho_L + L_u*s_u_L + L_p*s_p_L;
      }
    //determine a_R, b_R and d_R
      if(CRW[1]) //the 3-wave is a CRW
      {
	a_R = 1.0;
        b_R = -1.0 / rho_star_R / c_star_R;
	d_R = 0.5*(pow(c_star_R/c_R, 0.5/zeta)*(1.0+zeta) + pow(c_star_R/c_R, (1.0+zeta)/zeta)*zeta)/(0.5+zeta);
	d_R = d_R * (s_p_R - s_rho_R*c_R*c_R)/(gamma-1.0)/rho_R;
	d_R = d_R + c_R*pow(c_star_R/c_R, 0.5/zeta)*(s_u_R - (gamma*s_p_R/c_R - c_R*s_rho_R)/(gamma-1.0)/rho_R);
      }
      else //the 3-wave is a shock
      {
	H1 = 0.5*sqrt((1.0-zeta)/(rho_R*(p_star+zeta*p_R))) * (p_star + (1.0+2.0*zeta)*p_R)/(p_star+zeta*p_R);
	H2 = -0.5*sqrt((1.0-zeta)/(rho_R*(p_star+zeta*p_R))) * ((2.0+zeta)*p_star + zeta*p_R)/(p_star+zeta*p_R);
	H3 = -0.5*sqrt((1.0-zeta)/(rho_R*(p_star+zeta*p_R))) * (p_star-p_R) / rho_R;
	shk_spd = (rho_star_R*u_star - rho_R*u_R)/(rho_star_R - rho_R);

	a_R = 1.0 + rho_star_R*(shk_spd-u_star)*H1;
	b_R = (u_star - shk_spd)/rho_star_R/c_star_R/c_star_R - H1;

	L_rho = (shk_spd-u_R) * H3;
	L_u = shk_spd - u_R - rho_R*c_R*c_R*H2 - rho_R*H3;
	L_p = (shk_spd-u_R)*H2 - 1.0/rho_R;

	d_R = L_rho*s_rho_R + L_u*s_u_R + L_p*s_p_R;
      }

      p_t_mat = (d_L*a_R/a_L-d_R)/(b_L*a_R/a_L-b_R);
      u_t_mat = (d_L - b_L*p_t_mat)/a_L;

      //printf("a_L=%lf\tb_L=%lf\td_L=%lf\n", a_L, b_L, d_L);
      //printf("a_R=%lf\tb_R=%lf\td_R=%lf\n", a_R, b_R, d_R);
      //printf("p_t_Lag=%lf\tu_t_Lag=%lf\t", p_t_Lag, u_t_Lag);

      if(u_star < 0.0) //the t-axi is between the contact discontinuety and the 3-wave
      {
	mid[0] = rho_star_R;
	mid[1] =   u_star;
	mid[2] =   p_star;
        direvative[1] = u_t_mat + u_star*p_t_mat/rho_star_R/c_star_R/c_star_R;
        direvative[2] = p_t_mat + rho_star_R*u_star * u_t_mat;

	if(CRW[1]) //the 3-wave is a CRW
	{
	  direvative[0] = rho_star_R*u_star*pow(c_star_R/c_R, (1.0+zeta)/zeta)*(s_p_R - s_rho_R*c_R*c_R)/rho_R;
	  direvative[0] = (direvative[0] + direvative[2]) / c_star_R/c_star_R;
	}
	else //the 3-wave is a shock
	{
	  shk_spd = (rho_star_R*u_star - rho_R*u_R)/(rho_star_R - rho_R);
	  H1 = rho_R * p_R    * (1.0 - zts) / (p_R + zeta*p_star) / (p_R + zeta*p_star);
	  H2 = rho_R * p_star * (zts - 1.0) / (p_R + zeta*p_star) / (p_R + zeta*p_star);
	  H3 = (p_star + zeta*p_R) / (p_R + zeta*p_star);

	  g_rho = u_star-shk_spd;
	  g_u   = u_star*rho_star_R*(shk_spd-u_star)*H1;
	  g_p   = shk_spd/c_star_R/c_star_R - u_star*H1;
	  f = (shk_spd-u_R)*(H2*s_p_R + H3*s_rho_R) - rho_R*(H2*c_R*c_R+H3)*s_u_R;

	  direvative[0] = (f*u_star - g_p*p_t_mat - g_u*u_t_mat) / g_rho;
	}
      }
      else //the t-axi is between the 1-wave and the contact discontinuety
      {
	mid[0] = rho_star_L;
	mid[1] =   u_star;
	mid[2] =   p_star;
        direvative[1] = u_t_mat + u_star*p_t_mat/rho_star_L/c_star_L/c_star_L;
        direvative[2] = p_t_mat + rho_star_L*u_star * u_t_mat;
	if(CRW[0]) //the 1-wave is a CRW
	{
	  direvative[0] = rho_star_L*u_star*pow(c_star_L/c_L, (1.0+zeta)/zeta)*(s_p_L - s_rho_L*c_L*c_L)/rho_L;
	  direvative[0] = (direvative[0] + direvative[2]) / c_star_L/c_star_L;
	}
	else //the 1-wave is a shock
	{
	  shk_spd = (rho_star_L*u_star - rho_L*u_L)/(rho_star_L - rho_L);
	  H1 = rho_L * p_L    * (1.0 - zts) / (p_L + zeta*p_star) / (p_L + zeta*p_star);
	  H2 = rho_L * p_star * (zts - 1.0) / (p_L + zeta*p_star) / (p_L + zeta*p_star);
	  H3 = (p_star + zeta*p_L) / (p_L + zeta*p_star);

	  g_rho = u_star-shk_spd;
	  g_u   = u_star*rho_star_L*(shk_spd-u_star)*H1;
	  g_p   = shk_spd/c_star_L/c_star_L - u_star*H1;
	  f = (shk_spd-u_L)*(H2*s_p_L + H3*s_rho_L) - rho_L*(H2*c_L*c_L+H3)*s_u_L;

	  direvative[0] = (f*u_star - g_p*p_t_mat - g_u*u_t_mat) / g_rho;
	}
      }
    //--end of non-sonic case--
    }
  //----end of non-trivial case----
  }
}
