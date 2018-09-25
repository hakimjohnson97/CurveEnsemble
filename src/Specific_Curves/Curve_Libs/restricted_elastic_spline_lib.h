/*
 * Copyright (C) 2014, 2015, Hakim Johnson, Michael J. Johnson
 *
 * This file is part of Curve Ensemble.
 *
 * Curve Ensemble is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Curve Ensemble is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Curve Ensemble.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "common_curve_lib.h"
#include <fstream>
//#include <QDebug>

namespace ResElasticS
{

// Macros:
//
// If t belongs to [0,pi], then the following performs t=tau(t).  After exiting, it
// leaves c8=cos(t) and s8=sin(t)
#define TAU(t) cos_sin(t); 2*acos(c8*M_1_SQRT2)-M_PI_2;
struct piece { //FORM8 = [Form, flag_reversal, flag_reflection, flag_free] and flag_free
  short form;
  bool flag_reversal, flag_reflection;
  double p1,p2;
  complex c1,c2;
};

bool is_closed (double **H, int n){
// We'll make an educated guess as to whether of not H is closed.
// It is likely (but not guaranteed) that
if (H[0][STATE]<=0 || H[0][STATE]>100 || H[n-1][STATE]<=0 || H[n-1][STATE]>100)
  return true;
else
  return false;
}// is_closed----------------------------------------------------------------------------------------------

inline double tau(double t){
  //del=tau(t) returns the turning angle in E_[0,t] assuming t in [0,pi]
  cos_sin(t); s8=s8*s8; c8=c8*sqrt(1+s8);
  return cos_sin_inv ();
  /*if (t>=0 && t<0.7) return asin(s8*s8);
    //{t=sin(t); t=asin(t*t);}
  else if (t>2.4 && t<=M_PI) return M_PI-asin(s8*s8);
    //{t=sin(t); t=M_PI-asin(t*t);}
  else
    return 2*acos(c8*M_1_SQRT2)-M_PI_2; */
}// tau--------------------------------------------------------------------------------------------------------

inline double tau_inv(double t){
  // r=tau_inv(t) returns r in [0,pi] such that the turning angle of E_[0,r] equals r, where
  // t is assumed to be in [0,pi].  After exiting, I leave:
  // c8=cos(r), s8=sin(r), c_wig8=cos(t), s_wig8=sin(t)
  cos_sin(t); c_wig8=c8; s_wig8=s8;
  c8=c8/sqrt(1+s8); s8=sqrt(s8);
  return cos_sin_inv ();
}// tau_inv--------------------------------------------------------------------------------------------------------

double xi(double t,bool use_current_c8_s8=false){
  //y=xi(t) returns the function y=xi(t), employed in Bernoulli's rectangular elastica, defined by
  //dy/dt = sin^2(t)/sqrt(2-cos^2(t)), xi(0)=0.   (see Notebook April 16, 2010)
  // I call cos_sin(t) after running MOD2PI(t), so I require t in [-3pi,3pi]
  // After exiting, I lease and c8 = cos(t) and s8 = sin(t).
  // If c8 and s8 already hold cos(t) and sin(t), then the
  // y=xi(t,true)  will compute xi(t) using current values of c8 and s8.
  double c2,y;
  temp8=M_D_PI*t;
  if (!use_current_c8_s8) {MOD2PI(t); cos_sin(t);}
  c2=c8*c8;
   y=c2*1.00035293417090e-18 +  2.06221188091345e-18;
   y=c2*y+4.25387245342757e-18;
   y=c2*y+8.78047912720486e-18;
   y=c2*y+1.81362169223682e-17;
   y=c2*y+3.74871904367177e-17;
   y=c2*y+7.75426486345260e-17;
   y=c2*y+1.60522057722184e-16;
   y=c2*y+3.32568415385792e-16;
   y=c2*y+6.89598765055463e-16;
   y=c2*y+1.43119693457556e-15;
   y=c2*y+2.97309930362496e-15;
   y=c2*y+6.18226633414578e-15;
   y=c2*y+1.28687737418066e-14;
   y=c2*y+2.68165278917706e-14;
   y=c2*y+5.59462164114658e-14;
   y=c2*y+1.16861350359742e-13;
   y=c2*y+2.44418976348451e-13;
   y=c2*y+5.11914662197353e-13;
   y=c2*y+1.07373462108851e-12;
   y=c2*y+2.25566687632569e-12;
   y=c2*y+4.74654652205888e-12;
   y=c2*y+1.00059040759765e-11;
   y=c2*y+2.11333029715903e-11;
   y=c2*y+4.47274316493846e-11;
   y=c2*y+9.48738740321125e-11;
   y=c2*y+2.01727290148972e-10;
   y=c2*y+4.30048528481521e-10;
   y=c2*y+9.19404312329470e-10;
   y=c2*y+1.97173077912547e-09;
   y=c2*y+4.24300512071817e-09;
   y=c2*y+9.16511379379223e-09;
   y=c2*y+1.98801299577944e-08;
   y=c2*y+4.33238794956947e-08;
   y=c2*y+9.49097007045419e-08;
   y=c2*y+2.09154027544043e-07;
   y=c2*y+4.64038517339675e-07;
   y=c2*y+1.03756513778596e-06;
   y=c2*y+2.34099224789788e-06;
   y=c2*y+5.33832157495902e-06;
   y=c2*y+1.23291087891121e-05;
   y=c2*y+2.89186408612988e-05;
   y=c2*y+6.91481766888737e-05;
   y=c2*y+1.69459629368445e-04;
   y=c2*y+4.29044510876595e-04;
   y=c2*y+1.13658790424165e-03;
   y=c2*y+3.22078817165385e-03;
   y=c2*y+1.02025698829527e-02;
   y=c2*y+4.03745709937904e-02;
   y=c2*y+3.25726899435641e-01;
  return temp8-s8*c8*y;
  //y=0.381379881750907*t-0.5*sin(2*t).*polyval(E_coeff,cos(t).^2);
  return y;
}// xi--------------------------------------------------------------------------------------------------------

inline double xi_tau_inv(double t){
  // y=xi_tau_inv(t) returns y=xi(tau_inv(t)).
  // After exiting, I leave c8=cos(r), s8=sin(r), c_wig8=cos(t), s_wig8=sin(t),
  // where r denotes tau_inv(t).
  cos_sin(t); c_wig8=c8; s_wig8=s8;
  c8=c8/sqrt(1+s8); s8=sqrt(s8);
  return xi(cos_sin_inv (),true);
}//xi_tau_inv-------------------------------------------------------------------------------------------------

inline double G_val(double gamma){
  //v=G_val(gamma) returns G(gamma), where alpha,beta are global variables and gamma is a scalar.
  //I assume alpha and beta are canonical and gamma belongs to Gamma.
  //temp8=xi(tau_inv(alpha-gamma))+xi(tau_inv(beta-gamma));
  temp8=xi_tau_inv(alpha-gamma)+xi_tau_inv(beta-gamma);
  return -temp8*temp8/sin(gamma);
}// G_val--------------------------------------------------------------------------------------------------------

double sigma_wig(double gamma){
  //w=sigma_wig(gamma) returns sigma_wig(gamma), and if the global variable nargout==2,
  //we set df_dx=sigma_wig'(gamma). 
  //We use the global variables alpha,beta, which are assumed to be in canonical form,
  //and it is assumed that gamma belongs to the interval Gamma.
  //See page 10,11 of Notebook Dec. 24, 2011
  //After exiting, I leave c8=cos(gamma) and s8=sin(gamma).
  double y,ta,tb,ta_plus_tb,ca,cb;
  // The following four lines of optimized code result in:
  // ta=sqrt(sin(alpha-gamma)) and ca=cos(alpha-gamma))
  // tb=sqrt(sin(beta-gamma)) and cb=cos(beta-gamma))
  // y=xi(tau_inv(alpha-gamma))+xi(tau_inv(beta-gamma));
  // s8=sin(gamma) and c8=cos(gamma)
  y=xi_tau_inv(alpha-gamma); ta=sqrt(s_wig8); ca=c_wig8;
  y+=xi_tau_inv(beta-gamma); tb=sqrt(s_wig8); cb=c_wig8;
  cos_sin(gamma);
  // we are pledged to return 
  // sigma_wig(gamma)=cos(gamma)*y + sin(gamma)*(ta+tb)
  if (nargout==1) return c8*y+s8*(ta+tb);
  else {
    // if nargout>1 then we also write
    // df_dx = -sin(gamma)*(y+2*d2y)-cos(gamma)*dy, where dy=-0.5*(ta+tb),
    // d2y=0.25*(cos(alpha-gamma)/ta + cos(beta-gamma)/tb)
    ta_plus_tb=ta+tb;
    df_dx=-s8*(y+0.5*(ca/ta + cb/tb))+0.5*c8*ta_plus_tb;
    return c8*y+s8*ta_plus_tb; 
  } // else  
}// sigma_wig--------------------------------------------------------------------------------------------------------


// ************************************************************************************************
//  The above code has been optimized (except for xi).  The code below awaits optimizing.  ********
// ************************************************************************************************

double entry_angle(double t){
  //alf=entry_angle(t) returns the entry angle of E_[t,t2], where
  //t belongs to the interval [0,pi/2) and t2 in (t,pi] is such
  //that the turning angle of E_[t,t2] equals del (a global variable).
  //If the global variable nargout==2, then we set df_dx = (d/dt) alf
  double t2,x,x2,s,s2,u,c,c2,alf,tau1,dx,dx2,Del_s,Del_x,temp,du,du2,dr,dr2;
  t2=tau_inv(tau(t)+del); x=xi(t); x2=xi(t2); s=sin(t); s2=sin(t2); c=cos(t);
  u=direction_angle(s2-s,x2-x); dx=s*s/sqrt(1+s*s);
  tau1=direction_angle(c,dx); alf=tau1-u;
  if (nargout==2){
    c2=cos(t2); dx2=s2*s2/sqrt(1+s2*s2); Del_s=s2-s; Del_x=x2-x; temp=1/(Del_s*Del_s + Del_x*Del_x);
    du=temp*(-dx*Del_s+Del_x*c); du2=temp*(dx2*Del_s - Del_x*c2);
    dr=2*s/sqrt(1+s*s); dr2=2*s2/sqrt(1+s2*s2);
    df_dx=dr-du-du2*dr/dr2;
  } // if
  return alf;
}// entry_angle --------------------------------------------------------------------------------------

double E_theta(double t){
// th=E_theta(t) returns theta(t) and writes the first derivative 
// to the global variable df_dx if nargout>1.
double s,s2,xi0,w;
s=sin(t); xi0=xi(t);
if (nargout>1){
  s2=s*s; w=sqrt(1+s2);
  df_dx = 2*s/w-(s*s2/w-xi0*cos(t))/(s2+xi0*xi0);
}// if
return tau(t) - atan(xi0/(1e-30+s));
}// E_theta-----------------------------------------------------------------------------------------------------

double sqrt_E_theta(double t){
// sqrt_th=sqrt_E_theta(t) returns sqrt(theta(t)) and writes
// the first derivative to the global variable df_dx if nargout>1.
double s,s2,xi0,w,sqrt_th;
s=sin(t); xi0=xi(t);
sqrt_th=sqrt(tau(t) - atan(xi0/(1e-30+s)));
if (nargout>1){
  s2=s*s; w=sqrt(1+s2);
  df_dx = (2*s/w-(s*s2/w-xi0*cos(t))/(s2+xi0*xi0))/(2*sqrt_th);
}// if
return sqrt_th;
}// sqrt_E_theta-----------------------------------------------------------------------------------------------------

double E_theta_inv(double th){
// t=E_theta_inv(th) returns t in [0,t_max] such that E_theta(t) = th,
// provided th in [0,theta_max], where theta_max =  1.73779813442124 
// and t_max =  2.73347515529434 are found in arj069/Octave/exper7.m
// In case th >= theta_max, we return t_max and if th <= 0, we return 0;
// The global variable tol is read.
error_code8=0;
if (th<=0)
  return 0;
if (th>=1.73779813442124)
  return 2.73347515529434;
if (th>=0.25){
  a8=0.65029; b8=2.73347515529434; fa8=-1; fb8=1; N18=3; N28=16; N38=2;
  return bisection_Newton(E_theta,th);
}// if
else{
  a8=0; b8=0.6503; fa8=-1; fb8=1; N18=2; N28=16; N38=2;
  return bisection_Newton(sqrt_E_theta,sqrt(th));
}// else
}// E_theta_inv--------------------------------------------------------------------------------------------------------

int S_canonical(){
// error_code = S_canonical () 
// Inputs: (global variables) alpha, beta, FORM8[3], tol.
// Outputs: (global variables) BE8, FORM8[0], p18 and p28.
//
// Error Codes:
// 0. No errors.
// 1. alpha, beta do not satisfy canonical assumptions.
// 2. an error was encountered in Zone 1: entry_angle does not straddle target at a8 and b8.
// 3. an error was encountered in Zone 3: G_a < G_gamma even though alpha < alpha1_approx(beta).
//    This could mean bisection_method5 was unable to find an initial interval or the optimal choice
//    of gamma is G_a (neither of which is supposed to happen, according to current conjectures).
// Description in case FORM8[3]=0:
// Let u=[0,e^(I alpha)], v=[1,e^(I beta)] be a pair of s-feasible unit tangent vectors in canonical form: 
// alpha in [0,pi), |beta| <= alpha, beta >= alpha-pi. We find the distinguished curve in S(u,v) having 
// minimal bending energy and output Form, BE, p1, p2 to FORM8[0], BE8, p18, p28.
// Description in case FORM8[3]=1:
// Let u=[0,e^(I alpha)] be in canonical form for S(u,1)  with alpha in [0,pi).  We find the distinguished 
//curve in S(u,1) having minimal bending energy and output Form, BE, p1, p2 to FORM8[0], BE8, p18, p28.
// See Notebook Dec. 29, 2011, Notebook Feb. 13, 2012, and ../free_end_condition.tex for details.
double beta0,t1_b,y1_b,dy1_b,sigma_b,a,b,t2,t3,t4,target,temp1,temp2,gamma,y2,G_a,G_gamma;

// #################### the free case, see Notebook Feb. 13, 2012 ####################
if (FORM8[3]==1){

// test canonical assumptions
if (alpha<0 || alpha>M_PI_2+tol)
  {BE8=INFINITY; printf("S_canonical: error_code=1 in free case.\n"); return 1;}
// the trivial case
if (alpha==0)
  {BE8=0; FORM8[0]=0; p18=0; p28=0; return 0;}
else {
  FORM8[0]=1; temp1=E_theta_inv(alpha); p18=-temp1; p28=0;
  a=sin(temp1); b=xi(temp1); BE8=b*sqrt(a*a+b*b); return 0;
}// else

}// ####################  end of free case ##########################################

//##################### the not free case (ie s(u,v)) ###################################
// test canonical assumptions
if (alpha<0 || alpha>M_PI_2+tol || fabs(beta)>alpha+tol/2 || beta<alpha-M_PI-tol/2) {
  BE8=INFINITY; printf("S_canonical: error_code=1.\n");
  printf("alpha = %1.16lf \n",alpha); printf("beta = %1.16lf \n",beta); return 1;
}// if

// the trivial case alpha = 0
if (alpha<tol)
  {BE8=0; FORM8[0]=0; p18=0; p28=0; return 0;}

// The extreme case beta=alpha-pi
if (alpha-beta>M_PI-tol){
  BE8=M_D_SQ;
  FORM8[0]=2; p18=-M_D*cos(alpha)/sin(alpha); p28=0;
  return 0;
}// end of extreme case

// Compute beta0, sigma_b and, if beta<0, t1_b, y1_b, dy1_b and set a=alpha-pi, b=beta0
// NOTE:  The variable sigma here refers to sigma_wig in the documentation
if (beta<0){
  beta0=beta;
  t1_b=tau_inv(alpha-beta);  y1_b=xi(t1_b); dy1_b=sqrt(sin(alpha-beta));
  sigma_b=cos(beta)+sin(beta)*dy1_b/y1_b;
  }// if
else 
  {beta0=0; sigma_b=1; t1_b=0; y1_b=0;}
a=alpha-M_PI; b=beta0;
// end of that computation 

//                    Case 1: beta<0 and sigma_wig(beta)=0
if (fabs(sigma_b)<tol)
  {FORM8[0]=1; BE8=-y1_b*y1_b/sin(beta); p18=-t1_b; p28=0; return 0;}
//                    end of Case 1
//
//                    Case 2: beta<0 and sigma_wig(beta)<0
if (sigma_b<0){
  FORM8[0]=1;
  if (alpha+beta<tol) // ie beta=-alpha
    {t3=tau_inv(M_PI_2 - alpha); t4=M_PI-t3;}
  else{
    // b now takes on a different meaning, see page 2 of Notebook Dec. 29, 2011
    del = alpha-beta; b=tau_inv((M_PI-del)/2); target=beta; 
    a8=0; b8 = b;
    fa8=sign(entry_angle(a8)-target); fb8=sign(entry_angle(b8)-target);
    if (fa8*fb8!=-1)
      {printf("S_canonical: error_code=2.\n"); return 2;} // entry_angle does not straddle target at a8 and b8
    N18=8; N28=16; N38=2; tol2=tol;
    t3=bisection_Newton(entry_angle,target);
    //t3=solve_c_curve2(b,tol,target,del);
    t4=tau_inv(tau(t3)+del);
  }// else
  //z=E([t3,t4]); Dz=z(2)-z(1); BE8=abs(Dz)*imag(Dz); 
  temp1=sin(t4)-sin(t3); temp2=xi(t4)-xi(t3);
  BE8=sqrt(temp1*temp1 + temp2*temp2)*temp2;
  p18=-t4; p28=-t3;
  return 0;
}// if
//                    end of Case 2
//
//                    Case 3: beta>=0 or (beta<0 and sigma_wig(beta)>0)
nargout=1;
// Since this is the restricted elastic spline, we have seen that sigma(-pi/2)<0
// so we can start Newton bisection with a8=-pi/2 and b8=b.
  a8 = -M_PI_2; b8 = b; fa8=-1; fb8=1; N18=8; N28=16; N38=2; tol2=tol;
  gamma=bisection_Newton(sigma_wig,0);
  //gamma=solve_s_curve2(a,b,tol,alphabet);
  BE8=G_val(gamma); FORM8[0]=1;
  p18=-tau_inv(alpha-gamma); p28=tau_inv(beta-gamma);
  return 0;
//
}// S_canonical--------------------------------------------------------------------------------------------------------

double BE_quick (double *u, double *v,double psi,double recip_breadth){
// be=BE_quick(u,v,psi,recip_breadth) returns the energy of the configuration (u,v) if it is feasible and returns
// INFINITY otherwise.  For sake of speed, we assume that the calling function has pre-computed the direction angle
// psi (from u to v) and the reciprocal of the breadth recip_breadth.
// The global variables max_alphaQ8, lamQ8 and tol are read and the global variable flag_change8 is written; however
// flag_change8 is used here to indicate whether the configuration is feasible.
bool flag_free=false;
double u2,v2,alf,bet;
double be,ca,sa,ca2,sa2,ca3,sa3,ca4,cb,sb,cb2,sb2,cb3,sb3,cb4;
//
u2=u[DIR]; v2=v[DIR]-v[CORNER]; alf=0; bet=0;
if (v[STATE]==FREE_CLAMPED || v[STATE]==FREE) 
  flag_free=true;
else
  bet=mod2pi(v2*M_PI_180-psi);
if (u[STATE]>FREE_CLAMPED) // ie u[STATE] equals CLAMPED_FREE or FREE
  {flag_free=true; alf=bet; bet=0;}
else
  alf=mod2pi(u2*M_PI_180-psi);
//
// I am assuming that Energy(alpha,beta) to justify the above.  We now are in one of two cases:
// 1. the curve is free at v (ie flag_free==true) or
// 2. the curve is not free at both u and v.
// Note that if the curve is free at both u and v, then you'll be in case 1 with alf=bet=0.
u2=fabs(alf)-tol; v2=fabs(bet)-tolQ8; flag_change8=true;
if (u2 > alpha_max8 || u2>M_PI_2 || v2 > alpha_max8 || v2>M_PI_2 )
  {flag_change8=false;  return INFINITY;} // printf("BE_quick: be=inf, alf= %e, bet= %e \n",alf,bet);
// The above is probably method independent; the rest is method dependent
if (flag_free){
  ca=alf*alf;
  be=ca*(7.50338056832770e-01+ca*(-9.05618885092225e-02+ca*(3.25095423872991e-03-7.06130976566319e-04*ca)));
}// if
else{
Cos_Sin(alf,&ca,&sa);
ca2=ca*ca-sa*sa; sa2=2.0*ca*sa; // ca2 = cos(2*alpha)  ;  sa2 = sin(2*alpha)
ca3=ca*ca2 - sa*sa2; sa3 = ca*sa2 + ca2*sa; // ca3 = cos(3*alpha) ; sa3 = sin(3*alpha)
ca4=ca*ca3 - sa*sa3; // ca4 = cos(4*alpha)
Cos_Sin(bet,&cb,&sb);
cb2=cb*cb-sb*sb; sb2=2.0*cb*sb; // cb2 = cos(2*beta)   ;  sb2 = sin(2*beta)
cb3=cb*cb2 - sb*sb2; sb3=cb*sb2 + cb2*sb; // cb3 = cos(3*beta) ; sb3 = sin(3*beta)
cb4=cb*cb3 - sb*sb3; // cb4 = cos(4*beta)
be  =  1.44182008428724306 * 2.0;         // v_{0,0}
be +=  -0.86747266110211974 * (ca+cb);    // v_{1,0}
be +=  -0.48415933760982860 * 2.0*ca*cb;  // v_{1,1}
be +=  -0.25461111590405483 * (ca2+cb2);  // v_{2,0}
be +=   0.14881353072332951 * (ca2*cb+cb2*ca); // v_{2,1}
be +=  -0.00778958522977333 * 2.0*ca2*cb2;     // v_{2,2}
be +=   0.04271111564439376 * (ca3+cb3);       // v_{3,0}
be +=  -0.01617954671783946 * (ca3*cb+cb3*ca); // v_{3,1}
be +=  -0.00313248409135077 * (ca4+cb4);       // v_{4,0}
be +=   0.99446761418046825 * 2.0*sa*sb;       // w_{1,1}
be +=  -0.35272580461709380 * (sa2*sb+sb2*sa); // w_{2,1}
be +=   0.01460435032384606 * 2.0*sa2*sb2;     // w_{2,2}
be +=   0.05085553125277856 * (sa3*sb+sb3*sa); // w_{3,1}
}// else
return be*recip_breadth;
}// BE_quick----------------------------------------------------------------------------------------

int S_canonical_quick(){
// error_code = S_canonical_quick () 
// Inputs: (global variables) alpha, beta, FORM8[3], tol.
// Outputs: (global variables) BE8.
//
// Error Codes:
// 0. No errors.
// 1. alpha, beta do not satisfy canonical assumptions.
//
// Description in case FORM8[3]=0:
// Let u=[0,e^(I alpha)], v=[1,e^(I beta)] be a pair of s-feasible unit tangent vectors in canonical form: 
// alpha in [0,pi/2], |beta| <= alpha. We approximate the minimal bending energy in S(u,v) and
// output it to BE8.  The approximation is described in Notebook Jan. 2, 2014 on page 8 where the coefficients
// are computed using Mike/arj075/Octave/exper9.m with h=0.5, alpha_max=90 and N=4.
// 
// Description in case FORM8[3]=1:
// Let u=[0,e^(I alpha)] with alpha in [0,pi/2].  We approximate the minimal bending energy
// in S(u,1) and output it to BE8.  The approximation is described in Notebook Jan. 2, 2014 on page 8 
// where the coefficients are computed using Mike/arj075/Octave/exper9b.m with h=0.125, alpha_max=90 and N=4.
// 
double ca,sa,ca2,sa2,ca3,sa3,ca4,cb,sb,cb2,sb2,cb3,sb3,cb4;

// #################### the free case, see Notebook Feb. 13, 2012 ####################
if (FORM8[3]==1){
// test canonical assumptions
if (fabs(alpha)>M_PI_2+tol)
  {BE8=INFINITY; printf("S_canonical_quick: error_code=1 in free case.\n"); return 1;}
ca=alpha*alpha;
BE8=ca*(7.50338056832770e-01+ca*(-9.05618885092225e-02+ca*(3.25095423872991e-03-7.06130976566319e-04*ca)));
return 0;
}// ####################  end of free case ##########################################

//##################### the not free case (ie s(u,v)) ###################################
// test canonical assumptions
if (fabs(alpha)>M_PI_2+tol || fabs(beta)>M_PI_2+tol) {
  BE8=INFINITY; printf("S_canonical_quick: error_code=1.\n");
  printf("alpha = %1.16lf \n",alpha); printf("beta = %1.16lf \n",beta); return 1;
}// if
cos_sin(alpha); ca=c8; sa=s8; 
ca2=ca*ca-sa*sa; sa2=2.0*ca*sa; // ca2 = cos(2*alpha)  ;  sa2 = sin(2*alpha)
ca3=ca*ca2 - sa*sa2; sa3 = ca*sa2 + ca2*sa; // ca3 = cos(3*alpha) ; sa3 = sin(3*alpha)
ca4=ca*ca3 - sa*sa3; // ca4 = cos(4*alpha)
cos_sin(beta); cb=c8; sb=s8; 
cb2=cb*cb-sb*sb; sb2=2.0*cb*sb; // cb2 = cos(2*beta)   ;  sb2 = sin(2*beta)
cb3=cb*cb2 - sb*sb2; sb3=cb*sb2 + cb2*sb; // cb3 = cos(3*beta) ; sb3 = sin(3*beta)
cb4=cb*cb3 - sb*sb3; // cb4 = cos(4*beta)
BE8  =  1.44182008428724306 * 2.0;         // v_{0,0}
BE8 +=  -0.86747266110211974 * (ca+cb);    // v_{1,0}
BE8 +=  -0.48415933760982860 * 2.0*ca*cb;  // v_{1,1}
BE8 +=  -0.25461111590405483 * (ca2+cb2);  // v_{2,0}
BE8 +=   0.14881353072332951 * (ca2*cb+cb2*ca); // v_{2,1}
BE8 +=  -0.00778958522977333 * 2.0*ca2*cb2;     // v_{2,2}
BE8 +=   0.04271111564439376 * (ca3+cb3);       // v_{3,0}
BE8 +=  -0.01617954671783946 * (ca3*cb+cb3*ca); // v_{3,1}
BE8 +=  -0.00313248409135077 * (ca4+cb4);       // v_{4,0}
BE8 +=   0.99446761418046825 * 2.0*sa*sb;       // w_{1,1}
BE8 +=  -0.35272580461709380 * (sa2*sb+sb2*sa); // w_{2,1}
BE8 +=   0.01460435032384606 * 2.0*sa2*sb2;     // w_{2,2}
BE8 +=   0.05085553125277856 * (sa3*sb+sb3*sa); // w_{3,1}
//qDebug() << "hello quick world";
return 0; 
}// S_canonical_quick--------------------------------------------------------------------------------------------------

bool is_feasible_lite (double *u, double *v,double psi,double dummy_var=M_PI_PLUS){
// b=is_feasible_lite(u,v,psi) returns true if the configuration (u,v) is feasible and returns false otherwise.  Here
// u and v are nodes (see README). Additionally, the global variables associated with canonical form are written.
// Specifically, alpha, beta and FORM8[1,2,3], where FORM8 = [Form, flag_reversl, flag_reflection, flag_free] and flag_free
// equals 1 if one or both of the nodes is free and equals 0 otherwise.
// After running is_feasible_lite, S_canonical can be run.
bool u_is_free=false, v_is_free=false;
double u2,v2;
FORM8[1]=0; FORM8[2]=0;; FORM8[3]=0; // default options
if (u[STATE]>FREE_CLAMPED) // ie u[STATE] equals CLAMPED_FREE or FREE
  u_is_free=true;
if (v[STATE]==FREE_CLAMPED || v[STATE]==FREE)
  v_is_free=true;
u2=u[DIR]; v2=v[DIR]-v[CORNER];
if (!u_is_free && !v_is_free){
  alpha=mod2pi(u2*M_PI_180-psi); beta=mod2pi(v2*M_PI_180-psi);
  if (fabs(beta)>fabs(alpha)) // then we apply a direction reversal
    {u2=alpha; alpha=beta; beta=u2; FORM8[1]=1;}
  if (alpha<0) // then we apply a reflection
    {alpha=-alpha; beta=-beta;; FORM8[2]=1;}
  if (alpha>alpha_max8+tol || alpha>M_PI_2+tol)
    return false;
  else
    return true;
}// if
FORM8[3]=1;
if (u_is_free && v_is_free)
  {alpha=0; return true;}
// Now we are in the case where exactly one of u and v is free.
if (u_is_free)
  {alpha=mod2pi(v2*M_PI_180-psi); FORM8[1]=1;}
else
  alpha=mod2pi(u2*M_PI_180-psi);
if (alpha<0)
  {alpha=-alpha; FORM8[2]=1;}
// if (alpha<M_PI)
if (alpha <= alpha_max8+tol && alpha <= M_PI_2+tol)
  return true;
else
  return false;
}// is_feasible_lite----------------------------------------------------------------------------------------

bool is_feasible (double *u, double *v){
// b=is_feasible(u,v) returns true if the configuration (u,v) is feasible and returns false otherwise.  Here
// u and v are nodes (see README). Additionally, the global variables associated with canonical form are set.
// Specifically, alpha, beta, breadth8 and FORM8[1,2,3], where FORM8 = [Form, flag_reversl, flag_reflection, style] and style 
// equals 1 if one or both of the nodes is free and equals 0 otherwise.
// After running is_feasible, S_canonical can be run.
double psi;
psi=direction_angle(v[X]-u[X],v[Y]-u[Y]);
return is_feasible_lite(u,v,psi);
}// is_feasible--------------------------------------------------------------------------------------------------------

int S_general(double *u, double *v){
// error_code = S_general (*u, *v) does the same thing as S_general, where u and v are two nodes.
// Global variables: tol FORM8[0] (read after calling S_canonical)
//                   BE8, FORM8, p18 and p28 (written as outputs)
//                   alpha, beta             (written for S_canoncial to read)
// For convenience, we update u[BE] = BE8.
// The output FORM8 looks like
// FORM8=[Form, flag_reversal, flag_reflection,flag_free], where flag_reversal and flag_reflection
// are set to 1 if the reduction to canonical form employs a direction reversal or a reflection, respectively;
// and flag_free is set to 1 if at least one of u or v is free.
//
// S_general is essentially a front end to S_canonical.  The returned integer error_code actually comes
// from S_canonical and carries the following significance:
// Error Codes:
// 0. No errors.
// 1. alpha, beta do not satisfy canonical assumptions.
// 2. an error was encountered in Zone 1: entry_angle does not straddle target at a8 and b8.
// 3. an error was encountered in Zone 3: G_a < G_gamma even though alpha < alpha1_approx(beta).
//    This could mean bisection_method5 was unable to find a initial interval or the optimal choice
//    of gamma is G_a (neither of which is supposed to happen, according to current conjectures).
//
// Final Note: the case when (u1,u2) coincides with (v1,v2) is allowed and treated as having Form = 0 and
// bending energy 0. This may be useful for free end conditions.
int error_code;
double breadth;
if (is_feasible(u,v))
  {breadth=breadth8; error_code=S_canonical(); BE8=BE8/breadth; u[BE]=BE8; return error_code;}
else
  {BE8=INFINITY; u[BE]=BE8; return 1;}
}// S_general ----------------------------------------------------------------------------------------------

complex E(double t){
  complex z;
  z.real = sin(t);
  z.imag = xi(t);
  return z;
}// E------------------------------------------------------------------------------------------------------

int make_piece (double *u,double *v,piece *p){
// error_code = make_piece(u,v,p) makes the piece p representing the distinguished curve
// connecting node u to node v.  The returned value error code is the result of calling S_general(u,v).
// See Notebook Jan. 6, 2012, S_general.m and S_canonical.m for more details.
int error_code;
double u_save[6], v_save[6];
complex c1,c2,w1,w2,z1,z2;
piece q;
databall2node(u,u_save); databall2node(v,v_save);
error_code = S_general(u,v);

//printf("make_piece: u[BE]= %e \n", u[BE]);

if (error_code==0){
q.form = (short) FORM8[0];
q.flag_reversal = (bool) FORM8[1];
q.flag_reflection = (bool) FORM8[2];
q.p1 = p18;
q.p2 = p28;

if (q.form==0){ // Special case of a line segment. In this case, we arrange things so that
                 // the similarity transformation T(z) = c1 z + c2 maps the interval [0,1]
                 // onto the desired line segment.
  c1.real=v[X]-u[X]; c1.imag=v[Y]-u[Y];
  c2.real=u[X];      c2.imag=u[Y];
  q.c1=c1; q.c2=c2; *p=q; return error_code;
}// if

// the following two segments of code have been inserted to ensure that pieces of form=1
// do not have reflections or reversals.  We simply modify p1 and p2 to achieve this.
if (q.form==1 && q.flag_reflection)
  {q.p1=q.p1+M_PI; q.p2=q.p2+M_PI; q.flag_reflection=false;}
if (q.form==1 && q.flag_reversal)
  {temp8=q.p1; q.p1=-q.p2; q.p2=-temp8; q.flag_reversal=false;}

if (q.flag_reversal)
  {w1.real=v[X]; w1.imag=v[Y]; w2.real=u[X]; w2.imag=u[Y];}
else
  {w1.real=u[X]; w1.imag=u[Y]; w2.real=v[X]; w2.imag=v[Y];}

if (q.form==1)
  {z1=E(q.p1); z2=E(q.p2);}
else // ie q.form==2
  {z1.real=0; z1.imag=-M_D; z2=E(q.p2); z2.real += q.p1;}

if (q.flag_reflection)
  {z1=complex_conj(z1); z2=complex_conj(z2);}

c1=complex_div( complex_sub(w2,w1)  ,  complex_sub(z2,z1)  );
c2=complex_sub( w1   ,  complex_mult(c1,z1)  );
q.c1=c1; q.c2=c2;
*p=q;
}// if error_code==0
node2databall(u,u_save); node2databall(v,v_save);
return error_code;
}// make_piece---------------------------------------------------------------------------------

int piece_info (piece p, double *p_info, bool flag_full=false){
/* error_code=piece_info(p,p_info)
Given a piece p, we "return" a double array p_info holding the following fields:
p_info[0] = initial direction angle (in radians)
p_info[1] = terminal direction angle (in radians)
p_info[2] = initial signed curvature k
p_info[3] = terminal signed curvature k
p_info[4] = initial dk/ds (derivative of signed curvature w.r.t. arclength)
p_info[5] = terminal dk/ds 
   The following three fields relate to the segment lam*R_[t1,t2] of rectangular
   elastica which is directly congruent to our piece when p.form=1.  
   If p.form==0, then the piece is a line segment and we set p_info[6]=0,
   p_info[7]=|p.c1|, and p_info[8]=0.  Note that p.form==0 if and only if p_info[8]==0.
p_info[6] = p.p1 (the parameter t1)
p_info[7] = p.p2 (the parameter t2)
p_info[8] = lam (the scale factor)
The returned integer error_code is 0, if p.form = 0 or 1, and is 1 otherwise.
If the boolean flag_full==false (it's default), then only fields 0,1,2,3 are computed
(see Notebook May 22, 2016 for details)
*/
double lam;
complex c1,c;
c1=p.c1;
lam=complex_abs(c1);
if (p.form==0){
  p_info[0]=direction_angle(c1.real,c1.imag); p_info[1]=p_info[0]; 
  p_info[2]=0; p_info[3]=0; p_info[4]=0; p_info[5]=0; 
  //p_info[6]=0; p_info[7]=lam; p_info[8]=0; 
  return 0;}
if (p.form!=1)
  return 1;
// so it must be the case that p.form=1
// initial values
cos_sin(p.p1); 
p_info[2]=2.0*s8/lam;
p_info[4]=2.0*c8*sqrt(1.0+s8*s8)/(lam*lam);
c.real=c8*sqrt(1.0+s8*s8); c.imag=s8*s8; c=complex_mult(c,c1);
p_info[0]=direction_angle(c.real,c.imag);
// terminal
cos_sin(p.p2); 
p_info[3]=2.0*s8/lam;
p_info[5]=2.0*c8*sqrt(1.0+s8*s8)/(lam*lam);
c.real=c8*sqrt(1.0+s8*s8); c.imag=s8*s8; c=complex_mult(c,c1);
p_info[1]=direction_angle(c.real,c.imag);
return 0;
}// piece_info---------------------------------------------------------------------------------

double normal_forces (double **H, int n, complex *N, double *F){
/* Given a closed elastic spline represented by H (preferably optimized), we return a complex
array N and a double array F 
whereby N[i] = unit normal direction at node i,
        F[i] = the normal component of the supporting force at node i.
The returned double corresponds to ||F||
*/
int i,ip1,error_code,m;
double p_info[6],norm_F;
bool flag_closed=is_closed(H,n);
complex C;
piece p;
for (i=0;i<n;i++)
  F[i]=0;
m=n;
if (!flag_closed) m=n-1;
for (i=0;i<m;i++){
  ip1=(i+1) % n;
  error_code=make_piece(H[i],H[ip1],&p);
  error_code+=piece_info(p,p_info);
  cos_sin(p_info[0]+M_PI_2); C.real=c8; C.imag=s8; N[i]=C;
  F[i] += 0.5*p_info[4]; F[ip1] -= 0.5*p_info[5];
  //F[i]=p_info[4];
}// for i=0
if (!flag_closed)
  {cos_sin(p_info[1]+M_PI_2); C.real=c8; C.imag=s8; N[n-1]=C;}
norm_F=0;
for (i=0;i<n;i++)
  norm_F += F[i]*F[i];
return sqrt(norm_F);
}// normal_forces

int S_curve_lite (piece p,double *x, double *y,int D=257){
// k = S_curve_lite(p,x,y,D)
// Given a piece p, we return k points (x,y) along the distinguished curve, represented by p.
// The density of points is determined by the input D, where the guiding principal is that a 
// u-turn will receive D points.  Regardless of all else, we maintain 2 <= k <= 2D.
// See Notebook Jan. 6, 2012, S_general.m and S_canonical.m for more details.
int k,j,k1;
double h,t,c;
complex z;
//
// We first construct the curve K and store it in (x,y)
if (p.form==0)
  {x[0]=0; y[0]=0; x[1]=1; y[1]=0; k=2;}
else if (p.form==1){
  k=int(double(D)*(p.p2-p.p1)*M_1_PI);
  if (k<2)
     k=2;
  if (k>2*D)
     k=2*D;
  h=(p.p2-p.p1)/double(k-1); t=p.p1;
  for (j=0; j<k; j++)
    {z=E(t); x[j]=z.real; y[j]=z.imag; t+=h;}
}// else if
else{ // ie p.form=2
  c=p.p1; k1=D; h=M_PI/double(k1-1); t=-M_PI;
  //t=linspace(-pi,0,N); z=E(t);
  for (j=0;j<k1;j++)
    {z=E(t); x[j]=z.real; y[j]=z.imag; t+=h;}
  if (c>0 && p.p2==0)
    {x[k1]=c; y[k1]=0; k1++;}
  k=k1;
  if (p.p2>0){
    k=int(double(D)*p.p2*M_1_PI);
    if (k<2)
       k=2;
    h=p.p2/double(k-1); t=0; k=k1+k;
    //n=2+round(N*p2/pi); t=linspace(0,p2,n);  z=[z,c+E(t)];
    for (j=k1; j<k; j++)
      {z=E(t); x[j]=z.real; y[j]=z.imag; x[j] += c; t+=h;}
  }// if
}// else
//
// Now we apply the transformation T(z) = c1 z + c2 or
// T(z) = c1 conj(z) + c2 to our points.
if (p.flag_reflection)
  for (j=0; j<k; j++){
    z.real=x[j]; z.imag=-y[j];
    z=complex_add(complex_mult(p.c1,z),p.c2); x[j]=z.real; y[j]=z.imag;
  }// for
else
  for (j=0; j<k; j++){
    z.real=x[j]; z.imag=y[j];
    z=complex_add(complex_mult(p.c1,z),p.c2); x[j]=z.real; y[j]=z.imag;
  }// for
return k;
}// S_curve_lite---------------------------------------------------------------------------------

int S_curve (double *u,double *v,double *x, double *y,int D=257){
// k = S_curve(u,v,x,y,D)
// Given two nodes u and v, we return k points (x,y)
// along the distinguished curve in S(u,v), where the number of 
// returned points is the returned integer. The density of points is controlled 
// by the input D; the guiding principal is that a u-turn will receive D points.
// The procedure begins with a call to make_piece, whose output error_code is copied
// to the global variable error_code8.  If error_code8==0, then the procedure is completed
// with a call to S_curve_lite and the number of points k is returned; otherwise, k=0 is returned.
// Hence, if k=0, then you'll know that error_code8 is positive.
piece p;
int k;
error_code8=make_piece(u,v,&p);
if (error_code8==0)
  {k = S_curve_lite(p,x,y,D); return k;}
else
  return 0;
}//S_curve----------------------------------------------------------------------------------------

complex point_in_the_middle (piece p,double t=0.5){
// point_in_the_middle(p,t)
// Given a piece p, we return a complex point z, roughly at position t along the curve, where
// 0<=t<=1.
complex z;
//
// We first construct the complex point z in canonical form
if (p.flag_reversal)
  t=1-t;
if (p.form==0)
  {z.real=t; z.imag=0;}
else if (p.form==1)
  z=E(p.p1*(1-t)+p.p2*t);
else // ie p.form=2
  z=E(-M_PI*(1-t));
//
if (p.flag_reflection)
    z.imag=-z.imag;
// Now we apply the transformation T(z) = c1 z + c2
z=complex_add(complex_mult(p.c1,z),p.c2);
return z;
}// point_in_the_middle---------------------------------------------------------------------------------

int optimal_direction1(double *u,double *v,double *w,double *Theta,bool flag_offset){
// error_code = optimal_direction1(u,v,w,Theta,flag_offset) (something like that) :
// Let u,v,w be three consecutive nodes, with state of v regular.
// Theta should be an array of angles to be considered as candidates for v[DIR].  
// This procedure finds the optimal choice of v[DIR] and updates u,v,w accordingly.
// Global variables N8,BE8,tol are read, while global variables alpha,beta,flag_change8
// are written (the latter = 0 if H is unchanged and = 1 if H has changed).
// The calling routine is responsible to set tol (eg 1e-13) and N8 = length of Theta.
// If error_code > 0, then it means there was a problem.
// See page 1 of Notebook Feb. 2, 2012
// If flag_offset==true, then the angles in Theta are used to offset the current angle in v[DIR].
// Note: flag_change8 is set to true or false depending on whether v[DIR] has been changed.
double be_opt,psi_u,psi_v,v2_opt,be_u,be_v,breadth_u,breadth_v,alpha_u,beta_u,v2_save;
int k,error_code,FORM8_3;
be_opt = u[BE]+v[BE]; v2_save=v[DIR]; v2_opt=v[DIR]; error_code=0;
psi_u=direction_angle(v[X]-u[X],v[Y]-u[Y]); breadth_u=breadth8;
psi_v=direction_angle(w[X]-v[X],w[Y]-v[Y]); breadth_v=breadth8;
flag_change8=false;
//printf("opt_dir1: entering be = %e \n",be_opt);
for (k=0;k<N8;k++){
  if (flag_offset)
    v[DIR]=mod360(v2_save+Theta[k]);
  else
    v[DIR]=Theta[k];
  if (is_feasible_lite(u,v,psi_u)) {
    alpha_u=alpha; beta_u=beta; FORM8_3=FORM8[3];
    if (is_feasible_lite(v,w,psi_v)) { 
      error_code += S_canonical(); 
      be_v=BE8/breadth_v;
      if (be_v<be_opt){
        alpha=alpha_u; beta=beta_u; FORM8[3]=FORM8_3;
        error_code += S_canonical(); 
        be_u=BE8/breadth_u;
        if (be_u+be_v+tol/2<be_opt){
          u[BE]=be_u; v[BE]=be_v; be_opt=be_u+be_v; v2_opt=v[DIR]; flag_change8=true;
        }// if
      }// if
    }// if
  }// if
}// for (k=0
v[DIR]=v2_opt;
/*if (v2_opt==v2_save)
  flag_change8=false;
else{
  flag_change8=true;
  printf("opt_dir1 u: %e,  %e,  %e,  %e \n", u[X],u[Y],u[DIR],u[BE]);
  printf("opt_dir1 v: %e,  %e,  %e,  %e \n", v[X],v[Y],v[DIR],v[BE]);
  printf("opt_dir1 w: %e,  %e,  %e,  %e \n", w[X],w[Y],w[DIR],w[BE]);
}// else */
return error_code;
}// optimal_direction1--------------------------------------------------------------------------------------

int optimal_direction1_quick(double *u,double *v,double *w,double *Theta,bool flag_offset){
// error_code = optimal_direction1_quick(u,v,w,Theta,flag_offset) (something like that) :
// Let u,v,w be three consecutive nodes, with state of v regular.
// Theta should be an array of angles to be considered as candidates for v[DIR].  
// This procedure finds the optimal choice of v[DIR] and updates u,v,w accordingly.
// Global variable N8, tolQ8 are read, while the global variable flag_change8
// is written (the latter = 0 if H is unchanged and = 1 if H has changed).
// The calling routine is responsible to set tolQ8 (eg 1e-13) and N8 = length of Theta.
// I don't think the output error_code gets set in quick mode.
// See page 1 of Notebook Feb. 2, 2012
// If flag_offset==true, then the angles in Theta are used to offset the current angle in v[DIR].
// Note: flag_change8 is set to true or false depending on whether v[DIR] has been changed.
double be_opt,psi_u,psi_v,v2_opt,be_u,be_v,breadth_u_inv,breadth_v_inv,alpha_u,beta_u,v2_save,x,y;
int k,error_code;
bool flag_change;
// double tolQ8; tolQ8=tol;
be_opt = u[BE]+v[BE]; v2_save=v[DIR]; v2_opt=v[DIR]; error_code=0;
//psi_u=direction_angle(v[X]-u[X],v[Y]-u[Y]); breadth_u_inv=1/breadth8;
x=v[X]-u[X]; y=v[Y]-u[Y]; breadth_u_inv=1/(sqrt(x*x+y*y)+1e-200);
psi_u=Cos_Sin_inv(breadth_u_inv*x,breadth_u_inv*y);
//psi_v=direction_angle(w[X]-v[X],w[Y]-v[Y]); breadth_v_inv=1/breadth8;
x=w[X]-v[X]; y=w[Y]-v[Y]; breadth_v_inv=1/(sqrt(x*x+y*y)+1e-200);
psi_v=Cos_Sin_inv(breadth_v_inv*x,breadth_v_inv*y);
flag_change=false;
for (k=0;k<N8;k++){
  if (flag_offset)
    v[DIR]=mod360(v2_save+Theta[k]);
  else
    v[DIR]=Theta[k];
  be_u=BE_quick(u,v,psi_u,breadth_u_inv);
  if (flag_change8 && be_u<be_opt){
    be_v=BE_quick(v,w,psi_v,breadth_v_inv);
        if (flag_change8 && be_u+be_v+tolQ8/2<be_opt){
          u[BE]=be_u; v[BE]=be_v; be_opt=be_u+be_v; v2_opt=v[DIR]; flag_change=true;
        }// if
  }// if
}// for (k=0
v[DIR]=v2_opt;
flag_change8=flag_change;
/*if (v2_opt==v2_save)
  flag_change8=false;
else
  flag_change8=true; */
return error_code;
}// optimal_direction1_quick--------------------------------------------------------------------------------------

bool S_spline_feasible(double **H,int n,double *Theta,int N,long Max_trials=1000000){
// f=S_spline_feasible(H,n,Theta,N,Max_restarts) : The input n equals the number of
// nodes contained in H, and the input N equals then number of degree angles contained
// in Theta. Outgoing direction angles at regular nodes of H are initialized as mean
// stencil directions, but then a feasible offset from Theta is chosen randomly.  If
// Max_trials attempts are made, but no feasible offset is found, then we stick to the
// mean stencil angle.
// Afterwards the bending energies are updated and f=true is returned if all pieces are
// feasible; otherwise f=false is returned  (error_code8 is also written).
// 
int i,i0,im1,ip1,error_code=0;
long j;
bool feasible;
double v_save,phi,psi0,psi[n];
double u_save[6],up1_save[6],um1_save[6];

for (i = 0; i < n; i++){
  ip1=(i+1) % n;
  databall2node(H[i],u_save); databall2node(H[ip1],up1_save);
  psi[i] = direction_angle(H[ip1][X]-H[i][X], H[ip1][Y]-H[i][Y]);
  node2databall(H[i],u_save); node2databall(H[ip1],up1_save);
}// for
// first we initialize to mean stencil directions
for (i=0; i<n; i++){
  if (H[i][STATE]==REGULAR){
    ip1=(i+1) % n; im1=(i+n-1) % n;
    databall2node(H[im1],um1_save); databall2node(H[ip1],up1_save);
    psi0=psi[i];
    phi=direction_angle(H[i][X]-H[im1][X],H[i][Y]-H[im1][Y])+M_PI_180*H[i][CORNER];
    phi=phi-psi0; MOD2PI(phi);
    psi0+=0.5*phi; // Now psi is the desired direction angle in radians
                  // but we want to convert it to an angle in degrees
    psi0=M_180_PI*psi0;
    H[i][DIR]=mod360(psi0);
    node2databall(H[im1],um1_save); node2databall(H[ip1],up1_save);
  }// if
}// for
// next we search for feasible random offsets
for (i = 0; i < n; i++){
  if (H[i][STATE] == REGULAR){
    ip1=(i+1) % n; im1=(i+n-1) % n;
    databall2node(H[im1],um1_save); databall2node(H[ip1],up1_save);
    feasible=is_feasible_lite(H[(im1)], H[i], psi[im1],alpha_max8);
    feasible=feasible & is_feasible_lite(H[(i)], H[ip1], psi[i],alpha_max8);
    if (feasible){
      feasible=false;
      j=Max_trials;
      v_save=H[i][DIR];
         while (!feasible && j>0){
           j--;
           H[i][DIR] = mod360(v_save + Theta[rand() % N]);
           feasible=is_feasible_lite(H[(im1)], H[i], psi[im1],alpha_max8);
           feasible=feasible & is_feasible_lite(H[(i)], H[ip1], psi[i],alpha_max8);
         }// while
         if(!feasible)
             H[i][DIR]=v_save;
    }// if
    node2databall(H[im1],um1_save); node2databall(H[ip1],up1_save);
  }// if
}// for
// Finally we test feasibility of each piece and update bending energies
feasible = true;
for (i = 0; i < n; i++){
  ip1=(i+1) % n;
  databall2node(H[i],u_save); databall2node(H[ip1],up1_save);
  if (is_feasible_lite(H[(i)], H[ip1], psi[i],alpha_max8))
    error_code += S_general(H[i], H[ip1]);  // this command updates the bending energies
  else
    feasible=false;
  node2databall(H[i],u_save); node2databall(H[ip1],up1_save);
}// for
error_code8=error_code;
return feasible;
}// S_spline_feasible-------------------------------------------------------------------

int S_spline_check_feasible(double **H,int n){
// k=S_spline_check_feasible(H,n) : The input n equals the number of
// nodes contained in H.  If all pieces of H are feasible, then k=0; 
// otherwise piece k was found to be not feasible 
// NOTE: piece k connects node k-1 to node k.
// If all pieces are feasible, then we update the bending energies, where
// the global variable error_code8 will hold the accumulated error_code outputs
// of calls to S_general (error_code8=0 indicates everything went well).
int i=0,ip1,error_code=0;
bool feasible=true;
double up1_save[6],u_save[6],psi;
while (feasible && i<n){
  ip1=(i+1) % n;
  databall2node(H[i],u_save); databall2node(H[ip1],up1_save);
  feasible=is_feasible(H[i],H[ip1]);
  node2databall(H[i],u_save); node2databall(H[ip1],up1_save);
  i++;
}// while
if (feasible){ // the purpose here is to update the bending energies in H
  for (i = 0; i < n; i++){
    ip1 = (i+1) % n; 
    databall2node(H[i],u_save); databall2node(H[ip1],up1_save);
    if (flag_quick8){
      psi=direction_angle(H[ip1][X]-H[i][X],H[ip1][Y]-H[i][Y]);
      H[i][BE]=BE_quick(H[i],H[ip1],psi,1/breadth8); 
      if (H[i][BE]>1e6) printf("S_spline_check_feasible: H[i][BE] = %e and breadth8 = %e \n",H[i][BE],breadth8);
    }// if
    else
      error_code += S_general(H[i], H[ip1]);  
    node2databall(H[i],u_save); node2databall(H[ip1],up1_save);
  }                                          
  error_code8=error_code; i=0;
}// if
return i;
}// S_spline_check_feasible-------------------------------------------------------------------

int S_spline_feasible_stencil(double **H,int n,int N,bool flag_check_feasibility=true){
// k=S_spline_feasible_stencil(H,n,N) : The input n equals the number of
// nodes contained in H and N is leftover from elastic splines, so ignore it.
// Outgoing direction angles at regular nodes of H are chosen to equal the
// mean stencil direction.
// Afterwards, k=S_spline_check_feasible(double **H,int n) is run (this updates the
// bending energies if all pieces are feasible).  However, if you want to skip that, then simply
// give the command S_spline_feasible_stencil(H,n,N,false)
// 
int i,im1,ip1;
double psi,phi;
double up1_save[6],u_save[6],um1_save[6];
/* for (i=0; i<n; i++)
    if (H[i][STATE]<0)
      {ip1=(i+1) % n; H[i][NODEDIR]=direction_angle(H[ip1][X]-H[i][X],H[ip1][Y]-H[i][Y]);} */
for (i=0; i<n; i++){
  if (H[i][STATE]<=REGULAR || H[i][STATE]>100){
    im1=(i+n-1) % n; ip1=(i+1) % n;
    databall2node(H[im1],um1_save); databall2node(H[i],u_save); databall2node(H[ip1],up1_save);
    psi=direction_angle(H[ip1][X]-H[i][X],H[ip1][Y]-H[i][Y]);
    phi=direction_angle(H[i][X]-H[im1][X],H[i][Y]-H[im1][Y])+M_PI_180*H[i][CORNER];
    phi=phi-psi; MOD2PI(phi);
    psi+=0.5*phi; // Now psi is the desired direction angle in radians
                  // but we want to convert it to an angle in degrees
    psi=M_180_PI*psi;
    H[i][DIR]=mod360(psi);
    node2databall(H[im1],um1_save); node2databall(H[i],u_save); node2databall(H[ip1],up1_save);
  }// if
}// for
if (flag_check_feasibility)
  return S_spline_check_feasible(H,n);
else
  return 0;
}// S_spline_feasible_stencil--------------------------------------------------------------------------------------

void S_spline_feasible_stencil3(double **H,int n,int N,int i){
// k=S_spline_feasible_stencil3(H,n,N,i) : The input n equals the number of
// nodes contained in H and N is left over from elastic splines, so ignore it.
// Visiting only node i, the node prior to i and the node after i, if the node is regular, then
// the outgoing direction angle is chosen from Theta to approximate the average stencil direction.
// The caller should keep in mind that
// 1. the bending energies in H are no longer valid and
// 2. the feasibility of the resultant H has not been tested.
// It is expected that the caller will immediately plot the resultant ps-curve, thereby updating
// the bending energies and determining the feasibility of H.
// 
int k,im1,ip1;
double psi,phi;
double up1_save[6],um1_save[6];
i=(i+n-1) % n;
for (k=0; k<3; k++){
  if (H[i][STATE]==REGULAR){
    ip1=(i+1) % n; im1=(i+n-1) % n;
    databall2node(H[im1],um1_save); databall2node(H[ip1],up1_save);
    psi=direction_angle(H[ip1][X]-H[i][X],H[ip1][Y]-H[i][Y]);
    phi=direction_angle(H[i][X]-H[im1][X],H[i][Y]-H[im1][Y])+M_PI_180*H[i][CORNER];
    phi=phi-psi; MOD2PI(phi);
    psi+=0.5*phi; // Now psi is the desired direction angle in radians
                  // but we want to convert it to an angle in degrees
    psi=M_180_PI*psi;
    H[i][DIR]=mod360(psi);
    node2databall(H[im1],um1_save); node2databall(H[ip1],up1_save);
  }// if
  i=(i+1) % n;
}// for
}// S_spline_feasible_stencil3-------------------------------------------------------------------------------------

int Optimize(double **H,int n,double *Theta,int N,bool flag_coarse,bool *flag_opt=NULL,double *GAMMA=NULL,int M=0){
// error_code=Optimize(H,n,Theta,N,flag_coarse,NULL,GAMMA,M) : The input n equals the number of
// nodes contained in H and the inputs N,M equal the number of degree angles contained in 
// Theta,GAMMA (resp.) It is assumed that H is feasible and the bending energies in H are accurate.
// The optimization style depends on whether flag_coarse equals true or false.
// If flag_coarse==true, then outgoing direction angles at regular nodes of H
// are chosen optimally from Theta and H is updated accordingly. This mode of operation
// is suitable when Theta is an equi-angular partition of (-180,180].
// If flag_coarse==false, then Theta is treated as offset candidates for outgoing
// direction angles as regular nodes of H.  This mode of operation is suitable when the
// current outgoing direction angles are optimal on a coarse partition, and one seeks
// nearby outgoing direction angles which are optimal on a finer partition. In this
// case, one would probably have N rather small (eg 8 or 16) and Theta clustered near 0.
//
// The form error_code=Optimize(H,n,Theta,N,flag_coarse,flag_opt,GAMMA,M) can be used if flag_opt
// is a bool array of length n, where flag_opt[i] = false indicates that node i needs a
// visit on the first pass.  For example, if, after optimizing H, the user moves the position
// of node i, then only nodes i-1,i and i+1 are in need of a visit on the first pass, so
// one would send flag_opt with flag_opt[j]=false for j=i-1,i,i+1 and true otherwise.
//
// NOTE: flag_offset was motivated by elastic splines where u-turn detection was of critical importance.
// However, for restricted elastic splines and cubic splines Mike wants flag_offset=true regardless.
// Rather than bother Hakim with this, we'll simply set flag_offset=true on line 4 and be done with it.
//
// The inputs GAMMA and M have been introduced in order to allow data balls.  What this entails is that
// the nodes under consideration are those obtained by offsetting NODEDIR with angles from GAMMA.  GAMMA
// is generated in much the same way as Theta, except that M is probably a fraction of N when 
// flag_coarse==true (eg., N=360, M=30).
int i,j,im1,ip1,error_code=0,iter=0, num_visits;
N8=N;
bool all_optimal = false, flag_optimal[n], flag_offset;
double um1_save[6], v[6], up1_save[6];
double v_BE_opt, be_opt, be_opt_save, v_DIR_opt, v_NODEDIR_opt, v_NODEDIR, bee, be;
flag_offset=true; // originally it was flag_offset=!flag_coarse;
v[CORNER]=0; v[STATE]=0; // used in case of data balls; these need only be set once
if (flag_opt==NULL)
  for (i = 0; i < n; i++)
    flag_optimal[i] = false;
else
  for (i = 0; i < n; i++)
    flag_optimal[i] = flag_opt[i];
bee=0; 
for (i = 0; i < n; i++){
   bee += H[i][BE]; //printf("i= %d, be[i]= %e \n",i,H[i][BE]);
}//for
printf("be= %e \n",bee); be=bee;
if (flag_quick8)
   printf("quick mode! /n");

while (!all_optimal && iter<20000000)    {
  all_optimal = true;
  iter+=1;
  num_visits=0;
  for (i = 0; i < n; i++){
    if (!flag_optimal[i] && (H[i][STATE] <= 0 || H[i][STATE] > 100)){
      num_visits+=1;
      im1 = (n + i-1) % n; databall2node(H[im1],um1_save);
      ip1 = (i+1) % n; databall2node(H[ip1],up1_save);
      flag_optimal[i] = true;
      if (H[i][STATE] == REGULAR){
        if (flag_quick8)
          error_code+=optimal_direction1_quick(H[im1], H[i], H[ip1], Theta, flag_offset);
        else
          error_code+=optimal_direction1(H[im1], H[i], H[ip1], Theta, flag_offset);
      }// if
      else if (H[i][STATE] > 100){
        databall2node(H[i],v);
        if (flag_quick8)
          error_code+=optimal_direction1_quick(H[im1], H[i], H[ip1], Theta, flag_offset);
        else
          error_code+=optimal_direction1(H[im1], H[i], H[ip1], Theta, flag_offset);
        node2databall(H[i],v);
      }// else if
      else {  //(H[i][STATE] < 0){ // in here, we use the language of u,v,w to refer to
                            // H[im1], H[i], H[ip1]
        v_BE_opt=H[i][BE]; be_opt=H[im1][BE]+v_BE_opt; be_opt_save=be_opt; 
        v_DIR_opt=H[i][DIR]; v_NODEDIR_opt=H[i][NODEDIR]; v[BE]=H[i][BE];
        for (j = 0; j < M; j++){
          v_NODEDIR=mod360(H[i][NODEDIR]+GAMMA[j]);
          cos_sin(v_NODEDIR*M_PI_180);
          v[X]=H[i][X]-H[i][STATE]*c8; v[Y]=H[i][Y]-H[i][STATE]*s8;
          v[DIR]=H[i][DIR];  
          if (flag_quick8)
            error_code+=optimal_direction1_quick(H[im1], v, H[ip1], Theta, flag_offset);
          else
            error_code+=optimal_direction1(H[im1], v, H[ip1], Theta, flag_offset);
            //printf("OPT: j= %d, H[im1][BE]= %e \n", j, H[im1][BE]);
          if (flag_change8){
            v_DIR_opt=v[DIR]; v_NODEDIR_opt=v_NODEDIR; 
          }// if
        }// for j
        be_opt=H[im1][BE]+v[BE]; flag_change8=true;
        //printf("Hakim: i=%d, H[0][BE]=%e \n", i,H[0][BE]);
        if (be_opt+tol/2<be_opt_save){
           H[i][DIR]=v_DIR_opt; H[i][NODEDIR]=v_NODEDIR_opt; H[i][BE]=v[BE];
           be=be-be_opt_save+be_opt;
        }// if
        else
          {flag_change8=false;}
      }// else
      if (flag_change8)    {
          all_optimal = false;
          flag_optimal[im1] = false;
          flag_optimal[i]=flag_coarse;
          flag_optimal[ip1] = false;
      }// if
      node2databall(H[im1],um1_save); node2databall(H[ip1],up1_save);
    }// if
  }// for
  bee=0; for (i = 0; i < n; i++) bee += H[i][BE]; //qDebug() << "be: " << bee;
  printf("iter= %d, num_visits= %d,  be_pred= %e, be= %e \n",iter,num_visits,be,bee);
  //all_optimal=true;
}// while

/*ofstream file;
    file.open("trash0.txt", ios::out);
    file.precision(20);
    for (i = 0; i < n; i++)    {
        for (j = 0; j < 6; j++)    {
            file << H[i][j] << " ";
        }
        file << "\n";
    }


for (i = 0; i < n; i++)
    printf("Optimize: i= %d, H[i][BE]= %e \n", i, H[i][BE]);  */
if (flag_opt!=NULL)
  for (i = 0; i < n; i++)
    flag_opt[i]=flag_optimal[i];
return error_code; 
}// Optimize--------------------------------------------------------------------------------------


} //namespace
