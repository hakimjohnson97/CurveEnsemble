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

namespace ElasticS
{

// Macros:
//
// If t belongs to [0,pi], then the following performs t=tau(t).  After exiting, it
// leaves c8=cos(t) and s8=sin(t)
#define TAU(t) cos_sin(t); 2*acos(c8*M_1_SQRT2)-M_PI_2;
struct piece { //FORM8 = [Form, flag_reversl, flag_reflection, flag_free] and flag_free
  short form;
  bool flag_reversal, flag_reflection;
  double p1,p2;
  complex c1,c2;
};

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
  //y=xi(t) returns the function y=xi(t), employed in Euler's rectangular elastica, defined by
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

double alpha1_approx(double beta){
//alpha1=alpha1_approx(beta) : In Notebook Jan. 21, 2012 and in exper27*.m we discussed
//a function alpha1(beta) which is relevant to S_canonical.m in the case where sigma(beta0)>0 
//and sigma(alpha-pi)>0.  The conjecture has it that the distinguished
//curve s(u,v) is of first form iff alpha <= alpha1(beta).  The function alpha1 cannot
//be found precisely, but we provide a polynomial approximation here which is believed
//to deviate from the true function by no more than 0.001. In the Notebook, alpha1 is
//defined on [-pi/2,beta1], but we extend alpha1 by defining alpha1(beta)=beta1 if beta>beta1.
//Note: beta1=1.72011753853771
double y;
beta=0.607733924934515*beta+-0.0453737830442189;
// The mapping beta \mapsto 0.607733924934515*beta+-0.0453737830442189
// maps the interval [-pi/2,beta1] onto [-1,1] where our polynomial was constructed.
if (beta>1)
  beta=1;
y = beta * -4.30672575833856e+02 +  -6.14896099468632e+01;
y = beta*y +   2.13017461156870e+03;
y = beta*y +   2.72682793016909e+02;
y = beta*y +  -4.48008646618534e+03;
y = beta*y +  -5.08277218030759e+02;
y = beta*y +   5.21861360897100e+03;
y = beta*y +   5.16844475306459e+02;
y = beta*y +  -3.67475698608237e+03;
y = beta*y +  -3.11377853339595e+02;
y = beta*y +   1.60197465370037e+03;
y = beta*y +   1.12952378094633e+02;
y = beta*y +  -4.25278111431553e+02;
y = beta*y +  -2.39660349640972e+01;
y = beta*y +   6.48514052677347e+01;
y = beta*y +   2.72898451240471e+00;
y = beta*y +  -5.09852513800545e+00;
y = beta*y +   1.06906023641883e-02;
y = beta*y +  -2.35803030807646e-01;
y = beta*y +  -3.39446455839552e-02;
y = beta*y +   2.15964512678850e+00;
  return y;
}// alpha1_approx ------------------------------------------------------------------------------

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
if (alpha<0 || alpha>=M_PI)
  {BE8=INFINITY; printf("S_canonical: error_code=1 in free case.\n"); return 1;}
// the trivial case
if (alpha==0)
  {BE8=0; FORM8[0]=0; p18=0; p28=0; return 0;}
// Branch on whether alpha is <= theta_bar or not
if (alpha <= 1.73779813442124){
  FORM8[0]=1; temp1=E_theta_inv(alpha); p18=-temp1; p28=0;
  a=sin(temp1); b=xi(temp1); BE8=b*sqrt(a*a+b*b); return 0;
}// if
else 
  {FORM8[0]=2; a=sin(alpha); p18=-cos(alpha)*M_D/a; p28=0; BE8=M_D_SQ/a; return 0;}

}// ####################  end of free case ##########################################

//##################### the not free case (ie s(u,v)) ###################################
// test canonical assumptions
if (alpha<0 || alpha>=M_PI || fabs(beta)>alpha+tol/2 || beta<alpha-M_PI-tol/2) {
  BE8=INFINITY; printf("S_canonical: error_code=1.\n");
  printf("alpha = %1.16lf \n",alpha); printf("beta = %1.16lf \n",beta); return 1;
}// if

// the trivial case alpha = 0
if (alpha<tol)
  {BE8=0; FORM8[0]=0; p18=0; p28=0; return 0;}

// The extreme case beta=alpha-pi
if (alpha-beta>M_PI-tol){
  BE8=M_D_SQ/sin(alpha);
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
//                            Case 3a: sigma_wig(alpha-pi)<=0
if (sigma_wig(a)<=0){ // ie sigma(beta0)>0 and sigma(alpha-pi)<=0
  a8 = a; b8 = b; fa8=-1; fb8=1; N18=8; N28=16; N38=2; tol2=tol;
  gamma=bisection_Newton(sigma_wig,0);
  //gamma=solve_s_curve2(a,b,tol,alphabet);
  BE8=G_val(gamma); FORM8[0]=1;
  p18=-tau_inv(alpha-gamma); p28=tau_inv(beta-gamma);
  return 0;
}// if                        end of Case 3a
//
//                            Case 3b: sigma_wig(alpha-pi)>0
t2=tau_inv(beta-alpha+M_PI); y2=xi(t2); G_a=(M_D+y2)*(M_D+y2)/sin(alpha);
G_gamma=1024*G_a; gamma=alpha-M_PI;// At present G_a is smaller than G_gamma
if (!flag_shortcut8 || alpha<alpha1_approx(beta)+0.001){
  a8=a; b8=b; fa8=-1; fb8=1;
  temp1=bisection_method5(sigma_wig,0);
  nargout=1;
  if (sigma_wig(temp1)<0){
    a8=temp1; fa8=-1; fb8=1; N18=8; N28=16; N38=2; tol2=tol;
    gamma=bisection_Newton(sigma_wig,0);
    G_gamma=G_val(gamma);
  }// if
}// if
if (flag_shortcut8 && G_a<G_gamma && alpha<alpha1_approx(beta)-0.001)
  { printf("S_canonical: error_code=3.\n"); return 3;} // This may be a counter-example to Experimental Shortcut
if (G_a<G_gamma){
  FORM8[0]=2; BE8=G_a;
  p18=-cos(alpha)*(M_D+y2)/sin(alpha)-sqrt(sin(beta-alpha+M_PI)); p28=t2;
}// if
else
  {FORM8[0]=1; BE8=G_gamma; p18=-tau_inv(alpha-gamma); p28=tau_inv(beta-gamma);}
return 0;//                   end of Case 3b
//                    end of Case 3
}// S_canonical--------------------------------------------------------------------------------------------------------


bool is_feasible_lite (double *u, double *v,double psi){
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
  // if (alpha>M_PI_MINUS || beta<alpha-M_PI_PLUS)
  if (alpha > alpha_max8 || alpha>M_PI_MINUS || beta<alpha-M_PI_PLUS)  // changed March 2, 2013
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
if (alpha <= alpha_max8 && alpha<M_PI) // changed March 2, 2013
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
complex c1,c2,w1,w2,z1,z2;
piece q;
error_code = S_general(u,v);

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
return error_code;
}// make_piece---------------------------------------------------------------------------------

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

// TO DO: Create a high-level function Bending_Energy above.
// *********************************************************************************************************
// ********* THE DIVIDING LINE BETWEEN S_curve_lib.h and the C++ class Curve *******************************
// *********************************************************************************************************
// The functions below should only call high-level functions above they should not use any of 
// the above global variables.
// TO DO: 
// 1. Use the high-level function Bending Energy, rather than S_canonical in optimal_direction*
// 2. Ween optimal_direction* from its use of the global variable N8.

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
for (k=0;k<N8;k++){
  if (flag_offset)
    v[DIR]=mod360(v2_save+Theta[k]);
  else
    v[DIR]=Theta[k];
  if (is_feasible_lite(u,v,psi_u)) {
    alpha_u=alpha; beta_u=beta; FORM8_3=FORM8[3];
    if (is_feasible_lite(v,w,psi_v)) {
      error_code += S_canonical(); be_v=BE8/breadth_v;
      if (be_v<be_opt){
        alpha=alpha_u; beta=beta_u; FORM8[3]=FORM8_3;
        error_code += S_canonical(); be_u=BE8/breadth_u;
        if (be_u+be_v+tol/2<be_opt){
          u[BE]=be_u; v[BE]=be_v; be_opt=be_u+be_v; v2_opt=v[DIR];
        }// if
      }// if
    }// if
  }// if
}// for (k=0
v[DIR]=v2_opt;
if (v2_opt==v2_save)
  flag_change8=false;
else
  flag_change8=true;
return error_code;
}// optimal_direction1--------------------------------------------------------------------------------------

int optimal_direction2(double *u,double *v,double *w,double *x,double *Theta,bool flag_offset){
// error_code = optimal_direction2(u,v,w,x,Theta,flag_offset) :
// We are given four consecutive nodes u,v,w,x and it is assume that v,w are regular nodes
// with s(v,w) being a u-turn. Theta is an array of candidates for v[DIR].
// Maintaining s(v,w) as a u-turn, we find the optimal choice of v[DIR] and then
// update u,v,w accordingly.
// Global variables N8,BE8,tol are read, while global variables alpha,beta,flag_change8
// are written (the latter = 0 if H is unchanged and = 1 if H has changed).
// The calling routine is responsible to set tol (eg 1e-13) and N8 = length of Theta.
// If error_code > 0, then it means there was a problem.
// See pages 2,3 of Notebook Feb. 2, 2012
// If flag_offset==true, then the angles in Theta are used to offset the current angle in v[DIR].
// Note: flag_change8 is set to true or false depending on whether v[DIR] has been changed.
double be_opt,psi_u,psi_v,psi_w,v2_opt,v2_save,w2_opt,be_u,be_v,be_w;
double abs_sin_alf,breadth_u,breadth_v,breadth_w,alpha_u,beta_u,corner;
int k,error_code,FORM8_3;
be_opt = u[BE]+v[BE]+w[BE]; v2_opt=v[DIR]; v2_save=v[DIR]; w2_opt=w[DIR]; error_code=0; corner=w[CORNER];
psi_u=direction_angle(v[X]-u[X],v[Y]-u[Y]); breadth_u=breadth8;
psi_v=direction_angle(w[X]-v[X],w[Y]-v[Y]); breadth_v=breadth8;
psi_w=direction_angle(x[X]-w[X],x[Y]-w[Y]); breadth_w=breadth8;
for (k=0;k<N8;k++){
  if (flag_offset)
    v[DIR]=mod360(v2_save+Theta[k]);
  else
    v[DIR]=Theta[k];
  abs_sin_alf=fabs(sin(v[DIR]*M_PI_180-psi_v));
  if (M_D_SQ < abs_sin_alf * breadth_v * be_opt){
    be_v=M_D_SQ/(abs_sin_alf*breadth_v);
    w[DIR]=mod360(v[DIR]+180+corner);
    if (is_feasible_lite(u,v,psi_u)) {
      alpha_u=alpha; beta_u=beta; FORM8_3=FORM8[3];
      if (is_feasible_lite(w,x,psi_w)) {
        error_code += S_canonical(); be_w=BE8/breadth_w;
        if (be_v+be_w<be_opt){
          alpha=alpha_u; beta=beta_u; FORM8[3]=FORM8_3;
          error_code += S_canonical(); be_u=BE8/breadth_u;
          if (be_u+be_v+be_w+tol/2<be_opt){
            u[BE]=be_u; v[BE]=be_v; w[BE]=be_w; be_opt=be_u+be_v+be_w; v2_opt=v[DIR]; w2_opt=w[DIR];
          }// if
        }// if
      }// if
    }// if
  }// if
}// for (k=0
v[DIR]=v2_opt; w[DIR]=w2_opt;
if (v2_opt==v2_save)
  flag_change8=false;
else
  flag_change8=true;
return error_code;
}// optimal_direction2--------------------------------------------------------------------------------------

int optimal_direction3(double *u,double *v,double *v1,double *w,double *x,double *Theta,bool flag_offset){
// error_code = optimal_direction2(u,v,v1,w,x,Theta,flag_offset) :
// We are given five consecutive nodes u,v,v1,w,x and it is assume that v,v1,w are regular nodes
// with s(v,v1) and s(v1,w) being u-turns. Theta is an array of candidates for v[DIR].
// Maintaining s(v,v1) and s(v1,w) as u-turns, we find the optimal choice of v[DIR] and then
// update u,v,v1,w accordingly.
// Global variables N8,BE8,tol are read, while global variables alpha,beta,flag_change8
// are written (the latter = 0 if H is unchanged and = 1 if H has changed).
// The calling routine is responsible to set tol (eg 1e-13) and N8 = length of Theta.
// If error_code > 0, then it means there was a problem.
// See pages 2,3 of Notebook Feb. 2, 2012
// If flag_offset==true, then the angles in Theta are used to offset the current angle in v[DIR].
// Note: flag_change8 is set to true or false depending on whether v[DIR] has been changed.
double be_opt,psi_u,psi_v,psi_v1,psi_w,v2_opt,v12_opt,v2_save,w2_opt,be_u,be_v,be_v1,be_w;
double abs_sin_alf,abs_sin_alf1,breadth_u,breadth_v,breadth_v1,breadth_w,alpha_u,beta_u,corner,corner1;
int k,error_code,FORM8_3;
be_opt = u[BE]+v[BE]+v1[BE]+w[BE]; v2_opt=v[DIR]; v2_save=v[DIR]; v12_opt=v1[DIR]; w2_opt=w[DIR]; 
error_code=0; corner=v1[CORNER]; corner1=w[CORNER];
psi_u=direction_angle(v[X]-u[X],v[Y]-u[Y]); breadth_u=breadth8;
psi_v=direction_angle(v1[X]-v[X],v1[Y]-v[Y]); breadth_v=breadth8;
psi_v1=direction_angle(w[X]-v1[X],w[Y]-v1[Y]); breadth_v1=breadth8;
psi_w=direction_angle(x[X]-w[X],x[Y]-w[Y]); breadth_w=breadth8;
for (k=0;k<N8;k++){
  if (flag_offset)
    v[DIR]=mod360(v2_save+Theta[k]);
  else
    v[DIR]=Theta[k];
  v1[DIR]=mod360(v[DIR]+180+corner);
  w[DIR]=mod360(v1[DIR]+180+corner1);
  abs_sin_alf=fabs(sin(v[DIR]*M_PI_180-psi_v));
  abs_sin_alf1=fabs(sin(v1[DIR]*M_PI_180-psi_v1));
  if (abs_sin_alf>tol && abs_sin_alf1>tol){
    be_v=M_D_SQ/(abs_sin_alf*breadth_v);
    be_v1=M_D_SQ/(abs_sin_alf1*breadth_v1);
    if (be_v+be_v1 < be_opt && is_feasible_lite(u,v,psi_u)) {
      alpha_u=alpha; beta_u=beta; FORM8_3=FORM8[3];
      if (is_feasible_lite(w,x,psi_w)) {
        error_code += S_canonical(); be_w=BE8/breadth_w;
        if (be_v+be_v1+be_w<be_opt){
          alpha=alpha_u; beta=beta_u; FORM8[3]=FORM8_3;
          error_code += S_canonical(); be_u=BE8/breadth_u;
          if (be_u+be_v+be_v1+be_w+tol/2<be_opt){
            u[BE]=be_u; v[BE]=be_v; v1[BE]=be_v1; w[BE]=be_w; 
            be_opt=be_u+be_v+be_v1+be_w; v2_opt=v[DIR]; v12_opt=v1[DIR]; w2_opt=w[DIR];
          }// if
        }// if
      }// if
    }// if
  }// if
}// for (k=0
v[DIR]=v2_opt; v1[DIR]=v12_opt; w[DIR]=w2_opt;
if (v2_opt==v2_save)
  flag_change8=false;
else
  flag_change8=true;
return error_code;
}// optimal_direction3--------------------------------------------------------------------------------------

bool S_spline_feasible(double **H,int n,double *Theta,int N,long Max_restarts=10000000){
// f=S_spline_feasible(H,n,Theta,N,Max_restarts) : The input n equals the number of
// nodes contained in H, and the input N equals then number of degree angles contained
// in Theta. Outgoing direction angles at regular nodes of H
// are chosen randomly from Theta until feasibility is achieved, in which case the bending
// energies in H are updated and we return H and f=true (error_code8 is also written).
// On the other hand, if after Max_restarts restarts are made, feasibility
// is not attained, then the returned value is f=false.  For tuning purposes, the actual
// number of restarts used is written to the global variable Num_restarts8.
// 
int i,i0,ip1,error_code=0;
long j;
bool feasible;
double psi[n];
for (i = 0; i < n; i++){
  ip1=(i+1) % n;
  psi[i] = direction_angle(H[ip1][X]-H[i][X], H[ip1][Y]-H[i][Y]);
}// for
feasible = false; i0=0; j=Max_restarts;
while (feasible == false && j>0){
  j--; feasible = true; i=i0;
  if (H[i][STATE] == REGULAR)
    H[i][DIR] = Theta[rand() % N];
  i++;
  while (feasible && i < n)    {
    if (H[i][STATE] == REGULAR){
      H[i][DIR] = Theta[rand() % N];
      feasible=is_feasible_lite(H[(i-1)], H[i], psi[i-1]);
    }// if
    else{
      feasible=is_feasible_lite(H[(i-1)], H[i], psi[i-1]);
      if (feasible)
        i0=i+1;
    }// else
    i++;
  }// while 
  if (feasible)
    feasible=is_feasible_lite(H[(n-1)], H[0], psi[n-1]);
}// while
Num_restarts8=Max_restarts-j;
if (feasible)
    for (i = 0; i < n; i++)    {
        ip1 = (i+1) % n; 
        error_code += S_general(H[i], H[ip1]);  // the purpose here is to update
    }                                          // the bending energies in H
error_code8=error_code;
return feasible;
}// S_spline_feasible--------------------------------------------------------------------------------------

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
while (feasible && i<n){
  ip1=(i+1) % n;
  feasible=is_feasible(H[i],H[ip1]);
  i++;
}// while
if (feasible){
  for (i = 0; i < n; i++){
    ip1 = (i+1) % n; 
    error_code += S_general(H[i], H[ip1]);  // the purpose here is to update
  }                                          // the bending energies in H
  error_code8=error_code; i=0;
}// if
return i;
}// S_spline_check_feasible-------------------------------------------------------------------

int S_spline_feasible_stencil(double **H,int n,int N,bool flag_check_feasibility=true){
// k=S_spline_feasible_stencil(H,n,N) : The input n equals the number of
// nodes contained in H and N equals the number of angles in the equi-angular partition Theta.
// Outgoing direction angles at regular nodes of H are chosen from Theta to approximate the
// average stencil direction.
// Afterwards, k=S_spline_check_feasible(double **H,int n) is run (this updates the
// bending energies if all pieces are feasible).  However, if you want to skip that, then simply
// give the command S_spline_feasible_stencil(H,n,N,false)
// 
int i,im1,ip1;
double psi,phi,N360;
N360=360/double(N);
for (i=0; i<n; i++){
  if (H[i][STATE]==REGULAR){
    ip1=(i+1) % n; im1=(i+n-1) % n;
    psi=direction_angle(H[ip1][X]-H[i][X],H[ip1][Y]-H[i][Y]);
    phi=direction_angle(H[i][X]-H[im1][X],H[i][Y]-H[im1][Y])+M_PI_180*H[i][CORNER];
    phi=phi-psi; MOD2PI(phi);
    psi+=0.5*phi; // Now psi is the desired direction angle in radians
                  // but we want to convert it to an angle in Theta
    psi=ceil(M_1_PI*psi*(double(N)/2)-0.5);
    if (N!=360) psi*=N360;
    H[i][DIR]=mod360(psi);
  }// if
}// for
if (flag_check_feasibility)
  return S_spline_check_feasible(H,n);
else
  return 0;
}// S_spline_feasible_stencil--------------------------------------------------------------------------------------

void S_spline_rectify(double **H,int n,int N){
// k=S_spline_feasible_rectify(H,n,N) : The input n equals the number of
// nodes contained in H and N equals the number of angles in the equi-angular partition Theta.
// Outgoing direction angles at regular nodes of H are 'rectified' chosen from Theta to approximate the
// their current value.
// In order for Optimize (in coarse mode) to correctly identify u-turns when they arise, it is important that
// outgoing direction angles at regular nodes belong to Theta.  Consequently, if you run Optimize (coarse mode) starting
// from direction angles "AS IS", then it is strongly advisable to run S_spline_rectify first.  Regarding non-coarse
// mode, it is expected that Optimize is run from directions "AS IS" and so it is recommended to first run S_spline_rectify
// before running Optimize (non-coarse mode).
// 
int i;
double psi,three60_N,N_360;
three60_N=360/double(N);
N_360=double(N)/360;
for (i=0; i<n; i++){
  if (H[i][STATE]==REGULAR){
    psi=H[i][DIR];// Now psi is the desired direction angle (in degrees)
                  // but we want to convert it to an angle in Theta
    if (N!=360)
      {psi=ceil(N_360*psi-0.5);  psi*=three60_N;}
    else psi=ceil(psi-0.5);
    H[i][DIR]=mod360(psi);
  }// if
}// for
}// S_spline_rectify--------------------------------------------------------------------------------------

void S_spline_feasible_stencil3(double **H,int n,int N,int i){
// k=S_spline_feasible_stencil3(H,n,N,i) : The input n equals the number of
// nodes contained in H and N equals the number of angles in the equi-angular partition Theta.
// Visiting only node i, the node prior to i and the node after i, if the node is regular, then
// the outgoing direction angle is chosen from Theta to approximate the average stencil direction.
// The caller should keep in mind that
// 1. the bending energies in H are no longer valid and
// 2. the feasibility of the resultant H has not been tested.
// It is expected that the caller will immediately plot the resultant ps-curve, thereby updating
// the bending energies and determining the feasibility of H.
// 
int k,im1,ip1;
double psi,phi,N360;
N360=360/double(N);
i=(i+n-1) % n;
for (k=0; k<3; k++){
  if (H[i][STATE]==REGULAR){
    ip1=(i+1) % n; im1=(i+n-1) % n;
    psi=direction_angle(H[ip1][X]-H[i][X],H[ip1][Y]-H[i][Y]);
    phi=direction_angle(H[i][X]-H[im1][X],H[i][Y]-H[im1][Y])+M_PI_180*H[i][CORNER];
    phi=phi-psi; MOD2PI(phi);
    psi+=0.5*phi; // Now psi is the desired direction angle in radians
                  // but we want to convert it to an angle in Theta
    psi=ceil(M_1_PI*psi*(double(N)/2)-0.5);
    if (N!=360) psi*=N360;
    H[i][DIR]=mod360(psi);
  }// if
  i=(i+1) % n;
}// for
}// S_spline_feasible_stencil3-------------------------------------------------------------------------------------

int Optimize(double **H,int n,double *Theta,int N,bool flag_coarse,bool *flag_opt=NULL){
// error_code=Optimize(H,n,Theta,N,flag_coarse,NULL) : The input n equals the number of
// nodes contained in H, and the input N equals then number of degree angles contained
// in Theta. It is assumed that H is feasible and the bending energies in H are accurate.
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
// The form error_code=Optimize(H,n,Theta,N,flag_coarse,flag_opt) can be used if flag_opt
// is a bool array of length n, where flag_opt[i] = false indicates that node i needs a
// visit on the first pass.  For example, if, after optimizing H, the user moves the position
// of node i, then only nodes i-1,i and i+1 are in need of a visit on the first pass, so
// one would send flag_opt with flag_opt[j]=false for j=i-1,i,i+1 and true otherwise.
int i,im1,ip1,ip2,ip3,error_code=0;
N8=N;
bool all_optimal = false, flag_offset, flag_optimal[n];
flag_offset=!flag_coarse;
if (flag_opt==NULL)
  for (i = 0; i < n; i++)
    flag_optimal[i] = false;
else
  for (i = 0; i < n; i++)
    flag_optimal[i] = flag_opt[i];
while (!all_optimal)    {
  all_optimal = true;
  for (i = 0; i < n; i++){
    if ((!flag_optimal[i]) && (H[i][STATE] == REGULAR)){
      im1 = (n + i-1) % n;
      ip1 = (i+1) % n;
      ip2 = (i+2) % n;
      flag_optimal[i] = true;
      error_code+=optimal_direction1(H[im1], H[i], H[ip1], Theta, flag_offset);
      if (flag_change8)    {
          all_optimal = false;
          flag_optimal[im1] = false;
          flag_optimal[i]=flag_coarse;
          flag_optimal[ip1] = false;
      }// if
      if (n>2 && H[ip1][STATE] == REGULAR && fabs(mod360(H[i][DIR]-(H[ip1][DIR]-H[ip1][CORNER])-180)) < 1e-12){
        error_code+=optimal_direction2(H[im1], H[i], H[ip1], H[ip2], Theta, flag_offset);
        if (flag_change8)    {
            flag_optimal[im1] = false;
            flag_optimal[i] = false;
            flag_optimal[ip1] = false;
            flag_optimal[ip2] = false;
        }// if
        if (n>3 && H[ip2][STATE] == REGULAR && fabs(mod360(H[ip1][DIR]-(H[ip2][DIR]-H[ip2][CORNER])-180)) < 1e-12){
          ip3 = (i+3) % n;
          error_code+=optimal_direction3(H[im1], H[i], H[ip1], H[ip2], H[ip3], Theta, flag_offset);
          if (flag_change8)    {
              flag_optimal[im1] = false;
              flag_optimal[i] = false;
              flag_optimal[ip1] = false;
              flag_optimal[ip2] = false;
              flag_optimal[ip3] = false;
          }// if
        }// if
      }// if
    }// if
  }// for
}// while
if (flag_opt!=NULL)
  for (i = 0; i < n; i++)
    flag_opt[i]=flag_optimal[i];
return error_code;
}// Optimize--------------------------------------------------------------------------------------

} //namespace
