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


//#include <octave/oct.h>
#include "common_curve_lib.h"
#include <math.h>

double signum(double x)
{
    if (x > 0)
        return 1;
    else if (x == 0)
        return 0;
    else
        return -1;
}


double mod360(double a){
// a1=mod360(a) returns a1 in (-180,180] which is equivalent to a, modulo 360.
if (a>180)
  while (a>180)
    a-=360;
else
  while (a<= -180)
    a+=360;
return a;
}// mod360--------------------------------------------------------------------------------------------------------

void Cos_Sin(double t,double *c,double *s){
//  Assuming t in [-pi,3pi], cos_sin(t) sets c=cos(t) and s=sin(t).
bool s_neg=false, c_neg=false, cos_first=false;
if (t>M_PI) t-=M_2PI;
if (t<0) {t=-t; s_neg=true;}
if (t>M_PI_2) {t=M_PI-t; c_neg=true;}
if (t>M_PI_4) {cos_first=true; t=M_PI_2-t;}
// Our task here is to do s=sin(t), given that |t| <= pi/4 (nevermind that t>=0)
// c will hold t^2 while s will be used as an accumulator
*c=t*t; *s=(*c)*0.158885823618741691138890929053e-9;
*s=((*s) + -0.250506001755004223197431374134e-7)*(*c);
*s=((*s) + 0.275573125443566572888348008102e-5)*(*c);
*s=((*s) + -0.198412698259328449933580950308e-3)*(*c);
*s=((*s) + 0.833333333331647995167967482916e-2)*(*c);
*s=((*s) + -0.166666666666666003633159990497e0)*(*c);
*s=((*s) + 1)*t;
//
if (cos_first) {*c=(*s); *s=sqrt(1-(*s)*(*s));}
else *c=sqrt(1-(*s)*(*s));
if (s_neg) *s=-(*s);
if (c_neg) *c=-(*c);
}// Cos_Sin -----------------------------------------------------------------------------------------------------

double Cos_Sin_inv (double c, double s){
// t=Cos_Sin_inv() returns t in (-pi,pi] such that cos(t)=c and sin(t)=s,
bool reflect_x_axis=false, reflect_y_axis=false;
double t,y,w,temp;
if (s<0) {s=-s; reflect_x_axis=true;}
if (c<0) {c=-c; reflect_y_axis=true;}
//
if (s<M_SIN_PI_4)
  if (s<M_SIN_PI_8)
    if (s<M_SIN_PI_16)
     {t=M_PI_32; y=s*M_COS_PI_32-c*M_SIN_PI_32;}
    else
     {t=M_3_PI_32; y=s*M_COS_3_PI_32-c*M_SIN_3_PI_32;}
  else // ie M_SIN_PI_8 < s < M_SIN_PI_4
    if (s<M_SIN_3_PI_16)
     {t=M_5_PI_32; y=s*M_COS_5_PI_32-c*M_SIN_5_PI_32;}
    else
     {t=M_7_PI_32; y=s*M_COS_7_PI_32-c*M_SIN_7_PI_32;}
else // ie s > M_SIN_PI_4
  if (s<M_SIN_3_PI_8)
    if (s<M_SIN_5_PI_16)
     {t=M_9_PI_32; y=s*M_COS_9_PI_32-c*M_SIN_9_PI_32;}
    else
     {t=M_11_PI_32; y=s*M_COS_11_PI_32-c*M_SIN_11_PI_32;}
  else // ie M_SIN_3_PI_8 < s < M_SIN_PI_2
    if (s<M_SIN_7_PI_16)
     {t=M_13_PI_32; y=s*M_COS_13_PI_32-c*M_SIN_13_PI_32;}
    else
     {t=M_15_PI_32; y=s*M_COS_15_PI_32-c*M_SIN_15_PI_32;}
//t+=asin(y);
// Our task here is to do w=asin(y), given that |t| <= pi/4 (nevermind that t>=0)
// temp will hold y^2 while w will be used as an accumulator
temp=y*y; w=temp*0.0309791038033985847494561520516;
w=(w + 0.0446371117335725533083570850688)*temp;
w=(w + 0.0750000241324320065403450771277)*temp;
w=(w + 0.166666666625326475178247304944)*temp;
w=(w + 1.00000000000001982318723293716)*y;
t+=w;
//
if (reflect_y_axis) {t=M_PI-t; c=-c;}
if (reflect_x_axis) {t=-t; s=-s;}
return t;
}// Cos_Sin_inv -----------------------------------------------------------------------

void equi_angular_partition(double *Theta,int N){
// equi_angular_partition(Theta,N) : Assuming that Theta points to an array of
// double's of length N, we define Theta[i], for i=0,1,...,N-1 such that Theta will constitute
// an equi-angular partition of (-180,180].  For proper functioning of Optimize, N should be even.
int i;
for (i=0; i<N; i++)
  Theta[i]=(double(i+1)*double(360))/double(N)-double(180);
}//equi_angular_partition ----------------------------------------------------------------------------

void equi_angular_cluster(double *Theta,double h,int k){
// equi_angular_cluster(Theta,h,k) : Assuming that Theta points to an array of doubles's of length 2k,
// we define Theta[i], for i=0,1,...,2k-1 such that  Theta contains the set {-kh,...,-2h,-h,h,2h,...,kh}
int i,j;
j=0;
for (i=1;i<=k;i++)
  {Theta[j]=-double(i)*h; j++; Theta[j]=double(i)*h; j++;}
}// equi_angular_cluster--------------------------------------------------------------------------------

double bisection_method(double (*f) (double)  ,double target){
// c = bisection_method(f,target) solves f(x)=target using the bisection method, where
// f is a univariate function.  The initial and terminal interval is [a8,b8], where a8 and b8
// are global variables.  The iteration stops when |b8-a8| <= tol, where tol is also a global variable.
// The initial signs, sign(f(a8)-target) and sign(f(b8)-target), are held in fa8 and fb8, respectively.
double c;
int fc;
nargout=1; c=0.5*(a8+b8);
while (fabs(b8-a8)>tol){
  fc=sign((*f)(c)-target); // fc=(*f)(c) does fc=f(c)
  if (fc==fb8)
    b8=c;
  else if (fc==fa8)
    a8=c;
  else
    {a8=c; b8=c;}
  // if
  c=0.5*(a8+b8);
} // while
return c;
}// bisection_method--------------------------------------------------------------------------------------------------------

double bisection_method5(double (*f) (double)  ,double target){
//c=bisection_method5(f,target) attempts to find a feasible initial interval
//for solving f(x)=target using the bisection method, where f is a univariate function
//(see pages 4,5 of Notebook Dec. 29, 2011 for details).
//For simplicity, let's assume that f is continuously differentiable on [a,b].
//We assume both f(a) and f(b) are greater than target, while our goal is to find
//[c,d] subset of (a,b] such that f(c)<target<f(d). It is also assumed that f'(a)<0 and f'(b)>0, and
//KEY ASSUMPTION: there exists a unique x1 in (a,b) such that f'(x1)=0.
//It follows from the above that our problem has a solution if and only if
//f(x1)<0, in which case c=x1 will do. The outputs have the form C=[c,f(c)] and
//D=[d,f(d)] with the understanding the our problem has a solution if and only if f(c)<target.
//The role of the input tol is that |c-x1|<tol in case of NO SOLUTION (ie f(c)>=target).
//The function f is assumed to be implemented as [y,dy]=fname(x) where y=f(x) and
//dy=f'(x), x in (a,b).  The algorithm is designed so that f and f' never get evaluated
//at the endpoints a or b, so in practice, f needs only to be continuously differentiable
//on the open interval (a,b), where the meaning of f(a)>target, f(b)>target, f'(a)<0
//and f'(b)>0 is explained in the above cited Notebook entry.
//
// The initial and terminal interval is [a8,b8], where a8 and b8
// are global variables.  The iteration stops in one of three ways, indicated by global variable status8:
// 0. status8==0 the interval [a8,b8] has length <= tol (global variable) without finding f(c)<target.
// 1. status8==1 indicates f(c)<target and the desired interval is [c,b8].
// 2. status8==2 or status8==0 indicates f(c) >= target for all c.
// The initial signs, sign(f(a8)-target) and sign(f(b8)-target), are held in fa8 and fb8, respectively.
double c,y;
int fc,status;
nargout=2; status=0; c=0.5*(a8+b8);
while (status==0 && fabs(b8-a8)>tol){
 y=(*f)(c)-target;
 if (y<0)
   {status=1; a8=c;}
 else{
  fc=sign(y);
  if (fc==fb8)
    {b8=c; c=0.5*(a8+b8);}
  else if (fc==fa8)
    {a8=c; c=0.5*(a8+b8);}
  else
    status=2;
 }// else
} // while
return c;
}// bisection_method5--------------------------------------------------------------------------------------------------------

double bisection_Newton(double (*f) (double)  ,double target){
// c = bisection_Newton(f,target) uses the bisection method and Newton's method
// to solve the equation f(c)=target. Global variables N18, N28, N38 control the number of iterations
// in each phase (recommended values: N18=8, N28=16, N38=4. N1 is the number of initial bisection steps.
// The global variable tol is used to check the distance between
// two inputs to f while tol2 (global var.) is used to check the distance from an output of f to target.
// The initial and terminal interval is [a8,b8], where a8 and b8 are global variables.
// The initial signs, sign(f(a8)-target) and sign(f(b8)-target), are held in fa8 and fb8, respectively.
// After running N18 bisection steps, which produces a refined solution interval [a8,b8],
// the algorithm switches to Newton's method.  If a Newton iterate falls outside the
// refined interval (a8,b8), then the algorithm restarts (using recursion).  Only N38 restarts of
// bisection_Newton are allowed; the next (and final) restart will be a call to bisection_method.
// Also, N28 Newton steps will result in a call to bisection_method.
//
// Stopping Criterion: During the bisection phase, we stop if |b8-a8| <= tol.
// During the Newton phase, we maintain the variable flag_stop, initialize to 0. Each
// time a stopping criterion is satisfied, flag_stop is incremented by 1, and the
// iteration proceeds while flag_stop<3.  We respect two stopping criterion:
// 1. the distance between the previous and current iterate is < tol_in.
// 2. the distance between the current output and target is < tol_out.
// The variable flag_stop is also used to indicate that a restart, or call to bisection_method3, is needed:
// flag_stop is set to 4 if Newton wanders outside (a,b),
// flag_stop is set to 5 if N28 Newton iterations have been performed.
// For tuning purposes, the number of Newton iterations performed in the initial
// attempt is returned as iter8 (a global variable).
double c, Fc, c_prev;
int fc, iter, flag_stop;
//                                   bisection Phase
nargout=1; c=0.5*(a8+b8); iter=0;
while (iter<N18 && fabs(b8-a8)>tol){
  iter+=1;
  fc=sign((*f)(c)-target); // ie fc = sign( f(c) - target )
  if (fc==fb8)
    b8=c;
  else if (fc==fa8)
    a8=c;
  else
    {a8=c; b8=c;}
  // if
  c=0.5*(a8+b8);
} // while
if (fabs(b8-a8)>tol)
  flag_stop=0;
else
  flag_stop=3;
//                                 Newton Phase
// Note: flag_stop = 4 indicates that the Newton iteration has fallen
// outside the interval (a,b), while flag_stop=5 indicates that we
// have reached N28 Newton iterations.
iter=0; nargout=2;
while (iter<N28 && flag_stop<3){
   iter+=1; Fc=(*f)(c)-target; // ie Fc=f(c)-target and df_dx = f'(c)
   c_prev=c; c=c-Fc/df_dx;
   if (fabs(c-c_prev)<tol || fabs(Fc)<tol2)
     flag_stop+=1;
   if (c<=a8 || c>=b8)
      flag_stop=4;
}// while
if (flag_stop<3 && iter==N28)
  flag_stop=5;
//                                 Restart if necessary
if (flag_stop>3){
  if (flag_stop==4 && N38>0)
    {N38-=1; c=bisection_Newton((*f),target);}
  else
    {c=bisection_method((*f),target);}
}// if
iter8=iter;
return c;
}// bisection_Newton--------------------------------------------------------------------------------------------------------

//struct complex {double real,imag;};

complex complex_create(double x, double y){
    complex c;
    c.real = x;
    c.imag = y;
    return c;
}

complex complex_add(complex x,complex y){
// returns x + y
  x.real += y.real;
  x.imag += y.imag;
  return x;
}// complex_add--------------------------------------------------------------------------------------------

complex complex_sub(complex x,complex y){
// returns x - y
  x.real -= y.real;
  x.imag -= y.imag;
  return x;
}// complex_sub--------------------------------------------------------------------------------------------

complex complex_mult(complex x,complex y){
// returns x * y
  complex z;
  z.real = x.real * y.real - x.imag * y.imag;
  z.imag = x.real * y.imag + x.imag * y.real;
  return z;
}// complex_mult--------------------------------------------------------------------------------------------

complex complex_mult(double t,complex x){
// returns t*x
  x.real *= t;
  x.imag *= t;
  return x;
}// complex_mult--------------------------------------------------------------------------------------------

complex complex_conj(complex x){
// returns conj(x)
  x.imag = - x.imag;
  return x;
}// complex_conj--------------------------------------------------------------------------------------------

complex complex_inv(complex x){
// returns 1/x
  double z2;
  z2 = x.real*x.real + x.imag*x.imag;
  x.real = x.real/z2;
  x.imag = - x.imag/z2;
  return x;
}// complex_inv--------------------------------------------------------------------------------------------

double complex_abs(complex x){
// returns |x|
  return sqrt(x.real*x.real + x.imag*x.imag);
}// complex_inv--------------------------------------------------------------------------------------------

complex complex_polar(complex x){
// returns a complex number z=r+i*theta such that
// r = abs(x) and theta is the direction angle, in degrees, of x.
  x.imag=direction_angle(x.real,x.imag)*M_180_PI;
  x.real=breadth8;
  return x;
}// complex_polar--------------------------------------------------------------------------------------------

complex complex_div(complex x,complex y){
// returns x / y
  double z2;
  complex z;
  z2 = y.real*y.real + y.imag*y.imag;
  y.real = y.real/z2;
  y.imag = - y.imag/z2;
  z.real = x.real * y.real - x.imag * y.imag;
  z.imag = x.real * y.imag + x.imag * y.real;
  return z;
}// complex_div--------------------------------------------------------------------------------------------

void mark_all(double **H,int n){
int i;
for (i = 0; i < n; i++)
  if (H[i][STATE]<0)
     H[i][STATE]=100-H[i][STATE];
}// mark_all-----------------------------------------------------------------------------

void unmark_all(double **H,int n){
int i;
for (i = 0; i < n; i++)
  if (H[i][STATE]>100)
     H[i][STATE]=100-H[i][STATE];
}// mark_all------------------------------------------------------------------------------

int extract_unmarked(double **H,int n){
/* extract_unmarked(H,n, S) runs through H and deletes all marked data balls.
   The returned integer is the resultant number of rows in H. */
int i,j,m=0;
for (i = 0; i < n; i++)
  if (H[i][STATE]<100){
    for (j=0; j<6; j++) H[m][j]=H[i][j];
    m++;
  }// if
return m;
}// extract_unmarked----------------------------------------------------------------------

double closest_point_thresh(double *w,double *z,double *x,double *y,int n,double thresh){
/* Let z be a point in R2 
and let x and y be n x 1 vectors interpreted as a piecewise-linear curve S.  If the minimum
distance from S to z is more than thresh, then the returned point w is a point
on S which is closest to z; otherwise w is a point on S whose distance to z is <= thresh.
The returned value is the distance from z to w.  Note: w and z should be double arrays
of length at least 2. 
This procedure is adapted from arj029/Octave/clossest_points.cc
*/
  int i;
  double z1, z2, x_r, x_i, y_r, y_i, u_r, u_i, t, min_dist2, w1, w2, Y1, Y2;
  double thresh2=thresh*thresh;
  z1=z[0]; z2=z[1];
  w1=x[0]; w2=y[0]; min_dist2=(w1-z1)*(w1-z1)+(w2-z2)*(w2-z2); i=0;
  while (i<n-1 && min_dist2>thresh2){
    x_r=x[i];   //real(S(i,0));
    x_i=y[i];   //imag(S(i,0));
    y_r=x[i+1]; //real(S(i+1,0));
    y_i=y[i+1]; //imag(S(i+1,0));
    u_r=y_r-x_r;
    u_i=y_i-x_i;
    if (u_r*z1+u_i*z2 > u_r*x_r+u_i*x_i){
      if (u_r*z1+u_i*z2 < u_r*y_r+u_i*y_i){
        t=(u_r*(z1-x_r)+u_i*(z2-x_i))/(u_r*u_r+u_i*u_i);
        Y1=(1-t)*x_r+t*y_r; Y2=(1-t)*x_i+t*y_i;   //Y=(1-t)*S(i,0)+t*S(i+1,0);
      }/* if */
      else{
        Y1=y_r; Y2=y_i;  //Y=S(i+1,0);
      }/* else */
    }/* if */
    else{
      Y1=x_r; Y2=x_i;  //Y=S(i,0);
    }/* else */
    if ((z1-Y1)*(z1-Y1)+(z2-Y2)*(z2-Y2) < min_dist2){   //if (abs(z-Y)<min_dist){
       w1=Y1; w2=Y2;  //w=Y;
       min_dist2=(w1-z1)*(w1-z1)+(w2-z2)*(w2-z2);  // abs(z-w);
    }/* if */
  i++;
  }/* while */
  w[0]=w1; w[1]=w2;
  return sqrt(min_dist2);
}// closest_point_thresh------------------------------------------------------------------

void anchor_data_balls(double **H,int n){
double u_save[6];
for (int i=0; i<n; i++)
  databall2node(H[i],u_save);
}// anchor_data_balls--------------------------------------------------------------------

double max_radius(double **H,int n){
double r,R=0;
for (int i=0; i<n; i++){
  if (H[i][STATE]<0)
    r=-H[i][STATE];
  else if (H[i][STATE]>100)
    r=H[i][STATE]-100;
  else
    r=0;
  if (r>R)
    R=r;
}// for i
return R;
}// max_radius---------------------------------------------------------------------------

