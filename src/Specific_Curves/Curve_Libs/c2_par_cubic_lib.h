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


namespace C2ParCubicS
{

struct piece {
  double h,a3,a2,a1,a0,b3,b2,b1,b0;
  // x(t)=a0+t*(a1+t*(a2+t*a3)), y(t)=b0+t*(b1+t*(b2+t*b3)), 0 <= t <= h.
};

int make_piece (double *u,double *v,piece *q){
// error_code = make_piece(u,v,p) makes the piece p representing the distinguished curve
// connecting node u to node v.  The returned value error code is 0 for no problems and 1
// if there's a problem.
int error_code=0;
double h,p0,p1,q0,q1,h_inv,temp;  // as on page 5 of Notebook Sept. 14, 2015
piece p;
h=u[BE]; p.h=h; h_inv=1/h; 
p0=u[X]; p1=u[DIR]; q0=v[X]; q1=v[DIR];
temp=h_inv*(q0-p0);
p.a3=h_inv*h_inv*(q1+p1-2*temp); p.a2=h_inv*(3*temp-2*p1-q1); p.a1=p1; p.a0=p0;
p0=u[Y]; p1=u[CORNER]; q0=v[Y]; q1=v[CORNER];
temp=h_inv*(q0-p0);
p.b3=h_inv*h_inv*(q1+p1-2*temp); p.b2=h_inv*(3*temp-2*p1-q1); p.b1=p1; p.b0=p0;
*q=p;
return error_code;
}// make_piece---------------------------------------------------------------------------------

int S_curve_lite (piece p,double *x, double *y,int D=257){
// k = S_curve_lite(p,x,y,D)
// Given a piece p, we return k points (x,y) along the curve represented by p.
// The density of points is determined by the input D, where the guiding principal is that a 
// u-turn will receive D points.  Regardless of all else, we maintain 2 <= k <= 2D.
// See Notebook March 23, 2012 for more details.
int k,j;
double t,dt;
  k=2*D;
  if (k<2)
     k=2;
  if (k>2*D)
     k=2*D;
  dt=p.h/double(k-1); t=0;
  for (j=0; j<k; j++){
    x[j]=p.a0+t*(p.a1+t*(p.a2+t*p.a3));
    y[j]=p.b0+t*(p.b1+t*(p.b2+t*p.b3));
    t+=dt;
  }// for 
//k=2; x[0]=p.a1; y[0]=p.b1; x[1]=p.a2; y[1]=p.b2;
return k;
}// S_curve_lite---------------------------------------------------------------------------------

complex point_in_the_middle (piece p,double t=0.5){
   complex z; 
   t=t*p.h;
   z.real=p.a0+t*(p.a1+t*(p.a2+t*p.a3));
   z.imag=p.b0+t*(p.b1+t*(p.b2+t*p.b3));
   return z;
}//point_in_the_middle------------------------------------------------------------------------------


void time_increments(double **H,int n){
// this routine sets H[i][BE] = h for i=0,1,...,n-2, where
// h = abs(delta x) which is appropriate for classical cubic splines
int i,ip1;
double dx,dy,h;
if (flag_parametric8) {
  for (i=0; i<n; i++){
     ip1=(i+1) % n;
     dx=H[ip1][X]-H[i][X]; dy=H[ip1][Y]-H[i][Y]; h=sqrt(dx*dx+dy*dy);
     H[i][BE]=pow(h,lam8);
  }// for
}// if
else{
  for (i=0; i<n-1; i++){
    H[i][BE]=H[i+1][X]-H[i][X];
    if (fabs(H[i][BE])<tol)
      H[i][BE]=1;
  }// for
}// else
}// time_increments---------------------------------------------------------------------------------

void multi_step_classical_cubic(double **H,int n,int k){
// global variables: p18,p28 are read and written.  The inputs s'(x0) and s''(x0) are read from 
// p18, p28 and the outputs s'(xn) and s''(xn) are written to p18, p28 (see Notebook Sept 14, 2015).
// delta x values are found in H(:,BE) and y values are found in H(:,k).
// Resultant first derivative values are written to H(:,k+2).
int i,ip1,N;
double h,h_inv,q0,q1,q2;
N=n-1;
if (flag_closed8 && flag_parametric8)
  N=n;
for (i=0; i<N; i++){
  ip1=(i+1) % n;
  H[i][k+2]=p18;
  h=H[i][BE]; h_inv=1.0/h; q0=H[ip1][k]-H[i][k];
  q1=3*h_inv*q0 - 2*p18 -0.5*h*p28;
  q2=6*h_inv*(h_inv*q0-p18)-2*p28;
  p18=q1; p28=q2;
}// for
if (N<n)
  H[n-1][k+2]=p18;
}// multi_step_classical_cubic----------------------------------------------------------------------

int Optimize(double **H,int n,double *Theta,int N,bool flag_coarse,bool *flag_opt=NULL){
// Here we solve for the classical cubic spline.
// The only inputs are H,n and eventually we'll want flag_closed8.
if (n == 0)
    return 0;
double a,b,c,d,e,f,temp;
/*flag_closed8=true;
if (lam8==3)
  flag_parametric8=false;
else
  flag_parametric8=true;*/
time_increments(H,n);
if (flag_closed8){
  p18=0; p28=0; multi_step_classical_cubic(H,n,0); e=p18; f=p28;
  p18=1; p28=0; multi_step_classical_cubic(H,n,0); a=p18-e; c=p28-f;
  p18=0; p28=1; multi_step_classical_cubic(H,n,0); b=p18-e; d=p28-f;
  a=1-a; b=-b; c=-c; d=1-d; temp=1/(a*d-b*c);
  p18=temp*(d*e-b*f); p28=temp*(-c*e+a*f);
  multi_step_classical_cubic(H,n,0);  
  p18=0; p28=0; multi_step_classical_cubic(H,n,1); e=p18; f=p28;
  p18=1; p28=0; multi_step_classical_cubic(H,n,1); a=p18-e; c=p28-f;
  p18=0; p28=1; multi_step_classical_cubic(H,n,1); b=p18-e; d=p28-f;
  a=1-a; b=-b; c=-c; d=1-d; temp=1/(a*d-b*c);
  p18=temp*(d*e-b*f); p28=temp*(-c*e+a*f);
  multi_step_classical_cubic(H,n,1);  
}// if
else{
  p18=0; p28=0; multi_step_classical_cubic(H,n,0); b=p28;
  p18=1; p28=0; multi_step_classical_cubic(H,n,0); a=p28-b;
  p18=-b/a; p28=0; multi_step_classical_cubic(H,n,0);
  p18=0; p28=0; multi_step_classical_cubic(H,n,1); b=p28;
  p18=1; p28=0; multi_step_classical_cubic(H,n,1); a=p28-b;
  p18=-b/a; p28=0; multi_step_classical_cubic(H,n,1);
}// else
return 0;
}// Optimize--------------------------------------------------------------------------------------

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   THE FUNCTIONS BELOW ARE SOMEWHAT GENERIC (ie not method specific)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

int S_canonical(){return 0;}

bool is_feasible_lite (double *u, double *v,double psi,double alpha_max=M_PI_PLUS){return true;}

bool is_feasible (double *u, double *v){return true;}

int S_general(double *u, double *v){return 0;}

int optimal_direction1(double *u,double *v,double *w,double *Theta,bool flag_offset){return 0;}

bool S_spline_feasible(double **H,int n,double *Theta,int N,long Max_trials=1000000){return true;}

int S_spline_check_feasible(double **H,int n){return 0;}

int S_spline_feasible_stencil(double **H,int n,int N,bool flag_check_feasibility=true){return 0;}

void S_spline_rectify(double **H,int n,int N){int i; i=0;}

void S_spline_feasible_stencil3(double **H,int n,int N,int i){int j; j=0;}


}// namespace
