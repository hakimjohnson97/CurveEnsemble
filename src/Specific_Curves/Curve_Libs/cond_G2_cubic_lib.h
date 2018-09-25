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


namespace CondG2CubicS
{         // See Notebook Dec. 4, 2013

struct piece {

  short form;
  //bool flag_reversal, flag_reflection;
  double p1,p2,q1,q2,alpha,beta; 
  complex c1,c2;                 
};                               

int S_canonical(){
// error_code = S_canonical () 
// Inputs: (global variables) alpha, beta, FORM8[3], tol.
// Outputs: (global variables) BE8, FORM8[0], p18, p28, q18, q28
//
// Error Codes:
// 0. No errors.
// 1. alpha, beta do not satisfy canonical assumptions (see Notebook March 21, 2013)
//
// Description in case FORM8[3]=0:
// Let u=[0,e^(I alpha)], v=[1,e^(I beta)] be a pair of unit tangent vectors in canonical form: 
// |alpha|, |beta| <= pi/3. We find the optimal curve s and
// output Form=1, BE, p1, p2, q1, q2 to FORM8[0], BE8, p18, p28, q18, q28
//
// Description for other cases:
// FORM8[3]=1 indicates that the curve is clamped on the right but free on the left.
// FORM8[3]=2 indicates that the curve is clamped on the left but free on the right.
// FORM8[3]=3 indicates that the curve is free on both sides.  
//
// #################### see Notebook March 21, 2013 ####################
double shape_param;
shape_param=lam8;
if (shape_param>2.99) shape_param=2.99;
if (shape_param<0.01) shape_param=0.01;
FORM8[0]=1;
if (FORM8[3]==1) { // this indicates the curve is clampled on the right but free on left
  cos_sin(beta); alpha=atan((-shape_param*s8)/(3-shape_param*c8));} // if
else if (FORM8[3]==2) { // this indicates the curve is clampled on the left but free on the right
  cos_sin(alpha); beta=atan((-shape_param*s8)/(3-shape_param*c8));} // else if
else if (FORM8[3]==3) { // this is the case when the curve is free on both ends
  alpha=0; beta=0;} // else if
//
// test whether |alpha| <= pi/2
if (fabs(alpha)>M_PI_2+tol || fabs(beta)>M_PI_2+tol)
  {BE8=INFINITY; printf("S_canonical: error_code=1.\n"); return 1;}
cos_sin(beta); q18=c8; q28=s8;
cos_sin(alpha); p18=c8; p28=s8;
temp8=c8*q18 + s8*q28;  // temp8 holds cos(alpha-beta)
BE8=6+2*shape_param*(shape_param*(temp8+2)-3*(c8+q18)); 
return 0;
}// S_canonical--------------------------------------------------------------------------------------------------------

double BE_quick (double *u, double *v,double psi,double recip_breadth){
// be=BE_quick(u,v,psi,recip_breadth) returns the energy of the configuration (u,v) if it is feasible and returns
// INFINITY otherwise.  For sake of speed, we assume that the calling function has pre-computed the direction angle
// psi (from u to v) and the reciprocal of the breadth recip_breadth.
// The global variables max_alphaQ8, lamQ8 and tolQ8 are read and the global variable flag_change8 is written; however
// flag_change8 is used here to indicate whether the configuration is feasible.
bool flag_free=false;
double u2,v2,alf,bet,shape_param,q1,q2,p1,p2,temp;
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
u2=fabs(alf)-tolQ8; v2=fabs(bet)-tolQ8; flag_change8=true;
if (u2 > alpha_maxQ8 || u2>M_PI_2 || v2 > alpha_maxQ8 || v2>M_PI_2 )
  {flag_change8=false; return INFINITY;}
// The above is probably method independent; the rest is method dependent
shape_param=lamQ8;
if (shape_param>2.99) shape_param=2.99;
if (shape_param<0.01) shape_param=0.01;
Cos_Sin(alf,&p1,&p2);
if (flag_free)
  {bet=atan((-shape_param*p2)/(3-shape_param*p1));}
Cos_Sin(bet,&q1,&q2);
temp=p1*q1 + p2*q2;  // temp holds cos(alpha-beta)
temp=6.0+2.0*shape_param*(shape_param*(temp+2.0)-3.0*(p1+q1)); 
return temp*recip_breadth;
}// BE_quick----------------------------------------------------------------------------------------

bool is_feasible_lite (double *u, double *v,double psi, double dummy_var=M_PI_PLUS){
// b=is_feasible_lite(u,v,psi) returns true if the configuration (u,v) is feasible and returns false otherwise.  Here
// u and v are nodes (see README). Additionally, the global variables associated with canonical form are written.
// Specifically, alpha, beta and FORM8[3], where 
// FORM8[3]=0 if the curve is clamped on both sides,
// FORM8[3]=1 if the curve is clamped on the right but free on the left.
// FORM8[3]=2 if the curve is clamped on the left but free on the right.
// FORM8[3]=3 if the curve is free on both sides.  
//
// The optional input dummy_var is not used.  
// After running is_feasible_lite, S_canonical can be run.
bool u_is_free=false, v_is_free=false;
double u2,v2;
FORM8[3]=0; // default options
if (u[STATE]>FREE_CLAMPED) // ie u[STATE] equals CLAMPED_FREE or FREE
  u_is_free=true;
if (v[STATE]==FREE_CLAMPED || v[STATE]==FREE)
  v_is_free=true;
u2=u[DIR]; v2=v[DIR]-v[CORNER]; alpha=0; beta=0;
if (!u_is_free)
  alpha=mod2pi(u2*M_PI_180-psi);
else
  FORM8[3] += 1;
if (!v_is_free)
  beta=mod2pi(v2*M_PI_180-psi);
else
  FORM8[3] += 2;
u2=fabs(alpha)-tol; v2=fabs(beta)-tol;
if (u2 > alpha_max8 || u2>M_PI_2 || v2 > alpha_max8 || v2>M_PI_2 )
  return false;
else
  return true;
}// is_feasible_lite----------------------------------------------------------------------------------------

bool is_feasible (double *u, double *v){
// b=is_feasible(u,v) returns true if the configuration (u,v) is feasible and returns false otherwise.  Here
// u and v are nodes (see README). Additionally, the global variables associated with canonical form are set.
// Specifically, alpha, beta, breadth8 and FORM8[3].
// After running is_feasible, S_canonical can be run.
double psi;
psi=direction_angle(v[X]-u[X],v[Y]-u[Y]);
return is_feasible_lite(u,v,psi);
}// is_feasible--------------------------------------------------------------------------------------------------------

int S_general(double *u, double *v){
// error_code = S_general (*u, *v) does the same thing as S_general, where u and v are two nodes.
// Global variables: (read after calling S_canonical)
//                   BE8, FORM8[3], p18, p28, q18, q28 (written as outputs)
//                   alpha, beta             (written for S_canoncial to read)
// For convenience, we update u[BE] = BE8.
//
// S_general is essentially a front end to S_canonical.  The returned integer error_code actually comes
// from S_canonical and carries the following significance:
// Error Codes:
// 0. No errors.
// 1. alpha, beta do not satisfy canonical assumptions.
int error_code;
double breadth;
if (is_feasible(u,v))
  {breadth=breadth8; error_code=S_canonical(); BE8=BE8/breadth; u[BE]=BE8; return error_code;}
else
  {BE8=INFINITY; u[BE]=BE8; return 1;}
}// S_general ----------------------------------------------------------------------------------------------

int make_piece (double *u,double *v,piece *p){
// error_code = make_piece(u,v,p) makes the piece p representing the distinguished curve
// connecting node u to node v.  The returned value error code is the result of calling S_general(u,v).
// See Notebook March 23, 2013
double shape_param;
shape_param=lam8;
if (shape_param>2.99) shape_param=2.99;
if (shape_param<0.01) shape_param=0.01;
int error_code;
complex c1,c2,w1,w2,z1,z2;
piece q;
error_code = S_general(u,v);
if (error_code==0){
if (FORM8[3]==3)
  {q.form = (short) 0; q.p1=0; q.p2=0; q.q1=0; q.q2=0; q.alpha=0; q.beta=0;}
else {
  q.form = (short) 1; q.p1 = shape_param*p18; q.p2 = shape_param*p28; 
  q.q1=shape_param*q18; q.q2=shape_param*q28; q.alpha=alpha; q.beta=beta;
} // else
w1.real=u[X]; w1.imag=u[Y]; w2.real=v[X]; w2.imag=v[Y];
z1.real=0; z2.real=1; z1.imag=0; z2.imag=0;
c1=complex_div( complex_sub(w2,w1)  ,  complex_sub(z2,z1)  );
c2=complex_sub( w1   ,  complex_mult(c1,z1)  );
q.c1=c1; q.c2=c2;
*p=q;
}// if error_code==0
return error_code;
}// make_piece---------------------------------------------------------------------------------

int S_curve_lite (piece p,double *x, double *y,int D=257){
// k = S_curve_lite(p,x,y,D)
// Given a piece p, we return k points (x,y) along the curve represented by p.
// The density of points is determined by the input D, where the guiding principal is that a 
// u-turn will receive D points.  Regardless of all else, we maintain 2 <= k <= 2D.
// See Notebook March 23, 2012 for more details.
int k,j;
double h,t,b1,b2,b3;
complex z;
//
// We first construct the curve K and store it in (x,y)
if (p.form==0)
  {x[0]=0; y[0]=0; x[1]=1; y[1]=0; k=2;}
else {
  k=int(double(D)*(fabs(p.alpha)+fabs(p.beta))*M_1_PI); 
  if (k<2)
     k=2;
  if (k>2*D)
     k=2*D;
  h=1/double(k-1); t=0;
  for (j=0; j<k; j++)
    { b1=(-2*t+3)*t*t;
      b2=(1-t)*(1-t)*t;
      b3=(t-1)*t*t;
      x[j]=b1 + b2*p.p1 + b3*p.q1;
      y[j]=b2*p.p2 + b3*p.q2;
      t+=h;}
}// else
//
// Now we apply the transformation T(z) = c1 z + c2 or
for (j=0; j<k; j++){
    z.real=x[j]; z.imag=y[j];
    z=complex_add(complex_mult(p.c1,z),p.c2); x[j]=z.real; y[j]=z.imag;
}// for
return k;
}// S_curve_lite---------------------------------------------------------------------------------

complex point_in_the_middle (piece p,double t=0.5){
// point_in_the_middle(p,t)
// Given a piece p, we return a complex point z, roughly at position t along the curve, where
// 0<=t<=1.
double b1,b2,b3,x,y;
complex z;
//
// We first construct the point in canonical form
if (p.form==0)
  {x=t; y=0;}
else {b1=(-2*t+3)*t*t;
      b2=(1-t)*(1-t)*t;
      b3=(t-1)*t*t;
      x=b1 + b2*p.p1 + b3*p.q1;
      y=b2*p.p2 + b3*p.q2;
}// else
//
// Now we apply the transformation T(z) = c1 z + c2 or
    z.real=x; z.imag=y;
    z=complex_add(complex_mult(p.c1,z),p.c2);
return z;
}// point_in_the_middle---------------------------------------------------------------------------------


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   THE FUNCTIONS BELOW ARE SOMEWHAT GENERIC (ie not method specific)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
// double tolQ8; tolQ8=tol;
be_opt = u[BE]+v[BE]; v2_save=v[DIR]; v2_opt=v[DIR]; error_code=0;
//psi_u=direction_angle(v[X]-u[X],v[Y]-u[Y]); breadth_u_inv=1/breadth8;
x=v[X]-u[X]; y=v[Y]-u[Y]; breadth_u_inv=1/(sqrt(x*x+y*y)+1e-200);
psi_u=Cos_Sin_inv(breadth_u_inv*x,breadth_u_inv*y);
//psi_v=direction_angle(w[X]-v[X],w[Y]-v[Y]); breadth_v_inv=1/breadth8;
x=w[X]-v[X]; y=w[Y]-v[Y]; breadth_v_inv=1/(sqrt(x*x+y*y)+1e-200);
psi_v=Cos_Sin_inv(breadth_v_inv*x,breadth_v_inv*y);
for (k=0;k<N8;k++){
  if (flag_offset)
    v[DIR]=mod360(v2_save+Theta[k]);
  else
    v[DIR]=Theta[k];
  be_u=BE_quick(u,v,psi_u,breadth_u_inv);
  if (flag_change8 && be_u<be_opt){
    be_v=BE_quick(v,w,psi_v,breadth_v_inv);
        if (flag_change8 && be_u+be_v+tolQ8/2<be_opt){
          u[BE]=be_u; v[BE]=be_v; be_opt=be_u+be_v; v2_opt=v[DIR];
        }// if
  }// if
}// for (k=0
v[DIR]=v2_opt;
if (v2_opt==v2_save)
  flag_change8=false;
else
  flag_change8=true;
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
for (i = 0; i < n; i++){
  ip1=(i+1) % n;
  psi[i] = direction_angle(H[ip1][X]-H[i][X], H[ip1][Y]-H[i][Y]);
}// for
// first we initialize to mean stencil directions
for (i=0; i<n; i++){
  if (H[i][STATE]==REGULAR){
    ip1=(i+1) % n; im1=(i+n-1) % n;
    psi0=psi[i];
    phi=direction_angle(H[i][X]-H[im1][X],H[i][Y]-H[im1][Y])+M_PI_180*H[i][CORNER];
    phi=phi-psi0; MOD2PI(phi);
    psi0+=0.5*phi; // Now psi is the desired direction angle in radians
                  // but we want to convert it to an angle in degrees
    psi0=M_180_PI*psi0;
    H[i][DIR]=mod360(psi0);
  }// if
}// for
// next we search for feasible random offsets
for (i = 0; i < n; i++){
  if (H[i][STATE] == REGULAR){
    ip1=(i+1) % n; im1=(i+n-1) % n;
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
  }// if
}// for
// Finally we test feasibility of each piece and update bending energies
feasible = true;
for (i = 0; i < n; i++){
  ip1=(i+1) % n;
  if (is_feasible_lite(H[(i)], H[ip1], psi[i],alpha_max8))
    error_code += S_general(H[i], H[ip1]);  // this command updates the bending energies
  else
    feasible=false;
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
// nodes contained in H and N is leftover from elastic splines, so ignore it.
// Outgoing direction angles at regular nodes of H are chosen to equal the
// mean stencil direction.
// Afterwards, k=S_spline_check_feasible(double **H,int n) is run (this updates the
// bending energies if all pieces are feasible).  However, if you want to skip that, then simply
// give the command S_spline_feasible_stencil(H,n,N,false)
// 
int i,im1,ip1;
double psi,phi;
for (i=0; i<n; i++){
  if (H[i][STATE]==REGULAR){
    ip1=(i+1) % n; im1=(i+n-1) % n;
    psi=direction_angle(H[ip1][X]-H[i][X],H[ip1][Y]-H[i][Y]);
    phi=direction_angle(H[i][X]-H[im1][X],H[i][Y]-H[im1][Y])+M_PI_180*H[i][CORNER];
    phi=phi-psi; MOD2PI(phi);
    psi+=0.5*phi; // Now psi is the desired direction angle in radians
                  // but we want to convert it to an angle in degrees
    psi=M_180_PI*psi;
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
i=(i+n-1) % n;
for (k=0; k<3; k++){
  if (H[i][STATE]==REGULAR){
    ip1=(i+1) % n; im1=(i+n-1) % n;
    psi=direction_angle(H[ip1][X]-H[i][X],H[ip1][Y]-H[i][Y]);
    phi=direction_angle(H[i][X]-H[im1][X],H[i][Y]-H[im1][Y])+M_PI_180*H[i][CORNER];
    phi=phi-psi; MOD2PI(phi);
    psi+=0.5*phi; // Now psi is the desired direction angle in radians
                  // but we want to convert it to an angle in degrees
    psi=M_180_PI*psi;
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
//
// NOTE: flag_offset was motivated by elastic splines where u-turn detection was of critical importance.
// However, for restricted elastic splines and cubic splines Mike wants flag_offset=true regardless.
// Rather than bother Hakim with this, we'll simply set flag_offset=true on line 4 and be done with it.
int i,im1,ip1,error_code=0;
N8=N;
bool all_optimal = false, flag_offset, flag_optimal[n];
flag_offset=true; // originally it was flag_offset=!flag_coarse;
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
      flag_optimal[i] = true;
      if (flag_quick8)
      //if (1==2)
        error_code+=optimal_direction1_quick(H[im1], H[i], H[ip1], Theta, flag_offset);
      else
        error_code+=optimal_direction1(H[im1], H[i], H[ip1], Theta, flag_offset);
      if (flag_change8)    {
          all_optimal = false;
          flag_optimal[im1] = false;
          flag_optimal[i]=flag_coarse;
          flag_optimal[ip1] = false;
      }// if
    }// if
  }// for
}// while
if (flag_opt!=NULL)
  for (i = 0; i < n; i++)
    flag_opt[i]=flag_optimal[i];
return error_code;
}// Optimize--------------------------------------------------------------------------------------


}
