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


namespace LinesS
{         // See Notebook Dec. 4, 2013

struct piece {

  short form;
  //bool flag_reversal, flag_reflection;
  double p1,p2,q1,q2,alpha,beta; 
  complex c1,c2;                 
};

bool is_closed (double **H){
// We'll make an educated guess as to whether of not H is closed or not.
// Since H is a curve of type lines, it should be the case that
if (H[0][STATE]<=0 || H[0][STATE]>100)
  return true;
else
  return false;
}// is_closed----------------------------------------------------------------------------------------------

void BE_stencil(double **H,int n){
// m=data_balls_2_lines(H,n,B) : Given a curve matrix B of type data_balls containing n data balls,
// we convert this to a curve matrix H of type lines, where the fields X,Y are computed and the BE field
// corresponds to the absolute stencil angle.
// The remaining fields are set to 0.
int i,ip1,im1,i_first,i_last;
bool flag_closed=is_closed(H);
double psi_u,psi_v,u_save[6],v_save[6],w_save[6];
if (flag_closed)
  {i_first=0; i_last=n-1;}
else
  {i_first=1; i_last=n-2; H[0][BE]=0; H[n-1][BE]=0;}
for (i=i_first; i<=i_last; i++){
  ip1=(i+1) % n;
  im1=(i+n-1) % n;
  databall2node(H[im1],u_save); databall2node(H[i],v_save); databall2node(H[ip1],w_save);
  psi_u=direction_angle(H[i][X]-H[im1][X],H[i][Y]-H[im1][Y]);
  psi_v=direction_angle(H[ip1][X]-H[i][X],H[ip1][Y]-H[i][Y]);
  H[i][BE]=fabs(mod2pi(psi_v-psi_u));
  node2databall(H[im1],u_save); node2databall(H[i],v_save); node2databall(H[ip1],w_save);
}// for
}// BE_stencil----------------------------------------------------------------------------------

double BE_lines(double *u,double *v){
// be=BE_lines(u,v) : Let u and v be two data balls.  Then be is the distance from the node of u
// to the node of v.
double u1[6],u_save[6],v1[6];
for (int i=0; i<6; i++){
  u1[i]=u[i]; v1[i]=v[i];
}//for
databall2node(u1,u_save); databall2node(v1,u_save); 
return sqrt((u1[0]-v1[0])*(u1[0]-v1[0])+(u1[1]-v1[1])*(u1[1]-v1[1]));
}// BE_lines--------------------------------------------------------------------------------------

int S_general(double *u,double *v){
// S_general(u,v) : Let u and v be two data balls.  Then BE is the distance from the node of u
// to the node of v.  We set u[BE]=BE8=BE accordingly and then return 0.
BE8=BE_lines(u,v); u[BE]=BE8;
return 0;
}// S_general--------------------------------------------------------------------------------------

int make_piece (double *u,double *v,piece *p){
// error_code = make_piece(u,v,p) makes the piece p representing the distinguished curve
// connecting node u to node v.  The returned value error code is the result of calling S_general(u,v).
// See Notebook March 23, 2013
double u1[6],u_save[6],v1[6];
piece q;
for (int i=0; i<6; i++){
  u1[i]=u[i]; v1[i]=v[i];
}//for
databall2node(u1,u_save); databall2node(v1,u_save); 
q.p1=u1[0]; q.p2=u1[1]; q.q1=v1[0]; q.q2=v1[1];
*p=q;
//BE8=sqrt((u1[0]-v1[0])*(u1[0]-v1[0])+(u1[1]-v1[1])*(u1[1]-v1[1]));
//u[BE]=BE8;
return 0;
}// make_piece---------------------------------------------------------------------------------

int S_curve_lite (piece p,double *x, double *y,int D=257){
// k = S_curve_lite(p,x,y,D)
// Given a piece p, we return k points (x,y) along the curve represented by p.
// The density of points is determined by the input D, where the guiding principal is that a 
// u-turn will receive D points.  Regardless of all else, we maintain 2 <= k <= 2D.
// See Notebook March 23, 2012 for more details.
x[0]=p.p1; y[0]=p.p2; x[1]=p.q1; y[1]=p.q2;
return 2;
}// S_curve_lite---------------------------------------------------------------------------------

complex point_in_the_middle (piece p,double t=0.5){
// point_in_the_middle(p,t)
// Given a piece p, we return a complex point z, roughly at position t along the curve, where
// 0<=t<=1.
complex z;
z.real=(1-t)*(p.p1+p.q1); z.imag=t*(p.p2+p.q2);
return z;
}// point_in_the_middle---------------------------------------------------------------------------------

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
double u1[6],u_save[6],v1[6];
for (int i=0; i<6; i++){
  u1[i]=u[i]; v1[i]=v[i];
}//for
databall2node(u1,u_save); databall2node(v1,u_save); 
x[0]=u1[0]; y[0]=u1[1]; x[1]=v1[0]; y[1]=v1[1]; 
return 2;
}//S_curve----------------------------------------------------------------------------------------


bool S_spline_feasible(double **H,int n,double *Theta,int N,long Max_trials=1000000){
return true;
}// S_spline_feasible-------------------------------------------------------------------

int S_spline_check_feasible(double **H,int n){
// k=S_spline_check_feasible(H,n) : The input n equals the number of
// nodes contained in H.  If all pieces of H are feasible, then k=0; 
// otherwise piece k was found to be not feasible 
// NOTE: piece k connects node k-1 to node k.
// If all pieces are feasible, then we update the bending energies, where
// the global variable error_code8 will hold the accumulated error_code outputs
// of calls to S_general (error_code8=0 indicates everything went well).
BE_stencil(H,n);
return 0;
}// S_spline_check_feasible-------------------------------------------------------------------

int S_spline_feasible_stencil(double **H,int n,int N,bool flag_check_feasibility=true){
if (flag_check_feasibility)
  return S_spline_check_feasible(H,n);
else
  return 0;
}// S_spline_feasible_stencil--------------------------------------------------------------------------------------

void S_spline_rectify(double **H,int n,int N){
}// S_spline_rectify--------------------------------------------------------------------------------------

void S_spline_feasible_stencil3(double **H,int n,int N,int i){
}// S_spline_feasible_stencil3-------------------------------------------------------------------------------------

bool optimal_direction1_lines(double *u,double *v,double *w,double *Theta,int i,int ip1,bool flag_closed){
// flag_change = optimal_direction1_lines(u,v,w,Theta,flag_offset) (something like that) :
// Let u,v,w be three consecutive data balls.
// Theta should be an array of angles to be considered as candidates for v[NODEDIR].  
// This procedure finds the optimal choice of v[NODEDIR] and updates u,v,w accordingly.
// Global variables N8, tol are read, while the output flag_change indicates whether
// v[NODEDIR] has changed (true means it changed).
// The calling routine is responsible to set tol (eg 1e-13) and N8 = length of Theta.
double be_opt,v2_opt,be_u,be_v,v2_save;
int k;
be_opt = u[BE]+v[BE]; v2_save=v[NODEDIR]; v2_opt=v[NODEDIR];
for (k=0;k<N8;k++){
  v[NODEDIR]=mod360(v2_save+Theta[k]);
  if (flag_closed || i>0)  be_u=BE_lines(u,v);
  else be_u=0;
  if (be_u<be_opt){
    if (flag_closed || ip1>0) be_v=BE_lines(v,w);
    else be_v=0;
    if (be_u+be_v+tol/2<be_opt){
      u[BE]=be_u; v[BE]=be_v; be_opt=be_u+be_v; v2_opt=v[NODEDIR];
    }// if
  }// if
}// for (k=0
v[NODEDIR]=v2_opt;
if (v2_opt==v2_save)
  return false;
else
  return true;
}// optimal_direction1_lines--------------------------------------------------------------------------------------

void Optimize(double **H,int n,double *Theta,int N,bool flag_coarse,bool *flag_opt=NULL){
// Optimize_lines(H,n,Theta,N,flag_coarse) : The input n equals the number of data balls
// contained in H, and the input N equals then number of degree angles contained
// in Theta. It is not assumed that the bending energies in H are accurate; they are
// computed fresh before entering the main loop.
// Angles in Theta are treated as offset candidates for outgoing
// direction angles in H[i][NODEDIR].  If Theta runs all the way around the compass, then use
// flag_coarse==true; while if Theta is a small cluster near 0 degrees (fine optimization),
// use flag_coarse==false.
int i,im1,ip1;
N8=N;
bool all_optimal = false, flag_change, flag_optimal[n];
bool flag_closed=is_closed(H);
for (i = 0; i < n; i++){
  flag_optimal[i] = false;
  ip1 = (i+1) % n;
  H[i][BE]=BE_lines(H[i],H[ip1]);
  } // for
while (!all_optimal)    {
  all_optimal = true;
  for (i = 0; i < n; i++){
    if (H[i][STATE]<0 && !flag_optimal[i]){
      im1 = (n + i-1) % n;
      ip1 = (i+1) % n;
      flag_optimal[i] = true;
      flag_change=optimal_direction1_lines(H[im1], H[i], H[ip1], Theta,i,ip1,flag_closed);
      if (flag_change)    {
          all_optimal = false;
          flag_optimal[im1] = false;
          flag_optimal[i]=flag_coarse;
          flag_optimal[ip1] = false;
      }// if
    }// if
  }// for
}// while
// Now we fill H[*][BE] with absolute stencil angles
BE_stencil(H,n);
}// Optimize--------------------------------------------------------------------------------------

} // namespace
