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


#ifndef COMMON_CURVE_LIB_H
#define COMMON_CURVE_LIB_H

#include <math.h>

//          Constants:
//#define M_PI            3.14159265358979323846264338327950288
#define M_PI_PLUS       3.1415926535898932
#define M_PI_MINUS      3.1415926535896932
//#define M_1_PI          0.31830988618379067153776752674502872
#define M_2PI           6.28318530717958647692528676655900576
//#define M_PI_2          1.57079632679489661923132169163975144
#define M_PI_3_PLUS     1.04719755119661     
#define M_PI_180        0.0174532925199432957692369076848861271
#define M_180_PI        57.2957795130823208767981548141051704
#define M_1_SQRT2       0.70710678118654752440084436210484904
#define M_D             1.19814023473559220743992249228032388
#define M_D_SQ          1.43554002209225999564238644733155885
#define M_D_PI          0.381379881750906594031162990481807897
//
#define M_PI_32         0.0981747704246810387019576057274844651311615437
#define M_SIN_PI_32     0.0980171403295606019941955638886418458611366732
#define M_COS_PI_32     0.995184726672196886244836953109479921575474869
#define M_PI_16         0.196349540849362077403915211454968930262323087
#define M_SIN_PI_16     0.195090322016128267848284868477022240927691618
#define M_COS_PI_16     0.980785280403230449126182236134239036973933731
#define M_3_PI_32       0.294524311274043116105872817182453395393484631
#define M_SIN_3_PI_32   0.290284677254462367636192375817395274691476278
#define M_COS_3_PI_32   0.956940335732208864935797886980269969482849206
#define M_PI_8          0.392699081698724154807830422909937860524646175
#define M_SIN_PI_8      0.382683432365089771728459984030398866761344562
#define M_COS_PI_8      0.923879532511286756128183189396788286822416626
#define M_5_PI_32       0.490873852123405193509788028637422325655807719
#define M_SIN_5_PI_32   0.471396736825997648556387625905254377657460319
#define M_COS_5_PI_32   0.881921264348355029712756863660388349508442621
#define M_3_PI_16       0.589048622548086232211745634364906790786969262
#define M_SIN_3_PI_16   0.555570233019602224742830813948532874374937191
#define M_COS_3_PI_16   0.831469612302545237078788377617905756738560812
#define M_7_PI_32       0.687223392972767270913703240092391255918130806
#define M_SIN_7_PI_32   0.634393284163645498215171613225493370675687095
#define M_COS_7_PI_32   0.773010453362736960810906609758469800971041293
//#define M_PI_4          0.78539816339744830961566084581987572104929235
#define M_SIN_PI_4      0.707106781186547524400844362104849039284835938
#define M_COS_PI_4      0.707106781186547524400844362104849039284835938
#define M_9_PI_32       0.883572933822129348317618451547360186180453894
#define M_SIN_9_PI_32   0.773010453362736960810906609758469800971041293
#define M_COS_9_PI_32   0.634393284163645498215171613225493370675687095
#define M_5_PI_16       0.981747704246810387019576057274844651311615437
#define M_SIN_5_PI_16   0.831469612302545237078788377617905756738560812
#define M_COS_5_PI_16   0.555570233019602224742830813948532874374937191
#define M_11_PI_32      1.07992247467149142572153366300232911644277698
#define M_SIN_11_PI_32  0.881921264348355029712756863660388349508442621
#define M_COS_11_PI_32  0.471396736825997648556387625905254377657460319
#define M_3_PI_8        1.17809724509617246442349126872981358157393852
#define M_SIN_3_PI_8    0.923879532511286756128183189396788286822416626
#define M_COS_3_PI_8    0.382683432365089771728459984030398866761344562
#define M_13_PI_32      1.27627201552085350312544887445729804670510007
#define M_SIN_13_PI_32  0.956940335732208864935797886980269969482849206
#define M_COS_13_PI_32  0.290284677254462367636192375817395274691476278
#define M_7_PI_16       1.37444678594553454182740648018478251183626161
#define M_SIN_7_PI_16   0.980785280403230449126182236134239036973933731
#define M_COS_7_PI_16   0.195090322016128267848284868477022240927691618
#define M_15_PI_32      1.47262155637021558052936408591226697696742316
#define M_SIN_15_PI_32  0.995184726672196886244836953109479921575474869
#define M_COS_15_PI_32  0.0980171403295606019941955638886418458611366732

//         Node indeces
#define X 0
#define Y 1
#define DIR 2
#define RADIUS 2 // for data balls only (depricated)
#define CORNER 3
#define NODEDIR 3 // for data balls only
#define STATE 4 
#define STATUS 4  // synonym for STATE
#define TANDIR 4  // for data balls only (depricated)
#define BE 5
//         STATE values
#define REGULAR 0
#define CLAMPED 1
#define FREE_CLAMPED 3
#define CLAMPED_FREE 5
#define FREE 7
/*         Changes:
If u[STATE]<0, then u is a data ball with 
   center = u[X]+I*u[Y] and radius = -u[STATE].  The actual node is
   center + radius * exp(I*u[NODEDIR]*M_PI_180), with corner = 0, STATE = 0 and
   u[DIR] and u[BE] as is  */

// Global variables
extern double alpha,beta; // read by G_val and sigma_wig
extern double df_dx; // written to by sigma_wig and entry_angle
extern double del;   // read by entry_angle
extern int nargout;  // read by sigma_wig and entry_angle
extern double a8,b8; // read and written to by bisection_method, bisection_method5, bisection_Newton
extern int fa8,fb8;  // read and written to by bisection_method, bisection_method5, bisection_Newton
extern double tol;   // read by bisection_method, bisection_method5, bisection_Newton
extern int N18, N28, N38;  // read by bisection_Newton
extern int iter8;  // written to by bisection_Newton
extern double tol2;  // read by bisection_Newton
extern int Form8;    // written to by S_canonical
extern double penalty8; // read by S_canonical
extern bool flag_shortcut8; // read by S_canonical (if true, experimental shortcut will be used)
extern int FORM8[4]; // written to by S_general and is_feasible
extern double BE8,p18,p28,q18,q28; // written to by S_canonical and S_general
extern double breadth8; // written to by direction_angle and read by optimal_direction1, optimal_direction2
extern int N8;          // read by optimal_direction1 and optimal_direction2
extern int error_code8; // written by S_curve, E_theta_inv and S_spline_feasible
extern bool flag_change8; // written by optimal_direction1 and optimal_direction2
extern long Num_restarts8; // written by S_spline_feasible
extern double c8,s8,c_wig8,s_wig8; // written and read in connection with cos_sin and cos_sin_inv
extern double temp8;       // temporary variable for anyone to use
extern double alpha_max8; // read by is_feasible_lite whereby alpha > alpha_max8 indicates non-feasibility
extern double lam8;       // read by Conditionally G2 Cubic and C2_par_cubic
extern bool flag_quick8; //read by S_general only in restricted elastic spline to indicate optimization on the fly.
extern bool flag_closed8; // indicates a closed curve
extern bool flag_parametric8; // read by C2_par_cubic to decide when to do nonparametric classical cubic spline.
extern double alpha_maxQ8, lamQ8, tolQ8; //used in on the fly optimization (i.e quick)


//         Macros:
//  If a belongs to (-3pi,3pi], then after MOD2PI(a), a will belong to (-pi,pi] and
//  will be equivalent, mod 2pi, to the original value of a.  Apply MOD2PI(a) twice
//  if a belongs to (-5pi,5pi], etc.
#define MOD2PI(a) if (a>M_PI) a-=M_2PI; else if (a<= -M_PI) a+=M_2PI;

inline double mod2pi(double a){
// a1=mod2pi(a) returns a1 in (-pi,pi] which is equivalent to a, modulo 2pi.
if (a>M_PI)
  while (a>M_PI)
    a-=M_2PI;
else
  while (a<= -M_PI)
    a+=M_2PI;
return a;
}// mod2pi--------------------------------------------------------------------------------------------------------

inline int sign(double x){
if (x>0)
  return 1;
else if (x<0)
  return -1;
else
  return 0;
}// sign----------------------------------------------------------------------------------------------


double signum(double x);
double mod360(double a);
void Cos_Sin(double t,double *c,double *s);
double Cos_Sin_inv(double c,double s);
double cos_sin_inv ();
void equi_angular_partition(double *Theta,int N);
void equi_angular_cluster(double *Theta,double h,int k);
//inline double direction_angle(double x,double y);
//inline int sign(double x);
double bisection_method(double (*f) (double)  ,double target);
double bisection_method5(double (*f) (double)  ,double target);
double bisection_Newton(double (*f) (double)  ,double target);

inline void cos_sin(double t){
//  Assuming t in [-pi,3pi], cos_sin(t) sets global variables c8=cos(t) and s8=sin(t).
Cos_Sin(t,&c8,&s8);
}// cos_sin

inline double cos_sin_inv (){
// t=cos_sin_inv() returns t in (-pi,pi] such that cos(t)=c8 and sin(t)=s8,
// where c8 and s8 are global variables.
return Cos_Sin_inv(c8,s8);
}// cos_sin_inv -----------------------------------------------------------------------

inline double direction_angle(double x,double y){
  // alf=direction_angle(x,y) sets breadth8 = sqrt(x^2 + y^2) and returns alf = arg(x+iy).
  // In order to avoid divisions by zero, breadth8 is not allowed to fall below 1e-200.
  // After exiting, I leave c8=cos(alf) and s8=sin(alf).
  breadth8=sqrt(x*x+y*y)+1e-200;
  temp8=1/breadth8; c8=temp8*x; s8=temp8*y;
  return Cos_Sin_inv (c8,s8);
}// direction_angle--------------------------------------------------------------------------------

inline void databall2node(double *u, double *u_save){
// databall2node(u,u_save) : Given a data ball or node u, we first copy u[STATE] and u[BE]
// into u_save.
// If u[STATE]<0 or >100, we then copy the rest of u into u_save and then  convert u into its node 
// (with STATE=0 and CORNER=0).  Afterwards, you can see whether u has been converted by checking
// whether u_save[STATE]<0 or >100
u_save[STATE]=u[STATE]; u_save[BE]=u[BE];
if (u[STATE]>100)
   u[STATE]=100-u[STATE];
if (u[STATE]<0){
  u_save[0]=u[0]; u_save[1]=u[1]; u_save[2]=u[2];
  u_save[3]=u[3]; 
  cos_sin(u[NODEDIR]*M_PI_180);
  u[X]=u[X]-u[STATE]*c8; u[Y]=u[Y]-u[STATE]*s8;
  u[CORNER]=0; u[STATE]=0;
}// if
}// databall2node--------------------------------------------------------------------------------------

inline void node2databall(double *u, double *u_save){
// node2databall(u,u_save) : If u_save[STATE]<0, then we copy 
// u_save[X,Y,NODEDIR,STATE] into u
if (u_save[STATE]<0 || u_save[STATE]>100)
  {u[X]=u_save[X]; u[Y]=u_save[Y]; u[NODEDIR]=u_save[NODEDIR]; u[STATE]=u_save[STATE];}
}// node2databall--------------------------------------------------------------------------------------

struct complex {double real,imag;};

complex complex_create(double x, double y);
complex complex_add(complex x,complex y);
complex complex_sub(complex x,complex y);
complex complex_mult(complex x,complex y);
complex complex_mult(double t,complex x);
complex complex_conj(complex x);
complex complex_inv(complex x);
double complex_abs(complex x);
complex complex_polar(complex x);
complex complex_div(complex x,complex y);

void mark_all(double **H,int n);
void unmark_all(double **H,int n);
int extract_unmarked(double **H,int n);
double closest_point_thresh(double *w,double *z,double *x,double *y,int n,double thresh);
void anchor_data_balls(double **H,int n);
double max_radius(double **H,int n);



#endif // COMMON_CURVE_LIB_H
