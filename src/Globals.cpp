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


#include <math.h>
#include <QMutex>

// Global variables
double alpha,beta; // read by G_val and sigma_wig
double df_dx; // written to by sigma_wig and entry_angle
double del;   // read by entry_angle
int nargout;  // read by sigma_wig and entry_angle
double a8,b8; // read and written to by bisection_method, bisection_method5, bisection_Newton
int fa8,fb8;  // read and written to by bisection_method, bisection_method5, bisection_Newton
double tol=1e-13;   // read by bisection_method, bisection_method5, bisection_Newton
int N18, N28, N38;  // read by bisection_Newton
int iter8;  // written to by bisection_Newton
double tol2;  // read by bisection_Newton
int Form8;    // written to by S_canonical
double penalty8=0; // read by S_canonical
bool flag_shortcut8=true; // read by S_canonical (if true, experimental shortcut will be used)
int FORM8[4]; // written to by S_general and is_feasible
double BE8,p18,p28, q18, q28; // written to by S_canonical and S_general
double breadth8; // written to by direction_angle and read by optimal_direction1, optimal_direction2
int N8;          // read by optimal_direction1 and optimal_direction2
//int D8;          // read by S_curve_old
//complex c18,c28; // written by S_curve_old
int error_code8; // written by S_curve, E_theta_inv and S_spline_feasible
bool flag_change8; // written by optimal_direction1 and optimal_direction2
long Num_restarts8; // written by S_spline_feasible
double c8,s8,c_wig8,s_wig8; // written and read in connection with cos_sin and cos_sin_inv
double temp8;       // temporary variable for anyone to use
double alpha_max8=M_PI_2; // read by is_feasible_lite whereby alpha > alpha_max8 indicates non-feasibility
double lam8=1;

bool flag_quick8=false; //read by S_general only in restricted elastic spline to indicate optimization on the fly.
double alpha_maxQ8, lamQ8, tolQ8; //used in on the fly optimization (i.e quick)
bool flag_closed8 = false, flag_parametric8 = false;


//Hakims Stuff
QMutex mutex;
