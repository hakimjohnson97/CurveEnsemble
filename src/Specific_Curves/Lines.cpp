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


#include "Lines.h"
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include <QDebug>
#include "Curve_Libs/common_curve_lib.h"
#include "Curve_Libs/lines_lib.h"
using namespace std;
using namespace LinesS;

/*const int X = 0;
const int Y = 1;
const int DIR = 2;
const int CORNER = 3;
const int STATE = 4;
const int OBE = 5;

const int FREE = 7;
const int REGULAR = 0;
const int CLAMPED = 1;
const int FREE_CLAMPED = 3;
const int CLAMPED_FREE = 5; */

//const double INFINITY = 1e+200;

Lines::Lines()
{

    //Set Global variables
    tol = 1*pow(10, -13);
//    D8 = 257;

    //int i;
   // pieces = new piece[MAX_NODES];


    flag_closed = true;
   // selected = -1;
}

/* int Lines::getForm(int index)
{
   // return pieces[index].form;
} */


void Lines::S_splineFeasible()
{
    double theta[coarseRes];
    equi_angular_partition(theta, coarseRes);
    S_spline_feasible(H, count, theta, coarseRes);
}

void Lines::S_splineFeasibleSmart()
{
    S_spline_feasible_stencil(H, count, coarseRes);
}

void Lines::S_splineFeasibleSmart3(int index)
{
    S_spline_feasible_stencil3(H, count, coarseRes, index);
}

void Lines::optimizeCurve(double** nodes, int numNodes, InitDir initDir)
{
    double theta[coarseRes];
    equi_angular_partition(theta, coarseRes); printf("OptimizeCurve: hello \n");
    if (initDir == Smart)
        S_spline_feasible_stencil(nodes, numNodes, coarseRes);
    else if (initDir == Random)
        S_spline_feasible(nodes, numNodes, theta, coarseRes);
    Optimize(nodes, numNodes, theta, coarseRes, true, NULL);
}

void Lines::fineOptimizeCurve(double **nodes, int numNodes, int depth)
{
     int i;
     double h = 360/double(coarseRes);
     int k = 8; printf("fineOptimizeCurve: hello \n");
     double theta[2*k];
     for (i = 0; i < depth; i++)    {
         h /= 4;
         equi_angular_cluster(theta, h, k);
         Optimize(nodes, numNodes, theta, 2*k, false);
     }
}

/*int Lines::save(fstream &file)
{

}

int Lines::load(fstream &file)
{

}

int Lines::saveAsText(fstream &file)
{

}

int Lines::loadFromText(fstream &file)
{

} */

int Lines::S_Curve(double *a, double *b, double *rx, double *ry)
{

    piece p;
  //  if (pieceValid[index] == true)
  //      return S_curve_lite(pieces[index], rx, ry, 257);
    if (make_piece(a, b, &p) == 0)     {
     //   pieces[index] = p;
      //  pieceValid[index] = true;
        return S_curve_lite(p, rx, ry, 257);
    }
    else
        return 0;


}

Curve::Node Lines::pointInTheMiddle(double *a, double *b)
{
    piece p;
    if (make_piece(a, b, &p) == 0)    {
        complex c = point_in_the_middle(p);
        return Curve::Node(c.real, c.imag);
    }
    return Node();
}


double Lines::runCommand(int c, double c2)
{
    c2 = Curve::runCommand(c, c2);
    if (c==0){ // optimize dense data ball curve
       int i,j,k,j_opt,iter=0;
       for (i = 0; i < count; i++)
          for (j = 0; j < 6; j++)   Htemp[i][j] = 0;

       double **B = Htemp;
       int B_count=count; // usage, eg., curve2->count += 1;
       if (H[0][STATE]>100 || H[count-1][STATE]>100 || H[0][STATE]<0 || H[count-1][STATE]<0) return 0;
       j_opt=furthest_marked(); // find furthest marked data ball and unmark it
       while (j_opt>=0 && iter<count){
         iter++;
         H[j_opt][STATE]=100-H[j_opt][STATE]; 
         for (i = 0; i < count; i++) // next we extract unmarked nodes of active curve into curve 2
           for (j = 0; j < 6; j++)  B[i][j]=H[i][j];
         B_count = extract_unmarked(B,count);
         optimizeCurve(B,B_count,AsIs);  //coarse optimize NODEDIR in curve 2
         fineOptimizeCurve(B, B_count, 4);   //fine optimize NODEDIR in curve 2
         k=-1;  // we copy curve 2 into unmarked nodes of active curve
         for (i = 0; i < B_count; i++){
           k++; 
           while(H[k][STATE]>100)  k=(k+1) % count;
           for (j = 0; j < 6; j++)  H[k][j]=B[i][j];
         }// for 
         j_opt=furthest_marked(); // find furthest marked data ball and unmark it
       }// while
       for (i = 0; i < B_count; i++) // next we extract unmarked nodes of active curve into curve 2
         for (j = 0; j < 6; j++)  H[i][j]=B[i][j];
       count=B_count;
       c2=getCurveBE();
    }// if 
    return c2;
}
