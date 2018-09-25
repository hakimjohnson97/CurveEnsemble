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


#include "ResElasticSpline.h"
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include <QDebug>
#include "Curve_Libs/common_curve_lib.h"
#include "Curve_Libs/restricted_elastic_spline_lib.h"
using namespace std;
using namespace ResElasticS;



ResElasticSpline::ResElasticSpline()
{

    //Set Global variables
    tol = 1*pow(10, -13);


    int i;
   // pieces = new piece[MAX_NODES];



    flag_closed = true;
   // selected = -1;
}

int ResElasticSpline::getForm(int index)
{
   // return pieces[index].form;
}


void ResElasticSpline::S_splineFeasible()
{
    double theta[coarseRes];
    equi_angular_partition(theta, coarseRes);
    S_spline_feasible(H, count, theta, coarseRes);
}

void ResElasticSpline::S_splineFeasibleSmart()
{
    S_spline_feasible_stencil(H, count, coarseRes);
   // qDebug() << "g";
}

void ResElasticSpline::S_splineFeasibleSmart3(int index)
{
    S_spline_feasible_stencil3(H, count, coarseRes, index);
   // qDebug() << "f";
}

void ResElasticSpline::optimizeCurve(double **nodes, int numNodes, InitDir initDir)
{
    int i;
    double theta[coarseRes], GAMMA[coarseResGamma];
    equi_angular_partition(theta, coarseRes);
    equi_angular_partition(GAMMA, coarseResGamma);
    if (initDir == Smart)
        S_spline_feasible_stencil(nodes, numNodes, coarseRes);
    else if (initDir == Random)
        S_spline_feasible(nodes, numNodes, theta, coarseRes);
    S_spline_check_feasible(nodes, numNodes);        
    Optimize(nodes, numNodes, theta, coarseRes, true, NULL, GAMMA, coarseResGamma);
}

void ResElasticSpline::fineOptimizeCurve(double** nodes, int numNodes, int depth)
{
     int i;
     double h = 360/double(coarseRes), h_G = 360/double(coarseResGamma);
     int k = 8, k_G = 8;
     double theta[2*k+1], GAMMA[2*k_G];
     S_spline_check_feasible(nodes, numNodes);
     for (i = 0; i < depth; i++)    {
         h /= 4; h_G /= 4;
         equi_angular_cluster(theta, h, k); theta[2*k]=0;
         equi_angular_cluster(GAMMA, h_G, k_G);
         Optimize(nodes, numNodes, theta, 2*k+1, false, NULL, GAMMA, 2*k_G);
     }// for
}

/*int ResElasticSpline::save(fstream &file)
{

}

int ResElasticSpline::load(fstream &file)
{

}

int ResElasticSpline::saveAsText(fstream &file)
{

}

int ResElasticSpline::loadFromText(fstream &file)
{

} */

int ResElasticSpline::S_Curve(double *a, double *b, double *rx, double *ry)
{

    int temp = penalty8;
    penalty8 = 360; //Big number to relieve restriction

    piece p;
  //  if (pieceValid[index] == true)
  //      return S_curve_lite(pieces[index], rx, ry, 257);
    if (make_piece(a, b, &p) == 0)     {
     //   pieces[index] = p;
       // pieceValid[index] = true;
        return S_curve_lite(p, rx, ry, 257);
    }
    else
        return 0;

    penalty8 = temp;

}

Curve::Node ResElasticSpline::pointInTheMiddle(double *a, double *b)
{
    piece p;
    if (make_piece(a, b, &p) == 0)    {
        complex c = point_in_the_middle(p);
        return Curve::Node(c.real, c.imag);
    }
    else
        return NULL;
}


double ResElasticSpline::runCommand(int c, double c2)
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
         mark_all(B,B_count); // mark all databalls in curve 2
         flag_quick8=false; // printf("B_count = %d \n",B_count);
         optimizeCurve(B,B_count,Smart);  //coarse optimize DIR in curve 2
         fineOptimizeCurve(B, B_count, 4);   //fine optimize DIR in curve 2
         unmark_all(B,B_count); // unmark all databalls in curve 2
         flag_quick8=true; tolQ8= tol; alpha_maxQ8= alpha_max8;
         optimizeCurve(B,B_count,AsIs); //coarse optimize DIR and NODEDIR in curve 2
         //fineOptimizeCurve(B, B_count, 1);   //fine optimize DIR and NODEDIR in curve 2
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
       c2=2*sqrt(getCurveBE()*getCurveLength());
    }// if 
    return c2;
}

//Unused Functions

/*void ResElasticSpline::boundingRect(double &x, double &y, double &width, double &height)
{
    x = H[count-1][X]; y = H[count-1][Y]; width = 0; height = 0;
    int i, j;
    double pointsX[coarseRes*2+2];
    double pointsY[coarseRes*2+2];

    int numPieces;
    if (closed)
        numPieces = count;
    else
        numPieces = count-1;
   // for (i = 0; i < numPieces; i++)    {
          i = 0;
        pointsOnCurve(i, pointsX, pointsY);
        for (j = 0; j < coarseRes*2+2; j++)    {
            if (pointsX[j] < x)
                x = pointsX[j];
            if (pointsY[j] < y)
                y = pointsY[j];
            if (pointsX[j] > width)
                width = pointsX[j];
            if (pointsY[j] > height)
                height = pointsY[j];
        }
   // }
}*/

