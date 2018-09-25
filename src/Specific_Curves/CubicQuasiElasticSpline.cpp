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


#include "CubicQuasiElasticSpline.h"
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include <QDebug>
#include "Curve_Libs/common_curve_lib.h"
#include "Curve_Libs/cubic_quasi_elastic_spline_lib.h"
using namespace std;
using namespace CubicQuasiElasticS;

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

CubicQuasiElasticSpline::CubicQuasiElasticSpline()
{

    //Set Global variables
    tol = 1*pow(10, -13);
//    D8 = 257;

    int i;
   // pieces = new piece[MAX_NODES];


    flag_closed = true;
   // selected = -1;
}

int CubicQuasiElasticSpline::getForm(int index)
{
   // return pieces[index].form;
}


void CubicQuasiElasticSpline::S_splineFeasible()
{
    double theta[coarseRes];
    equi_angular_partition(theta, coarseRes);
    S_spline_feasible(H, count, theta, coarseRes);
}

void CubicQuasiElasticSpline::S_splineFeasibleSmart()
{
    S_spline_feasible_stencil(H, count, coarseRes);
}

void CubicQuasiElasticSpline::S_splineFeasibleSmart3(int index)
{
    S_spline_feasible_stencil3(H, count, coarseRes, index);
}

void CubicQuasiElasticSpline::optimizeCurve(double** nodes, int numNodes, InitDir initDir)
{
    int i;
    double theta[coarseRes];
    equi_angular_partition(theta, coarseRes);
    if (initDir == Smart)
        S_spline_feasible_stencil(nodes, numNodes, coarseRes);
    else if (initDir == Random)
        S_spline_feasible(nodes, numNodes, theta, coarseRes);
    Optimize(nodes, numNodes, theta, coarseRes, true, NULL);
}

void CubicQuasiElasticSpline::fineOptimizeCurve(double **nodes, int numNodes, int depth)
{
     int i;
     double h = 360/double(coarseRes);
     int k = 8;
     double theta[2*k];
     for (i = 0; i < depth; i++)    {
         h /= 4;
         equi_angular_cluster(theta, h, k);
         Optimize(nodes, numNodes, theta, 2*k, false);
     }
}

/*int CubicQuasiElasticSpline::save(fstream &file)
{

}

int CubicQuasiElasticSpline::load(fstream &file)
{

}

int CubicQuasiElasticSpline::saveAsText(fstream &file)
{

}

int CubicQuasiElasticSpline::loadFromText(fstream &file)
{

} */

int CubicQuasiElasticSpline::S_Curve(double *a, double *b, double *rx, double *ry)
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

Curve::Node CubicQuasiElasticSpline::pointInTheMiddle(double *a, double *b)
{
    piece p;
    if (make_piece(a, b, &p) == 0)    {
        complex c = point_in_the_middle(p);
        return Curve::Node(c.real, c.imag);
    }
}


void CubicQuasiElasticSpline::boundingRect(double &x, double &y, double &width, double &height)
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
}
