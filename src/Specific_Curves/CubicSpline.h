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


#ifndef CUBICSPLINE_H
#define CUBICSPLINE_H

#include <cmath>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <list>
#include "Curve.h"

using namespace std;

namespace CubicS {

struct piece;

}

class CubicSpline : public Curve
{
  public:
      CubicSpline();

      //Optimization
      void S_splineFeasible();
      void S_splineFeasibleSmart();
      void S_splineFeasibleSmart3(int index);
      void optimizeCurve(double** nodes, int numNodes, InitDir initDir);
      void fineOptimizeCurve(double** nodes, int numNodes, int depth);
      int S_Curve(double *a, double *b, double *rx, double *ry);
      Curve::Node pointInTheMiddle(double* a, double* b);


      void boundingRect(double &x, double &y, double &width, double &height);

      bool *flag_feasible, *flag_optimal, flag_closed, *pieceValid, *respectCorner;

};


#endif // CUBICSPLINE_H
