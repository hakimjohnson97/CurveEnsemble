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


#include "Curve.h"
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <time.h>
#include <QDebug>
#include "Curve_Libs/common_curve_lib.h"

using namespace std;

#define M_PI_180        0.0174532925199432957692369076848861271
#define MAX_NODES 1000


double bufx[600];
double bufy[600];

double DistBetweenPoints2(double x1, double y1, double x2, double y2)
{
       return sqrt(pow(x2-x1, 2) + pow(y2-y1, 2));
}


Curve::Node::Node()
{
    mX = 0;
    mY = 0;
    mDir = 0;
    mCorner = 0;
    mState = REGULAR;
    mRespectCorner = 0;
}

Curve::Node::Node(double x, double y, bool regular)
{
    mX = x;
    mY = y;
    mDir = 0;
    mCorner = 0;
    if (regular == true)
        mState = REGULAR;
    else
        mState = FREE;
    mRespectCorner = 0;
}

Curve::Node::Node(double* node, bool arg_respectCorner)
{
    mX = node[X];
    mY = node[Y];
    mDir = node[DIR];
    mCorner = node[CORNER];
    mState = node[STATE];
    mRespectCorner = arg_respectCorner;
}

void Curve::Node::setPos(double x, double y)
{
     mX = x;
     mY = y;
}

double Curve::Node::x()
{
    return mX;
}

double Curve::Node::y()
{
    return mY;
}


void Curve::Node::setRegular(bool regular)
{
    if (regular == true)
        mState = REGULAR;
    else
        mState = FREE;
}

bool Curve::Node::isRegular()
{
     if (mState == REGULAR)
         return true;
     return false;
}

double Curve::Node::dataBallRadius()
{
    if (mState < 0)
        return -mState;
    else
        return mState-100;
}

void Curve::Node::toggleTick()
{
     mState = 100-mState;
}

bool Curve::Node::isTicked()
{
    if (mState < 0)
        return false;
    else
        return true;
}

//Regular Nodes

void Curve::Node::setCorner(double corner)
{
    mCorner = mod360(corner);
}

double Curve::Node::corner()
{
    return mCorner;
}

double Curve::Node::incDir()
{
    return mod360(mDir - mCorner);
}

double Curve::Node::outDir()
{
        return mod360(mDir);
}

//Irregular Nodes

int Curve::Node::setIncClamp(double incDir)
{
    if (mState == REGULAR)
        return false;
    if (mRespectCorner)
        mDir = mod360(incDir + mCorner);
    else
       mCorner = mod360(mDir - incDir);
       setIncState(INC_CLAMPED);

    return true;
}

int Curve::Node::setOutClamp(double outDir)
{
    if (mState == REGULAR)
        return false;
    if (!mRespectCorner)
        mCorner = mod360(mCorner + (outDir - mDir));
    mDir = mod360(outDir);
    setOutState(OUT_CLAMPED);
    return true;
}

int Curve::Node::removeIncClamp()
{
    setIncState(INC_FREE);
    return true;
}

int Curve::Node::removeOutClamp()
{
    setOutState(OUT_FREE);
    return true;
}

int Curve::Node::respectCorner(bool respect)
{
    mRespectCorner = respect;
    return true;
}

bool Curve::Node::hasIncClamp()
{
    if (mState == CLAMPED_FREE || mState == CLAMPED)
        return true;
    return false;
}

bool Curve::Node::hasOutClamp()
{
    if (mState == FREE_CLAMPED || mState == CLAMPED)
        return true;
    return false;
}

bool Curve::Node::isCornerRespected()
{
    return mRespectCorner;
}

void Curve::Node::setIncState(int state)
{
          if (mState == REGULAR)    {
              if (state == INC_CLAMPED)
                  mState = CLAMPED;
              else if (state == INC_FREE)
                  mState = FREE;
          }
          else if (mState == CLAMPED)    {
              if (state == INC_CLAMPED)
                  mState = CLAMPED;
              else if (state == INC_FREE)
                  mState = FREE_CLAMPED;
          }
          else if (mState == FREE)    {
              if (state == INC_CLAMPED)
                  mState = CLAMPED_FREE;
              else if (state == INC_FREE)
                  mState = FREE;
          }
          else if (mState == FREE_CLAMPED)    {
              if (state == INC_CLAMPED)
                  mState = CLAMPED;
              else if (state == INC_FREE)
                  mState = FREE_CLAMPED;
          }
          else if (mState == CLAMPED_FREE)    {
              if (state == INC_CLAMPED)
                  mState = CLAMPED_FREE;
              else if (state == INC_FREE)
                  mState = FREE;
          }
}

void Curve::Node::setOutState(int state)
{
          if (mState == REGULAR)    {
              if (state == OUT_CLAMPED)
                  mState = CLAMPED;
              else if (state == OUT_FREE)
                  mState = FREE;
          }
          else if (mState == CLAMPED)    {
              if (state == OUT_CLAMPED)
                  mState = CLAMPED;
              else if (state == OUT_FREE)
                  mState = CLAMPED_FREE;
          }
          else if (mState == FREE)    {
              if (state == OUT_CLAMPED)
                  mState = FREE_CLAMPED;
              else if (state == OUT_FREE)
                  mState = FREE;
          }
          else if (mState == CLAMPED_FREE)    {
              if (state == OUT_CLAMPED)
                  mState = CLAMPED;
              else if (state == OUT_FREE)
                  mState = CLAMPED_FREE;
          }
          else if (mState == FREE_CLAMPED)    {
              if (state == OUT_CLAMPED)
                  mState = FREE_CLAMPED;
              else if (state == OUT_FREE)
                  mState = FREE;
          }
}

void Curve::Node::setDataBall(bool x, double radius)
{
     if (x == 1)
         mState = -radius;
     else if (x == 0)    {
         mState = REGULAR;
         mCorner = 0;
     }
         
}

bool Curve::Node::isDataBall()
{
    if (mState < 0 || mState > 100)
        return true;
    else
        return false;
}






Curve::Curve(string method)
{
    curveMethod = method;


    int i;
    H = new double*[MAX_NODES];
    for (i = 0; i < MAX_NODES; i++)
        H[i] = new double [6];
    Htemp = new double*[MAX_NODES];
    for (i = 0; i < MAX_NODES; i++)
        Htemp[i] = new double [6];
    Htemp2 = new double*[MAX_NODES];
    for (i = 0; i < MAX_NODES; i++)
        Htemp2[i] = new double [6];
    phant = new double*[100];
    for (i = 0; i < 100; i++)
        phant[i] = new double [6];

    for (i = 0; i < 100; i++)
        H[i][CORNER] = 0;

    for (i = 0; i < 100; i++)
        H[i][STATE] = REGULAR;



    respectCorner = new bool[MAX_NODES];
    for (i = 0; i < MAX_NODES; i++)
        respectCorner[i] = false;

    count = 0;
    phantCount = 0;
    closed = false;



    subCount = 0;

}

void Curve::addNode(Node node)
{
     H[count][X] = node.mX;
     H[count][Y] = node.mY;
     H[count][DIR] = node.mDir;
     H[count][STATE] = node.mState;
     H[count][CORNER] = node.mCorner;
     respectCorner[count] = 0;
     count++;

     if (closed == false)    {
         if (count >= 3)    {
             setState(count-1, OUT_FREE);

            H[count-2][STATE] = REGULAR;
         }
         else     {

             setState(count-1, INC_FREE);
         }
     }

     pointChanged(count-1);
}

void Curve::insertNode(int index, Node node)
{
     if (index > count-1 || index < -1)
         return;
     if (index == count-1)    {
         addNode(node);
         return;
     }
     int i;
     for (i = count; i > index+1; i--)    {
         H[i][X] = H[i-1][X];
         H[i][Y] = H[i-1][Y];
         H[i][DIR] = H[i-1][DIR];
         H[i][STATE] = H[i-1][STATE];
         H[i][CORNER] = H[i-1][CORNER];
     }
     H[index+1][X] = node.mX;
     H[index+1][Y] = node.mY;
     H[index+1][DIR] = node.mDir;
     H[index+1][STATE] = node.mState;
     H[index+1][CORNER] = node.mCorner;
     count++;
     pointChanged(index+1);
     S_splineFeasibleSmart();
}

void Curve::deleteNode(int index)
{
    if (count <= 0)
        return;

    int i;
    pointChanged(index);
    for (i = index; i < count-1; i++)    {
        H[i][X] = H[i+1][X];
        H[i][Y] = H[i+1][Y];
        H[i][DIR] = H[i+1][DIR];
        H[i][STATE] = H[i+1][STATE];
        H[i][CORNER] = H[i+1][CORNER];
        H[i][BE] = H[i+1][BE];
    }
    if (closed == false)    {
        if (index == 0)
            H[0][STATE] = FREE;
        if (index == count-1)
            if (count-2 >= 0)
                H[count-2][STATE] = FREE;
    }

   count--;

}

void Curve::moveCurve(double arg_x, double arg_y)
{
    int i;
    for (i = 0; i < count; i++)    {
        H[i][X] += arg_x;
        H[i][Y] += arg_y;
    }
}

void Curve::movePoint(int index, double arg_x, double arg_y)
{
     H[index][X] = arg_x;
     H[index][Y] = arg_y;
     pointChanged(index);
}

void Curve::setCorner(int index, double corner)
{
     H[index][CORNER] = corner;
     pointChanged(index);
}

void Curve::setDir(int index, double dir)
{
     H[index][DIR] = dir;
     pointChanged(index);
}

void Curve::setState(int index, int arg_state)
{
     if (arg_state >= REGULAR && arg_state <= FREE)
         H[index][STATE] = arg_state;
     else if (arg_state == INC_CLAMPED || arg_state == INC_FREE)    {
          if (H[index][STATE] == REGULAR)    {
              if (arg_state == INC_CLAMPED)
                  H[index][STATE] = CLAMPED;
              else if (arg_state == INC_FREE)
                  H[index][STATE] = FREE;
          }
          else if (H[index][STATE] == CLAMPED)    {
              if (arg_state == INC_CLAMPED)
                  H[index][STATE] = CLAMPED;
              else if (arg_state == INC_FREE)
                  H[index][STATE] = FREE_CLAMPED;
          }
          else if (H[index][STATE] == FREE)    {
              if (arg_state == INC_CLAMPED)
                  H[index][STATE] = CLAMPED_FREE;
              else if (arg_state == INC_FREE)
                  H[index][STATE] = FREE;
          }
     }
     else if (arg_state == OUT_CLAMPED || arg_state == OUT_FREE)    {
          if (H[index][STATE] == REGULAR)    {
              if (arg_state == OUT_CLAMPED)
                  H[index][STATE] = CLAMPED;
              else if (arg_state == OUT_FREE)
                  H[index][STATE] = FREE;
          }
          else if (H[index][STATE] == CLAMPED)    {
              if (arg_state == OUT_CLAMPED)
                  H[index][STATE] = CLAMPED;
              else if (arg_state == OUT_FREE)
                  H[index][STATE] = CLAMPED_FREE;
          }
          else if (H[index][STATE] == FREE)    {
              if (arg_state == OUT_CLAMPED)
                  H[index][STATE] = FREE_CLAMPED;
              else if (arg_state == OUT_FREE)
                  H[index][STATE] = FREE;
          }
     }

}

int Curve::getOutState(int index)
{
    if (H[index][STATE] == REGULAR)
        return REGULAR;
    if (H[index][STATE] == CLAMPED || H[index][STATE] == FREE_CLAMPED)
        return OUT_CLAMPED;
    if (H[index][STATE] == FREE || CLAMPED_FREE)
        return OUT_FREE;
}

void Curve::setIncState(int index, int arg_state)
{

}

void Curve::setOutState(int index, int arg_state)
{

}

void Curve::setClosed(bool arg_closed)
{
     if (arg_closed == false)    {

         Node node;
         node = getNode(0);
         node.setRegular(false);
         node.removeIncClamp();
         setNode(0, node);
         node = getNode(count-1);
         node.setRegular(false);
         node.removeOutClamp();
         setNode(count-1, node);
     }
     else    {

         H[0][STATE] = REGULAR;
         H[count-1][STATE] = REGULAR;
     }
     closed = arg_closed;
     flag_closed8 = arg_closed;
     // S_splineFeasibleSmart(); I removed it cuz of optimize on fly was implemented
}


//Transformations

complex Complex(double real, double imag)
{
    complex c;
    c.real = real;
    c.imag = imag;
    return c;
}

void Curve::transform(double c1real, double c1imag, double c2real, double c2imag, bool reflection)
{
     int i;
     complex c1, c2;
     c1.real = c1real;
     c1.imag  = c1imag;
     c2.real = c2real;
     c2.imag = c2imag;

     complex z;
     complex polar;
     polar = complex_polar(c1);
     for (i = 0; i < count; i++)    {
         z.real = H[i][X];
         z.imag = H[i][Y];
         if (reflection)
             z = complex_conj(z);
         z = complex_add(complex_mult(z, c1), c2);


         H[i][X] = z.real;
         H[i][Y] = z.imag;
         H[i][BE] = H[i][BE]/polar.real;

         if (H[i][STATE] >= 0 && H[i][STATE] < 100) //Data Ball
             H[i][DIR] = mod360(H[i][DIR] + polar.imag);
         else    {
             H[i][DIR] = mod360(H[i][DIR] + polar.imag);
             H[i][CORNER] = mod360(H[i][CORNER] + polar.imag);
             double r;
             if (H[i][STATE] < 0)    {
                 r = complex_abs(c1)*(-H[i][STATE]);
                 H[i][STATE] = -r;
             }
             else    {
                 r = complex_abs(c1)*(H[i][STATE] - 100);
                 H[i][STATE] = 100 + r;
             }
             
         }

     }


}

void Curve::transformWithPivot(int index, double newPosX, double newPosY, double pivotX, double pivotY)
{
     complex oldPos = Complex(H[index][X], H[index][Y]);
     complex newPos = Complex(newPosX, newPosY);
     complex pivot = Complex(pivotX, pivotY);

     complex c1, c2;
     c1 = complex_div(complex_sub(newPos, pivot), complex_sub(oldPos, pivot));
     c2 = complex_sub(pivot, complex_mult(c1, pivot));
     transform(c1.real, c1.imag, c2.real, c2.imag, false);

}

void Curve::rotate(int index, double newPosX, double newPosY, double pivotX, double pivotY)
{
     complex oldPos = Complex(H[index][X], H[index][Y]);
     complex newPos = Complex(newPosX, newPosY);
     complex pivot = Complex(pivotX, pivotY);

     complex c1, c2;
     c1 = complex_div(complex_sub(newPos, pivot), complex_sub(oldPos, pivot));
     double mag = complex_abs(c1);
     c1.real = c1.real / mag;
     c1.imag = c1.imag / mag;

     c2 = complex_sub(pivot, complex_mult(c1, pivot));
     transform(c1.real, c1.imag, c2.real, c2.imag, false);

}

void Curve::rotate(double angle, double pivotX, double pivotY)
{
     double radians = angle*M_PI_180;
     complex pivot = Complex(pivotX, pivotY);
     complex c1, c2;
     c1 = Complex(cos(radians), sin(radians));
     c2 = complex_sub(pivot, complex_mult(c1, pivot));
     transform(c1.real, c1.imag, c2.real, c2.imag, false);

}

void Curve::resize(int index, double newPosX, double newPosY, double pivotX, double pivotY)
{
     complex oldPos = Complex(H[index][X], H[index][Y]);
     complex newPos = Complex(newPosX, newPosY);
     complex pivot = Complex(pivotX, pivotY);

     complex c1, c2;
     c1 = complex_div(complex_sub(newPos, pivot), complex_sub(oldPos, pivot));
     c1.real = complex_abs(c1);
     c1.imag = 0;
     c2 = complex_sub(pivot, complex_mult(c1, pivot));
     transform(c1.real, c1.imag, c2.real, c2.imag, false);

}

void Curve::resize(double scaleFactor, double pivotX, double pivotY)
{
     complex pivot = Complex(pivotX, pivotY);
     complex c1, c2;
     c1.real = scaleFactor;
     c1.imag = 0;
     c2 = complex_sub(pivot, complex_mult(c1, pivot));
     transform(c1.real, c1.imag, c2.real, c2.imag, false);
}

void Curve::reflect(double pos1X, double pos1Y, double pos2X, double pos2Y)
{
     complex pos1 = Complex(pos1X, pos1Y);
     complex pos2 = Complex(pos2X, pos2Y);
     complex c1, c2;
     c1 = complex_div(complex_sub(pos1, pos2), complex_sub(complex_conj(pos1), complex_conj(pos2)));
     c2 = complex_sub(pos2, complex_mult(c1, complex_conj(pos2)));
     transform(c1.real, c1.imag, c2.real, c2.imag, true);
}


int Curve::save(fstream &file, bool old)
{
     char emptySpace[35];
     int i;
     for (i = 0; i < 35; i++)
         emptySpace[i] = 0;
     file.write((char*)&count, sizeof(long int));
     for (i = 0; i < count; i++)    {
         file.write((char*)H[i], sizeof(double)*6);
     }
     file.write((char*)&closed, sizeof(bool));

     file.write(emptySpace, sizeof(int)); //coarse res is now const
     for (i = 0; i < count; i++)
         file.write((char*)&(respectCorner[i]), sizeof(bool));
     if (old)    {
         if ((35 - count) >= 0)
             file.write(emptySpace, 35 - count);
     }
     else
         file.write(emptySpace, 35);

     return true;
}

int Curve::load(fstream &file, bool old)
{
     char emptySpace[35];
     int i;
     file.read((char*)&count, sizeof(long int));
     for (i = 0; i < count; i++)    {
         file.read((char*)H[i], sizeof(double)*6);
     }
     file.read((char*)&closed, sizeof(bool));

     file.read(emptySpace, sizeof(int));
     for (i = 0; i < count; i++)
         file.read((char*)&(respectCorner[i]), sizeof(bool));

     if (old)    {
         if ((35 - count) >= 0)
             file.read(emptySpace, 35 - count);
     }
     else
         file.read(emptySpace, 35);


     return true;
}

int Curve::saveAsText(fstream &file)
{
    int i, j;
    file.precision(20);
    for (i = 0; i < count; i++)    {
        for (j = 0; j < 6; j++)    {
            file << H[i][j] << " ";
        }
        file << "\n";
    }
}

int Curve::loadFromText(fstream &file)
{

    int i, j;
    bool lineRead, three_col = true;
    string s;



    i = 0;

    while (!file.eof())    {
        getline(file, s);
        stringstream line (s);
        if (line.peek() == '#')
            continue;
        j = 0;
        lineRead = false;
        while (j < 6 && (!line.eof()))    {
            line >> H[i][j];
            lineRead = true;
            if (j > 2)
                three_col = false;
            j++;
        }
        if (j < 6)    {
            for (int k = j; k < 6; k++)
                H[i][k] = 0;
        }
        if (lineRead)
            i++;
    }
    count = i;

    if ((H[0][4] == FREE || H[0][4] == FREE_CLAMPED) && (H[count-1][4] == FREE || H[count-1][4] == CLAMPED_FREE))     {
        closed = false;
    }
    else
        closed = true;

    if (curveMethod == "DataBalls" && three_col){
        for (i=0; i<count; i++)
            {H[i][STATE]=-H[i][DIR]; H[i][DIR]=0;}
    }
    return 1;
}




int Curve::pointChanged(int index)
{

}

Curve::Node Curve::getNode(int index)
{
     if (index < count && index > -1)    {
         Node node(H[index], respectCorner[index]);
         return node;
     }
     else if (index == count)    {
         Node node(H[0], respectCorner[0]);
         return node;
     }
     else
         return Node();

}

void Curve::setNode(int index, Node node)
{
    if (index >= count || index <= -1)
         return;
    if (H[index][X] == node.mX  && H[index][Y] == node.mY && H[index][DIR] == node.mDir
          &&  H[index][CORNER] == node.mCorner && H[index][STATE] == node.mState &&  respectCorner[index] == node.mRespectCorner)

        return;


     H[index][X] = node.mX; 
     H[index][Y] = node.mY;
     H[index][DIR] = node.mDir;
     H[index][CORNER] = node.mCorner;
     respectCorner[index] = node.mRespectCorner;

      if (closed == false)    {
          if (node.mState != REGULAR)    {
              if (index == 0 && node.hasIncClamp() == false)
                  H[index][STATE] = node.mState;
              if (index == count-1 && node.hasOutClamp() == false)
                  H[index][STATE] = node.mState;
          }
          if (index != count-1 && index != 0)
              H[index][STATE] = node.mState;

      }
      else
          H[index][STATE] = node.mState;


}

void Curve::setAlphaMax(double alphaMax)
{
   alpha_max8 = alphaMax*M_PI_180;
}

double Curve::alphaMax()
{
   return alpha_max8/M_PI_180;
}


void Curve::setLambda(double arg_lambda)
{
   lam8 = arg_lambda;
}

double Curve::lambda()
{
   return lam8;
}


void Curve::copy(Curve* curve)
{
    count = curve->count;

    closed = curve->closed;

    phantCount = curve->phantCount;
    subCount = curve->subCount;


        for (int i = 0; i < count; i++)
            for (int j = 0; j < 6; j++)
                H[i][j] = curve->H[i][j];
        for (int i = 0; i < phantCount; i++)
            for (int j = 0; j < 6; j++)
                phant[i][j] = curve->phant[i][j];
        for (int i = 0; i < count; i++)
            respectCorner[i] = curve->respectCorner[i];
}

int Curve::pointsOnCurve(int index, double *rx, double *ry, int subdividePieceNum)
{
    if (count < 2)
        return 0;

    int indexp1 = (index+1) % count;
    int subNum = getPhantNodesPerPiece();

    if (getSubdivideCount() == 0)
        return S_Curve(H[index], H[indexp1], rx, ry);
    else    {
        if (subdividePieceNum == 1)
            return S_Curve(H[index], phant[index*subNum], rx, ry);
        else if (subdividePieceNum < (subNum+1))
            return S_Curve(phant[index*subNum+subdividePieceNum-2], phant[index*subNum+subdividePieceNum-1], rx, ry);
        else if (subdividePieceNum == (subNum+1))
            return S_Curve(phant[index*subNum+subNum-1], H[indexp1], rx, ry);
    }
    return 0;
}

int Curve::S_Curve(double *a, double *b, double *rx, double *ry)
{ //Default treats the curve as Lines
    rx[0] = a[0];
    rx[1] = b[0];
    ry[0] = a[1];
    ry[1] = b[1];
    return 2;
}

void Curve::subdivide()
{
    int subNum = getPhantNodesPerPiece();
    int k = 0;
    int r = 0;


    for (int i = 0; i < getPieceCount(); i++)    {
        Node a;
        int ip1 = (i+1) % count;
        qDebug() << "Ip1: " << ip1;
        if (subCount == 0)    {
            a = pointInTheMiddle(H[i], H[ip1]);
            Htemp[k][X] = a.mX;
            Htemp[k][Y] = a.mY;
            k++;
        }
        else    {
            a = pointInTheMiddle(H[i], phant[r]);    
                Htemp[k][X] = a.mX;
                Htemp[k][Y] = a.mY;
                k++;

                Htemp[k][X] = phant[r][X];
                Htemp[k][Y] = phant[r][Y];
                k++; r++;

                for (int j = 0; j < subNum-1; j++)    {
                    a = pointInTheMiddle(phant[r-1], phant[r]);
                    Htemp[k][X] = a.mX;
                    Htemp[k][Y] = a.mY;
                    k++;

                    Htemp[k][X] = phant[r][X];
                    Htemp[k][Y] = phant[r][Y];
                    k++; r++;
                }

                a = pointInTheMiddle(phant[r-1], H[ip1]);
                Htemp[k][X] = a.mX;
                Htemp[k][Y] = a.mY;
                k++;
       }

    }


    phantCount = 0;

    for (int i = 0; i < k; i++) {
        phant[i][X] = Htemp[i][X];
        phant[i][Y] = Htemp[i][Y];
        phant[i][STATE] = 0;
        phant[i][CORNER] = 0;
        phant[i][DIR] = 0;
        phant[i][BE] = 0;
        phantCount++;
    }

    subCount++;

    optimize(Smart);
    fineOptimize();

    qDebug() << "PhantCount: " << phantCount;
}

void Curve::endSubdivide()
{
    phantCount = 0;
    subCount = 0;

    optimize(Smart);
}

void Curve::anchorPhantomCurve()
{
    int n = createPhantomCurve();
    for (int i = 0; i < n; i++)
         for (int j = 0; j < 6; j++)
             H[i][j] = Htemp[i][j];
    count = n;
    endSubdivide();
}

QList<Curve::Node> Curve::getPhantNodes(int index)
{
    int subNum =  getPhantNodesPerPiece();
    QList<Node> x;
    for (int i = 0; i < subNum; i++)
        x.append(phant[subNum*index + i]);
    return x;
}

void copyNode(double* to, double* from)
{
    to[X] = from[X];
    to[Y] = from[Y];
    to[DIR] = from[DIR];
    to[STATE] = from[STATE];
    to[CORNER] = from[CORNER];
    to[BE] = from[BE];
}

int Curve::createPhantomCurve()
{
        int k = 0;
        int subNum = getPhantNodesPerPiece();
        for (int i = 0; i < getPieceCount(); i++)    {
            copyNode(Htemp[k], H[i]);
            k++;
            for (int j = 0; j < subNum; j++)    {
                copyNode(Htemp[k], phant[i*subNum + j]);
                k++;
            }
        }
        if (!closed)    {
            copyNode(Htemp[k], H[count-1]);
            k++;
        }
        return k;

}

void Curve::destroyPhantomCurve()
{
    int k = 0;
    int subNum = getPhantNodesPerPiece();
    for (int i = 0; i < getPieceCount(); i++)    {
        copyNode(H[i], Htemp[k]);
        k++;
        for (int j = 0; j < subNum; j++)    {
            copyNode(phant[i*subNum + j], Htemp[k]);
            k++;
        }
    }
    if (!closed)    {
        copyNode(H[count-1], Htemp[k]);
        k++;
    }
}

void Curve::optimize(InitDir initDir, bool onFly)
{
    if (count <= 0)
        return;

    if (subCount == 0)    {
        if (onFly == false)
            optimizeCurve(H, count, initDir);
        else    {
            flag_quick8 = true;
            alpha_maxQ8=alpha_max8; lamQ8=lam8; tolQ8=tol;

            optimizeCurve(H, count, initDir);
            flag_quick8 = false;
        }
    }
    else    {
        int n = createPhantomCurve();
        optimizeCurve(Htemp, n, initDir);
        destroyPhantomCurve();
    }
}

void Curve::copyToTemp()
{
    for (int i = 0; i < count; i++)
        for (int j = 0; j < 6; j++)
            Htemp2[i][j] = H[i][j];
    count_temp = count;

    alpha_maxQ8=alpha_max8; lamQ8=lam8; tolQ8=tol;
}

void Curve::optimizeToTemp(InitDir initDir)
{
    flag_quick8 = true;

    optimizeCurve(Htemp2, count_temp, initDir);
    flag_quick8 = false;
}

void Curve::copyFromTemp()
{
    for (int i = 0; i < count; i++)
        for (int j = 0; j < 6; j++)
            H[i][j] = Htemp2[i][j];
}

void Curve::fineOptimize(bool onFly)
{
    alpha_maxQ8=alpha_max8; lamQ8=lam8; tolQ8=tol;
    flag_quick8 = onFly;
    if (subCount == 0)
        fineOptimizeCurve(H, count);
    else    {
        int n = createPhantomCurve();
        fineOptimizeCurve(Htemp, n);
        destroyPhantomCurve();
    }
    flag_quick8 = false;
}

void Curve::convertCirclesToDataBalls()
{
    double **B = Htemp;
    int n = count;

    int i,j;
    double x0,y0,x1,y1;
    i=0; j=0;
    while (i+1<n){
      x0=H[i][X]; y0=H[i][Y]; x1=H[i+1][X]; y1=H[i+1][Y];
      B[j][X]=0.5*(x0+x1); B[j][Y]=0.5*(y0+y1);
      B[j][STATE]=-0.5*sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
      B[j][NODEDIR]=0; B[j][DIR]=0; B[j][BE]=0;
      i+=2; j++;
    }// while
    count = j;


    for (int i = 0; i < count; i++)
        for (int j = 0; j < 6; j++)
            H[i][j] = Htemp[i][j];

    closed = true;
}

double Curve::createDataBallCurve(Curve* dataBallCurve)
{
     count = 0;

    copy(dataBallCurve);
    closed = false;

    H[0][STATE] = FREE;
    H[count-1][STATE] = FREE;
    for (int i = 1; i < count-1; i++)     {//ASSUMING ITS NOT CLOSED
        if (H[i][STATE] < 0)
            H[i][STATE] = 100 - H[i][STATE];
    }


    return runCommand(0);

 
}

double Curve::getCurveBE()
{
    double total = 0;
    for (int i = 0; i < count; i++)
        total += H[i][BE];
    return total;
}

double Curve::getCurveLength()
{
    double total = 0;
    int lastIndex;
    if (closed == true)
        lastIndex = count;
    else
        lastIndex = count-1;

    for (int i = 0; i < lastIndex; i++)    {
        int k = pointsOnCurve(i, bufx, bufy);
        for (int j = 1; j < k; j++)
            total += sqrt((bufx[j]-bufx[j-1])*(bufx[j]-bufx[j-1]) + (bufy[j]-bufy[j-1])*(bufy[j]-bufy[j-1]));
    }
    return total;
}

void Curve::clearExceptPoints()
{
    for (int i = 0; i < count; i++)    {
        for (int j = 2; j < 6; j++)    {
            if (j == STATE)
                continue;
            H[i][j] = 0;
        }
    }
    for (int i = 0; i < count; i++)     {
        if (i == 0 || i == count-1)   {
             if (closed == false)
                 H[i][STATE] = FREE;
             else
                 H[i][STATE] = REGULAR;
        }
        else
            H[i][STATE] = REGULAR;
    }
}

int Curve::furthest_marked(){
// We assume the first and last nodes of H are unmarked (either normal nodes or unmarked
// data balls).  Let S be the curve through the unmarked nodes of H.  In this procedure,
// we find the marked data ball which lies furthest from S.  Note that the distance in
// question refers to the distance from the boundary of the data ball to S.  
// BTW, the above specified distance is denoted min_dist in the code.
// Presently, initial values are min_dist = 1e-3 and j_opt=-1.  
// The returned integer is j_opt, whereby H[j_opt] is the marked data ball furthest from S.
// For your convenience, we unmark H[j_opt] and we set its NODEDIR pointing the the nearest
// point in S.
// Exception: If all data balls are within 1e-3 distance from S, j_opt=-1 is returned and
// H is unchanged (this indicates that smoothest path operation has terminated). 
int i,ip1,j,j_opt,k=0;
double md,r,min_dist,x[600],y[600],w[2],w_opt[2]; //min_dist generically denotes the
ip1=0; j_opt=-1; min_dist=max_radius(H,count)/25000.0;
while (ip1<count-1){
  i=ip1;
  ip1=i+1;
  while (H[ip1][STATE]>100) ip1++;
  // At this point the basic curve runs from node i to node ip1
  k=S_Curve (H[i],H[ip1],x,y);
  //printf("num points in S_curve = %d \n",k);
  for (j=i+1; j<ip1; j++){
    r=H[j][STATE]-100; md=closest_point_thresh(w,H[j],x,y,k,min_dist+r); 
//    printf("j=%d, min_dist=%e \n",j,md);
    if (md-r>min_dist){
      min_dist=md-r; j_opt=j; w_opt[0]=w[0]; w_opt[1]=w[1];
    }// if
  }//for j
}//while
if (j_opt >= 0)
  H[j_opt][NODEDIR]=M_180_PI*direction_angle(w_opt[0]-H[j_opt][X],w_opt[1]-H[j_opt][Y]);
printf("j_opt = %d and w= (%e, %e) \n",j_opt,w_opt[0],w_opt[1]);
return j_opt;
}// furthest_marked-----------------------------------------------------------------------

double Curve::runCommand(int c1, double c2)
{

    int j;

    double R;
    if (c1==1) // mark all databalls in selected curve
       mark_all(H,count);
    else if (c1==2) // unmark all databalls in selected curve
       unmark_all(H,count);
    else if (c1==3) // extract unmarked data balls from selected curve
       count=extract_unmarked(H,count);
    else if (c1==4){ // unmark the furthest marked data ball of selected curve
       j=furthest_marked(); 
       if (j>=0) H[j][STATE]=100-H[j][STATE];
    }//else if
    else if (c1==5) // anchor all data balls
       anchor_data_balls(H,count);
    else if (c1==6) // set DIR to mean stencil directions
       S_splineFeasibleSmart();
    else if (c1==7) // report maximum data ball radius
       c2=max_radius(H,count);
    else if (c1==8){// scale curve so max data ball radius = c2
      R=max_radius(H,count);
      if (R>0){
        // void transform(double c1x, double c1y, double c2x, double c2y, bool reflection);
        transform(c2/R,0,0,0,false);
        // now maximum radius equals c2 (input value) and current c2 equals multiplication factor
      } // if
      transform(1,0,-H[0][0],-H[0][1],false);  // now first node is (0,0)
    }// else if
    else if (c1==9){ // scale radii of all data balls by c2 (if c2>0)
      for (int i=0; i<count; i++){
        if (c2>0 && H[i][STATE]<0)
           H[i][STATE] *= c2;
        else if (c2>0 && H[i][STATE]>100)
           H[i][STATE]=100+c2*(H[i][STATE]-100);
      }// for
    }// else if
    return c2;
}// runCommand





//Not used anmymore
/*void Curve::optimizeDataBalls()
{
    double theta[coarseRes];
    equi_angular_partition(theta, coarseRes); printf("Curve::optimizeDataBalls: hello \n");
    DataBalls::Optimize_lines(H, count, theta, coarseRes, true);

     int i;
     double h = 360/double(coarseRes);
     int k = 8;
     double theta2[2*k];
     for (i = 0; i < 10; i++)    {
         h /= 4;
         equi_angular_cluster(theta2, h, k);
         DataBalls::Optimize_lines(H, count, theta2, 2*k, false);
     }
}

void Curve::createShortestPath(Curve* curve)
{
    curve->closed = true;
    for (int i = 0; i < count; i++)
        curve->addNode(Node());

    DataBalls::data_balls_2_lines(H, count, curve->H, 1); 

    curve->closed = true;
    for (int i = 0; i < count; i++)
        curve->addNode(Node());

    DataBalls::data_balls_2_lines(H, count, curve->H, 1); 
}


void Curve::createSmoothestCurve(Curve* curve)
{
    curve->closed = true;
    for (int i = 0; i < count; i++)
        curve->addNode(Node(H[i][X], H[i][Y]));
}

void Curve::createSmoothestCurveStep(Curve* curve, double h)
{

} */







