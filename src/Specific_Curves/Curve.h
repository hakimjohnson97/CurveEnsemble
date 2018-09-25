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


#ifndef CURVE_H
#define CURVE_H

#include <fstream>

#include <QList>

using namespace std;


#define INC_CLAMPED 8
#define INC_FREE    9

#define OUT_CLAMPED 10
#define OUT_FREE    11

const int X = 0;
const int Y = 1;
const int DIR = 2;
const int CORNER = 3;
const int STATE = 4;
const int BE = 5;

const int FREE = 7;
const int REGULAR = 0;
const int CLAMPED = 1;
const int FREE_CLAMPED = 3;
const int CLAMPED_FREE = 5;

//struct complex {double real,imag;};


/*
  Curve is an abstract class representing a curve. Curve by itself represents just the nodes with
  no curve method associated with them (i.e no way of drawing the curve) For a curve method to be used
  in this project a new class must be created that inherits Curve and specifies a way to draw and optimize
  the curve by overwriting the virtual methods.
  The base class provides the interface for GraphArea to communicate with these curves.
*/
class Curve
{
public:
    enum InitDir {Smart, Random, AsIs};

    class Node
    {
      public:
          Node();
          Node(double x, double y, bool regular = true);
          Node(double* node, bool arg_respectCorner = 0);
          void setPos(double x, double y);
          void setCorner(double corner);
          void setRegular(bool regular);

          double x();
          double y();
          double outDir();
          double incDir();
          double corner();
          bool isRegular();


          //Irregular Nodes
          int setIncClamp(double dir);
          int setOutClamp(double dir);
          int removeIncClamp();
          int removeOutClamp();
          int respectCorner(bool respect);

          bool hasIncClamp();
          bool hasOutClamp();
          bool isCornerRespected();

     // private:
          void setIncState(int state);
          void setOutState(int state);

          void setDataBall(bool x, double radius = 25);
          bool isDataBall();
          double dataBallRadius();
          void toggleTick();
          bool isTicked();

          double mX, mY;
          double mDir;
          double mCorner;
          double mState;
          bool mRespectCorner;

      friend class Curve;
    };


    Curve(string method = "");

          double **H, **phant, **Htemp, **Htemp2;

    void addNode(Node node);
    void insertNode(int index, Node node);
    void deleteNode(int index);
    //Setters
    void movePoint(int index, double x, double y);
    void moveCurve(double x, double y);
    void setCorner(int index, double corner);
    void setDir(int index, double dir);
    void setState(int index, int state);
    void setIncState(int index, int state);
    void setOutState(int index, int state);
    void setClosed(bool closed);
    void setAlphaMax(double alphaMax);
    void setLambda(double arg_lambda);

    int getOutState(int index);

    //Transformations
    void transform(double c1x, double c1y, double c2x, double c2y, bool reflection);
    void transformWithPivot(int index, double newPosX, double newPosY, double pivotX, double pivotY);
    void rotate(int index, double newPosX, double newPosY, double pivotX, double pivotY);
    void rotate(double angle, double pivotX, double pivotY);
    void resize(int index, double newPosX, double newPosY, double pivotX, double pivotY);
    void resize(double scaleFactor, double pivotX, double pivotY);
    void reflect(double pos1X, double pos1Y, double pos2X, double pos2Y);

    //Getters
    int getCount() {return count;};
    int getPieceCount() {if (closed) return count; else return count-1;}
    bool isClosed() {return closed;};
    double alphaMax();
    double lambda();
    Node getNode(int index);
    QList<Node> getPhantNodes(int index);
    void anchorPhantomCurve();
    int getSubdivideCount() {return subCount;}
    int getPhantNodesPerPiece() {return phantCount/getPieceCount();}
    void setNode(int index, Node node);


    int save(fstream &file, bool old = 0);
    int load(fstream &file, bool old = 0);
    int saveAsText(fstream &file);
    int loadFromText(fstream &file);

    int pointsOnCurve(int index, double *rx, double *ry, int subdividePieceNum = 1);
    double getCurveLength();
    void optimize(InitDir initDir = Smart, bool onFly = false);

    void copyToTemp();
    void optimizeToTemp(InitDir initDir = Smart);
    void copyFromTemp();

    int createPhantomCurve();
    void destroyPhantomCurve();
    void fineOptimize(bool onFly = false);

    virtual void S_splineFeasible() {}
    virtual void S_splineFeasibleSmart() {}
    virtual void S_splineFeasibleSmart3(int index) {(void)index;}
    virtual void optimizeCurve(double** nodes, int numNodes,InitDir initDir = Smart) {(void)nodes; (void)numNodes; (void)initDir;}
    virtual void fineOptimizeCurve(double** nodes, int numNodes, int depth=10) {(void)nodes; (void)numNodes; (void)depth;}
    virtual int S_Curve(double *a, double *b, double *rx, double *ry);
    virtual Node pointInTheMiddle(double*, double*) {return Node();}
    virtual int getForm(int) {return 0;}

    void subdivide();
    void endSubdivide();

    void copy(Curve* curve);

    double createDataBallCurve(Curve* databallCurve);


    void convertCirclesToDataBalls();

    double getCurveBE();

    void clearExceptPoints();

    int furthest_marked();
    virtual double runCommand(int c1, double c2 = 0);

    static const int coarseRes = 360;
    static const int coarseResGamma = 120;
    int pointChanged(int index);
    bool *respectCorner;
    int count, phantCount, count_temp;
    bool closed;
    int subCount;

    string curveMethod;

};

#endif // CURVE_H
