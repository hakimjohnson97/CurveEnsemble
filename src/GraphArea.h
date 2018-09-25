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


#ifndef GRAPH_AREA_H
#define GRAPH_AREA_H

#include <QWidget>
#include <QFrame>
#include <QScrollArea>
#include <QContextMenuEvent>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QFormLayout>
#include <QPushButton>
#include <QComboBox>
#include <QStackedWidget>
#include <QLabel>
#include <QImage>
#include "Specific_Curves/Curve.h"
#include "Specific_Curves/ResElasticSpline.h"
#include "Specific_Curves/CubicQuasiElasticSpline.h"
#include "Specific_Curves/ElasticSpline.h"
#include "Specific_Curves/CubicSpline.h"
#include "Specific_Curves/ScottSpline1.h"
#include "Specific_Curves/CondG2Cubic.h"
#include "Specific_Curves/JZSpline.h"
#include "Specific_Curves/C2ParCubic.h"
#include "Specific_Curves/Lines.h"
#include "RectWidget.h"
#include "OptimizeWorker.h"
#include <QThread>
#include <QMutex>
#include <QDebug>


extern QMutex mutex;

//Determines different ways an indivdual node will show itself and the arrow coming out of it
enum PointStyle
{
    DisplayPoint = 0,
    DisplayPointWithArrow,
    DisplayPointWithTail,
    DisplayPointWithArrowAndTail,
    DisplayArrow,
    DisplayReverseArrow,
    NoDisplay,
    DisplayTick,
    PointStyleLength,
};


//A simple struct that contains the data related to the display of a node and its piece (curve).
struct NodeDisplay
{
      NodeDisplay() { pieceWidth = 1; pieceColor = Qt::red; showPiece = true;
       pointRadius = 5; pointColor = Qt::blue; pointStyle = DisplayPoint;};
      double pieceWidth;
      QColor pieceColor;
      bool showPiece;
      double pointRadius;
      QColor pointColor;
      PointStyle pointStyle;
      void operator= (NodeDisplay d2)
      {
          pieceWidth = d2.pieceWidth;
          pieceColor = d2.pieceColor;
          showPiece = d2.showPiece;
          pointRadius = d2.pointRadius;
          pointColor = d2.pointColor;
          pointStyle = d2.pointStyle;
      };
};

/*
  This class is the display space that draws the curves and allows the user to interact with the curves themselves.
  Users can select nodes and their corresponding pieces and move them, add new nodes and much more.
  GraphArea holds an array of curves of type Curve. (note this is the abtract class; each curve is a specific curve method)
*/
class GraphArea : public QFrame
{
  Q_OBJECT
  public:
      enum GraphMode {Normal = 0, Zoom, Transform, ClipImageRect};
      enum ClickMode {Inactive = 0, MovePoint, MovePivot, MoveArrow, MoveTail, DragScreen, InsertPoint, MoveReflectLine1, MoveReflectLine2, MoveClipRect1, MoveClipRect2, MoveClipRect3, MoveClipRect4, MoveClipRect5};
      enum CurveMethod {MethodResElasticSpline, MethodCubicQuasiElasticSpline, MethodElasticSpline, MethodCubicSpline,  MethodScottSpline1, MethodCondG2Cubic, MethodJZSpline, MethodSeparator, MethodLines, MethodCircles, MethodDataBalls, MethodC2ParCubic, MethodClassicalCubicSpline};

      GraphArea (QWidget *parent = 0);
      ~GraphArea();
      double getCurveCount() {return count;};
      double getPointCount(int index) {return curve[index]->getCount();};
      double getScale() {return scale;};
      int getSelectedPointIndex() {return p_selected;};
      int getSelectedCurveIndex() {return selected;};
      Curve* getSelectedCurve() {return curve[selected];};
      Curve* getCurve(int i) {if (i<count)return curve[i];};
      Curve::Node getNode(int index) {return curve[selected]->getNode(index);};
      NodeDisplay getSelectedNodeDisplay() {return nodeDisplay[selected] [p_selected];};

      NodeDisplay getNodeDisplay(int index) {return nodeDisplay[selected] [index];};
      NodeDisplay getActualNodeDisplay(int curveIndex, int pointIndex);
      void setSelectedNodeDisplay(NodeDisplay display) { nodeDisplay[selected] [p_selected] = display; update();};
      void setNodeDisplay(int index, NodeDisplay display) { nodeDisplay[selected] [index] = display; update();};
      void setUseNodeDisplay(bool use) {useNodeDisplay = use;};
      bool isNodeDisplayUsed() {return useNodeDisplay;};
      void setSelection(int);
      void centerPoint(int x, int y);
      void centerNode(int index);
      void centerCurve(int index);
      QRect getClipRect() {return clipRect;};
      void setGraphMode(GraphMode);
      GraphMode getGraphMode();
      void setCurveMethod(CurveMethod, int i = -1);
      CurveMethod getCurveMethod() {return (CurveMethod)curveMethod[selected];}

      void save(QString fileName);
      void open(QString fileName);

      void setBackgroundColor(QColor color) {setPalette(QPalette(color));};
      QColor getBackgroundColor () {return palette().color(QPalette::Window);};

      QColor backgroundColor;
      QVector< QList<NodeDisplay> > nodeDisplay;


      //Paint Stuff
      void drawCurves(QPainter&, bool translate = true);
      void saveToImage(QString fileName);
      void keyPress(QKeyEvent*);

      void setOnFlyOptimization(bool fly) {onFly = fly;}
      void loadBackgroundImage();

      void anchorSubdivision();

  protected:
      void mousePressEvent(QMouseEvent*);
      void mouseMoveEvent(QMouseEvent*);
      void mouseReleaseEvent(QMouseEvent*);
      void mouseDoubleClickEvent(QMouseEvent*);

      void paintEvent(QPaintEvent*);
      void keyPressEvent(QKeyEvent*);
  signals:
      void selectionChanged();
      void mouseCoordsChanged(int, int);
      void directionChanged(double);
      void positionChanged(double, double);
      void nodeChanged();
      void cleared();
      void transformWindowClosed();
      void optimizeOnFlyThread();
      void nodeAdded();
  public slots:
      void addCurve(CurveMethod method = MethodResElasticSpline, Curve* curve = NULL);
      void changeCurve(int);
      void removeCurve(int);
      void clear();
      void addNode();
      void setScale(double);
      void saveAsText();
      void loadFromText();
      void saveToClipboard();
      void loadFromClipboard();
      void optimize(Curve::InitDir);
      void fineOptimize();
      void transform();

      void setReflectLine1();
      void setReflectLine2();
      void centerClipRect();
      void endTransformMode();
      void optimizeOnFly();
      void finishOptimizeOnFly();

  private:
      void createTransformWindow();
      void setDefaultCursor();
      void newCurve(int i, CurveMethod method = MethodResElasticSpline);

      Curve** curve;
      QString savefile;
      int count, selected, p_selected;
      CurveMethod *curveMethod;

      GraphMode graphMode;
      ClickMode clickMode;
      QPoint clickPoint;
      int clickIndex;
      QPoint pivot, reflectLinePos1, reflectLinePos2;
      bool useNodeDisplay;

      double scale;
      QPoint viewOffset;

      //Transform Window
      QDialog *transformDialog;
      QComboBox *transformType;
      //Translate
      QDoubleSpinBox *translateX, *translateY;
      //Rotate
      QDoubleSpinBox *rotateAngle;
      //Resize
      QDoubleSpinBox *scaleFactor;
      //Pivot
      QDoubleSpinBox *pivotX, *pivotY;
      //Reflect
      QDoubleSpinBox *reflectLine1X, *reflectLine1Y, *reflectLine2X, *reflectLine2Y;
      //Transform
      QDoubleSpinBox *c1real, *c1imag;
      QDoubleSpinBox *c2real, *c2imag;
      QCheckBox *reflection;

      QRect clipRect;

      QPushButton *transformButton;

      QThread optimizeThread;
      OptimizeWorker *optimizeWorker;

      QImage *backgroundImage;

      bool onFly;

};



#endif
