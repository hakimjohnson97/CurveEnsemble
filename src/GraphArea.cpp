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


#include <QApplication>
#include <QWidget>
#include <QMouseEvent>
#include <QPaintEvent>
#include <QPainter>
#include <QPushButton>
#include <QGridLayout>
#include <QTime>
#include <cmath>
#include <ctime>
#include <QFrame>
#include <QKeyEvent>
#include <QFileDialog>
#include <QScrollArea>
#include <QPicture>
#include <fstream>
#include <QDebug>
#include "GraphArea.h"
using namespace std;

#define ARROW_LENGTH 8
#define PREAMBLE "Curve Ensemble Project File. The version number is stored as a long int multiplied by a thousand. For example, version 1.075 corresponds to 1075."
#define VERSION 1.100


const int POINT_DIAMETER = 10;
const int POINT_RADIUS = POINT_DIAMETER/2;
const int SCREEN_WIDTH = 1000;
const int SCREEN_HEIGHT = 1000;

double test;
fstream file2;

#define X 0
#define Y 1
#define DIR 2
#define CORNER 3
#define NODEDIR 3
#define STATE 4
#define BE 5


const double pi=3.14159265358979323846; 
const double degrees = pi/180;  



double cpointsx[600];
double cpointsy[600];
QPointF cpoints[600];

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
  u[X]=u[X]-u[STATE]*cos(u[NODEDIR]*degrees); u[Y]=u[Y]-u[STATE]*sin(u[NODEDIR]*degrees);
  u[CORNER]=0; u[STATE]=0;
}// if
}// databall2node--------------------------------------------------------------------------------------

inline void node2databall(double *u, double *u_save){
// node2databall(u,u_save) : If u_save[STATE]<0, then we copy 
// u_save[X,Y,NODEDIR,STATE] into u
if (u_save[STATE]<0 || u_save[STATE]>100)
  {u[X]=u_save[X]; u[Y]=u_save[Y]; u[NODEDIR]=u_save[NODEDIR]; u[STATE]=u_save[STATE];}
}// node2databall--------------------------------------------------------------------------------------



double DistBetweenPoints(double x1, double y1, double x2, double y2)
{
       return sqrt(pow(x2-x1, 2) + pow(y2-y1, 2));
}

void PolarProjection(double* rx, double *ry, double x, double y, int dist, double dir)
{
      *rx = x + dist*cos(dir);
      *ry = y + dist*sin(dir);
}

QPoint PolarProjection(double x, double y, int dist, double dir)
{
      return QPoint(x + dist*cos(dir), y + dist*sin(dir));
}

GraphArea::GraphArea (QWidget *parent) : QFrame(parent)
{
                     
    nodeDisplay.append(QList<NodeDisplay>());

    setGeometry(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
    file2.open("stderrr.txt", ios::out);
    curve = new Curve*[1000];
    
    setPalette(QPalette(QColor(255, 255, 255)));
    setAutoFillBackground(true);
    setFocusPolicy(Qt::StrongFocus);   
    setMouseTracking(true);

    setFrameStyle(QFrame::Box | QFrame::Raised);
    setLineWidth(2);
    setFrameRect(geometry());
    scale = 1;
    viewOffset = QPoint(0, 0);
    count = 1;
    graphMode = Normal;
    clickMode = Inactive;
    selected = 0;
    p_selected = -1;
    useNodeDisplay = false;

    curveMethod = new CurveMethod[100];
    curveMethod[0] = MethodElasticSpline;
    
    createTransformWindow();
    
    clipRect.setCoords(-200, -200, 200, 200);

    //Optimize on the fly business
    optimizeWorker = new OptimizeWorker;
    optimizeWorker->moveToThread(&optimizeThread);
    connect(&optimizeThread, SIGNAL(finished()), optimizeWorker, SLOT(deleteLater()));
    connect(this, SIGNAL(optimizeOnFlyThread()), optimizeWorker, SLOT(optimizeOnFly()));
    connect(optimizeWorker, SIGNAL(optimizeFinished()), this, SLOT(finishOptimizeOnFly()));
    connect(this, SIGNAL(positionChanged(double, double)), this, SLOT(optimizeOnFly()));
    optimizeThread.start();

    onFly = true;
    backgroundImage = new QImage();

    newCurve(0);
}

void GraphArea::createTransformWindow()
{
     transformDialog = new QDialog(this);
     
     transformType = new QComboBox();
     transformType->addItem("Translate");
     transformType->addItem("Rotate");
     transformType->addItem("Resize");
     transformType->addItem("Rotate and Resize");
     transformType->addItem("Reflect");
     transformType->addItem("Similarity Transform");      
     
     QStackedWidget *transformViews = new QStackedWidget(); 
     connect(transformType, SIGNAL(currentIndexChanged(int)), transformViews, SLOT(setCurrentIndex(int))); 
     connect(transformType, SIGNAL(currentIndexChanged(int)), this, SLOT(update())); 
     
     //Translate
     QLabel *translateLabel = new QLabel("Translate By:");
     translateX = new QDoubleSpinBox();
     translateY = new QDoubleSpinBox(); 
     
     QHBoxLayout* translateLayout = new QHBoxLayout();
     translateLayout->addWidget(translateLabel); 
     translateLayout->addWidget(translateX); 
     translateLayout->addWidget(translateY); 
     QWidget *translateView = new QWidget();
     translateView->setLayout(translateLayout);
     transformViews->addWidget(translateView); 
     
     //Pivot
     QLabel *pivotLabel = new QLabel("Pivot:");
     pivotX = new QDoubleSpinBox();
     pivotY = new QDoubleSpinBox();
     QHBoxLayout* pivotLayout = new QHBoxLayout();
     pivotLayout->addWidget(pivotLabel);
     pivotLayout->addWidget(pivotX);
     pivotLayout->addWidget(pivotY);
     
     //Rotate
     QLabel *rotateLabel1 = new QLabel("Rotate by:");
     QLabel *rotateLabel2 = new QLabel("degrees");
     rotateAngle = new QDoubleSpinBox();
     QHBoxLayout *rotateLayout = new QHBoxLayout();
     rotateLayout->addWidget(rotateLabel1); 
     rotateLayout->addWidget(rotateAngle); 
     rotateLayout->addWidget(rotateLabel2);
     QWidget *rotateView = new QWidget();
     rotateView->setLayout(rotateLayout);
     transformViews->addWidget(rotateView); 
     
     //Resize
     QLabel *resizeLabel = new QLabel("Resize by scale factor:");
     scaleFactor = new QDoubleSpinBox();
     QHBoxLayout *resizeLayout = new QHBoxLayout();
     resizeLayout->addWidget(resizeLabel); 
     resizeLayout->addWidget(scaleFactor); 
     QWidget *resizeView = new QWidget();
     resizeView->setLayout(resizeLayout);
     transformViews->addWidget(resizeView); 
     
     //Rotate and Resize
     QWidget *rotateResizeView = new QWidget();
     transformViews->addWidget(rotateResizeView); 
     
     //Reflect
     QLabel *reflectLabel = new QLabel("Reflect about line:");
     QLabel *reflectLineLabel1 = new QLabel("Pos1:");
     reflectLine1X = new QDoubleSpinBox();
     reflectLine1Y = new QDoubleSpinBox();
     QLabel *reflectLineLabel2 = new QLabel("Pos2:");
     reflectLine2X = new QDoubleSpinBox();
     reflectLine2Y = new QDoubleSpinBox();
     QHBoxLayout *reflectLayout = new QHBoxLayout();
     reflectLayout->addWidget(reflectLabel);
     reflectLayout->addWidget(reflectLineLabel1);
     reflectLayout->addWidget(reflectLine1X);
     reflectLayout->addWidget(reflectLine1Y);
     reflectLayout->addWidget(reflectLineLabel2);
     reflectLayout->addWidget(reflectLine2X);
     reflectLayout->addWidget(reflectLine2Y);
     QWidget *reflectView = new QWidget();
     reflectView->setLayout(reflectLayout);
     transformViews->addWidget(reflectView); 
     
     //Similarity Transform
     c1real = new QDoubleSpinBox();
     c1imag = new QDoubleSpinBox();
     c2real = new QDoubleSpinBox();
     c2imag = new QDoubleSpinBox();
     c1real->setRange(-1000, 1000);
     c1imag->setRange(-1000, 1000);
     c2real->setRange(-1000, 1000);
     c2imag->setRange(-1000, 1000);
     c1real->setValue(1);
     c1imag->setValue(0);
     reflection = new QCheckBox();
     
     QHBoxLayout *c1Layout = new QHBoxLayout();
     QHBoxLayout *c2Layout = new QHBoxLayout();
     c1Layout->addWidget(c1real);
     c1Layout->addWidget(c1imag);
     c2Layout->addWidget(c2real);
     c2Layout->addWidget(c2imag);
     
     QFormLayout *transformLayout = new QFormLayout();
     transformLayout->addRow("C1:", c1Layout);
     transformLayout->addRow("C2:", c2Layout);
     transformLayout->addRow("With Reflection:", reflection);
     QWidget *transformView = new QWidget();
     transformView->setLayout(transformLayout);
     transformViews->addWidget(transformView);
     
     transformButton = new QPushButton("Transform");
     connect(transformButton, SIGNAL(clicked()), this, SLOT(transform()));
     connect(transformDialog, SIGNAL(rejected()), this, SLOT(endTransformMode()));
     connect(transformDialog, SIGNAL(rejected()), this, SIGNAL(transformWindowClosed()));
     
     QVBoxLayout *layout = new QVBoxLayout();
     layout->addWidget(transformType);
     layout->addWidget(transformViews);
     layout->addWidget(transformButton);
     transformDialog->setLayout(layout);  
}

void GraphArea::setDefaultCursor()
{
     if (graphMode == Normal || graphMode == Transform || graphMode == ClipImageRect)
         setCursor(Qt::ArrowCursor);
     else
         setCursor(Qt::CrossCursor);
}

void GraphArea::transform()
{
    if (transformType->currentIndex() == 0)
        curve[selected]->moveCurve(translateX->value(), translateY->value());
    else if (transformType->currentIndex() == 1)    {
        curve[selected]->rotate(rotateAngle->value(), pivot.x(), pivot.y());
        if (curveMethod[selected] == MethodClassicalCubicSpline || curveMethod[selected] == MethodC2ParCubic)
            optimizeOnFly();
    }
    else if (transformType->currentIndex() == 2)
        curve[selected]->resize(scaleFactor->value(), pivot.x(), pivot.y());
    else if (transformType->currentIndex() == 4)   {
        curve[selected]->reflect(reflectLinePos1.x(), reflectLinePos1.y(), reflectLinePos2.x(), reflectLinePos2.y());
        optimizeOnFly();
    }
    else if (transformType->currentIndex() == 5)    {
        curve[selected]->transform(c1real->value(), c1imag->value(), c2real->value(), c2imag->value(), reflection->isChecked());
        if (reflection->isChecked())
            optimizeOnFly();
    }
    update();
}

void GraphArea::setPivot()
{
}

void GraphArea::setReflectLine1()
{
    reflectLinePos1 = QPoint(reflectLine1X->value(), reflectLine1Y->value());
}

void GraphArea::setReflectLine2()
{
     reflectLinePos2 = QPoint(reflectLine2X->value(), reflectLine2Y->value());
}

void GraphArea::mousePressEvent(QMouseEvent* event)
{
     Curve::Node node;
     NodeDisplay display;
     QPoint pos;
     int eventX = event->x()/scale + viewOffset.x();
     int eventY = event->y()/scale + viewOffset.y();
     int i;
     if (event->button() == Qt::LeftButton)    {
        // if (graphMode == Normal)    {
             if (clickMode == Inactive)    {
                 if (graphMode == ClipImageRect)    {
                         if (DistBetweenPoints(clipRect.x(), clipRect.y(), eventX, eventY) <= 5)    {
                             clickMode = MoveClipRect1;
                             setCursor(Qt::SizeFDiagCursor);
                         }
                         else if (DistBetweenPoints(clipRect.x() + clipRect.width(), clipRect.y(), eventX, eventY) <= 5)    {
                             clickMode = MoveClipRect2;
                             setCursor(Qt::SizeBDiagCursor);
                         }
                         else if (DistBetweenPoints(clipRect.x(), clipRect.y() + clipRect.height(), eventX, eventY) <= 5)    {
                             clickMode = MoveClipRect3;
                            setCursor(Qt::SizeBDiagCursor);
                         }
                         else if (DistBetweenPoints(clipRect.x() + clipRect.width(), clipRect.y() + clipRect.height(), eventX, eventY) <= 5)    {
                             clickMode = MoveClipRect4;
                             setCursor(Qt::SizeFDiagCursor);
                         }
                 }  
                 for (i = 0; i < curve[selected]->count; i++)    {
                     node = getNode(i);
                     display = getActualNodeDisplay(selected, i);
                     
                     if (DistBetweenPoints(node.x(), node.y(), eventX, eventY) <= display.pointRadius)    {
                         clickIndex = i;
                         clickMode = MovePoint;
                         setSelection(i);
                     }
                     pos = PolarProjection(node.x(), node.y(), display.pointRadius*ARROW_LENGTH, node.outDir()*degrees);
                     if ((display.pointStyle == DisplayPointWithArrow || display.pointStyle == DisplayPointWithArrowAndTail) 
                       && DistBetweenPoints(pos.x(), pos.y(), eventX, eventY) <= display.pointRadius*3)    {
                          clickIndex = i;
                          clickMode = MoveArrow;
                          setSelection(i);
                     }
                     pos = PolarProjection(node.x(), node.y(), display.pointRadius*ARROW_LENGTH, (node.outDir()-node.corner()-180)*degrees);
                     if ((display.pointStyle == DisplayPointWithTail || display.pointStyle == DisplayPointWithArrowAndTail) 
                       && DistBetweenPoints(pos.x(), pos.y(), eventX, eventY) <= display.pointRadius)    {
                          clickIndex = i;
                          clickMode = MoveTail;
                          setSelection(i);
                     }
                     
                     if (graphMode == Transform)    {
                         if (DistBetweenPoints(pivot.x(), pivot.y(), eventX, eventY) <= 5)    {
                             clickMode = MovePivot;
                         }
                         if (DistBetweenPoints(reflectLinePos1.x(), reflectLinePos1.y(), eventX, eventY) <= 5)    {
                             clickMode = MoveReflectLine1;
                         }
                         if (DistBetweenPoints(reflectLinePos2.x(), reflectLinePos2.y(), eventX, eventY) <= 5)    {
                             clickMode = MoveReflectLine2;
                         }
                     }
                 }
             }
                 if (clickMode == Inactive)
                         clickMode = DragScreen;
         //}
     }
     else
         if (graphMode != Zoom)
             clickMode = InsertPoint;
     clickPoint = QPoint(eventX, eventY);

}

void GraphArea::mouseMoveEvent(QMouseEvent* event)
{
     Curve::Node node;
     int eventX = event->x()/scale + viewOffset.x();
     int eventY = event->y()/scale + viewOffset.y();
     
     setDefaultCursor();
     
     //All graph modes
     if (clickMode == Inactive)    {
         if (graphMode == ClipImageRect)    {
             if (DistBetweenPoints(clipRect.x(), clipRect.y(), eventX, eventY) <= 5)    {
                 setCursor(Qt::SizeFDiagCursor);
             }
             else if (DistBetweenPoints(clipRect.x() + clipRect.width(), clipRect.y(), eventX, eventY) <= 5)    {
                 setCursor(Qt::SizeBDiagCursor);
             }
             else if (DistBetweenPoints(clipRect.x(), clipRect.y() + clipRect.height(), eventX, eventY) <= 5)    {
                 setCursor(Qt::SizeBDiagCursor);
             }
             else if (DistBetweenPoints(clipRect.x() + clipRect.width(), clipRect.y() + clipRect.height(), eventX, eventY) <= 5)    {
                 setCursor(Qt::SizeFDiagCursor);
             }
         }
     }
     else if (clickMode == DragScreen)     {
         viewOffset.setX(viewOffset.x() - (eventX - clickPoint.x()));
         viewOffset.setY(viewOffset.y() - (eventY - clickPoint.y()));
     }
     if (graphMode == Normal)    {
         if (clickMode == MovePoint)    {
             node = getNode(clickIndex);
             node.setPos(eventX, eventY);
             curve[selected]->setNode(clickIndex, node);
             emit positionChanged(eventX, eventY);
         }
         else if (clickMode == MoveArrow)    {
             node = getNode(clickIndex);
             node.setOutClamp(-QLineF(node.x(), node.y(), eventX, eventY).angle());

             curve[selected]->setNode(clickIndex, node);
             emit directionChanged(-QLineF(node.x(), node.y(), eventX, eventY).angle());
         }
         else if (clickMode == MoveTail)    {
             node = getNode(clickIndex);
             if (node.isRegular())
                 node.setCorner(node.outDir() - (-QLineF(node.x(), node.y(), eventX, eventY).angle()) - 180);
             else
                 node.setIncClamp((-QLineF(node.x(), node.y(), eventX, eventY).angle()) - 180);                         

             curve[selected]->setNode(clickIndex, node);
             emit selectionChanged();
         }
     }
     else if (graphMode == Transform)    {
         if (clickMode == MovePoint)    {
             node = getNode(clickIndex);
             if (transformType->currentIndex() == 0)
                 curve[selected]->moveCurve(eventX - node.x(), eventY - node.y());
             else if (transformType->currentIndex() == 1)    {
                 curve[selected]->rotate(clickIndex, eventX, eventY, pivot.x(), pivot.y());
             }
             else if (transformType->currentIndex() == 2)
                 curve[selected]->resize(clickIndex, eventX, eventY, pivot.x(), pivot.y());
             if (transformType->currentIndex() == 3)
                 curve[selected]->transformWithPivot(clickIndex, eventX, eventY, pivot.x(), pivot.y());

             if (curveMethod[selected] == MethodClassicalCubicSpline || curveMethod[selected] == MethodC2ParCubic)
                 optimizeOnFly();
            // setSelection(clickIndex);
         }
         else if (clickMode == MovePivot)    {
             QPoint pos = QPoint(eventX, eventY);
             Curve::Node node;
             for (int i = 0; i < curve[selected]->getCount(); i++)    {
                 node = getNode(i);
                 if (DistBetweenPoints(eventX, eventY, node.x(), node.y()) <= 5)
                      pos = QPoint(node.x(), node.y());
             }
             pivot = pos;
         }
         else if (clickMode == MoveReflectLine1)    {
             QPoint pos = QPoint(eventX, eventY);
             Curve::Node node;
             for (int i = 0; i < curve[selected]->getCount(); i++)    {
                 node = getNode(i);
                 if (DistBetweenPoints(eventX, eventY, node.x(), node.y()) <= 5)
                      pos = QPoint(node.x(), node.y());
             }
             reflectLinePos1 = pos;
         }
         else if (clickMode == MoveReflectLine2)    {
             QPoint pos = QPoint(eventX, eventY);
             Curve::Node node;
             for (int i = 0; i < curve[selected]->getCount(); i++)    {
                 node = getNode(i);
                 if (DistBetweenPoints(eventX, eventY, node.x(), node.y()) <= 5)
                      pos = QPoint(node.x(), node.y());
             }
             reflectLinePos2 = pos;
         }
     }
     else if (graphMode == ClipImageRect)    {
         if (clickMode == MoveClipRect1)    {
             clipRect.setTopLeft(QPoint(eventX, eventY));
         }
         else if (clickMode == MoveClipRect2)    {
             clipRect.setTopRight(QPoint(eventX, eventY));
         }
         else if (clickMode == MoveClipRect3)    {
             clipRect.setBottomLeft(QPoint(eventX, eventY));
         }
         else if (clickMode == MoveClipRect4)    {
             clipRect.setBottomRight(QPoint(eventX, eventY));
         }
     }
         
     if (clickMode != DragScreen)
         emit mouseCoordsChanged(eventX, eventY);
     update();
}

void GraphArea::mouseReleaseEvent(QMouseEvent* event)
{
     Curve::Node node;
     int i;
     int eventX = event->x()/scale + viewOffset.x();
     int eventY = event->y()/scale + viewOffset.y();
     if (event->button() == Qt::LeftButton)    {
         if (clickMode == DragScreen)
             clickMode = Inactive;
         if (graphMode == Normal)    {
             if (clickMode == MoveArrow)    {
                 mouseMoveEvent(event);
                 clickMode = Inactive;
             }
             else if (clickMode == MoveTail)    {
                  node = getNode(clickIndex);
                 if (node.isRegular())
                     node.setCorner(node.outDir() - (-QLineF(node.x(), node.y(), eventX, eventY).angle()) - 180);
                 else
                    node.setIncClamp((-QLineF(node.x(), node.y(), eventX, eventY).angle()) - 180);
                 curve[selected]->setNode(clickIndex, node);
                 curve[selected]->S_splineFeasibleSmart();
                 clickMode = Inactive;
             }
             else if (clickMode == MovePoint)    {
                 node = getNode(clickIndex);
                 if (DistBetweenPoints(eventX, eventY, node.x(),node.y()) > POINT_RADIUS)
                     node.setPos(eventX, eventY);
                 setSelection(clickIndex);
                 curve[selected]->setNode(clickIndex, node);
                 clickMode = Inactive;
             }
         }
         else if (graphMode == Transform)    {
             if (clickMode == MovePoint)    {

                 clickMode = Inactive;
             }
             else if (clickMode == MovePivot)    {
                 clickMode = Inactive;
             }
             else if (clickMode == MoveReflectLine1)    {
                 clickMode = Inactive;
             }
             else if (clickMode == MoveReflectLine2)    {
                 clickMode = Inactive;
            }
         }
         else if (graphMode == Zoom)    {
             scale *= 2;
             centerPoint(eventX, eventY);
         }
         else if (graphMode == ClipImageRect)    {
             if (clickMode == MoveClipRect1)    {
                 clipRect.setTopLeft(QPoint(eventX, eventY));
                 clickMode = Inactive;
             }
             else if (clickMode == MoveClipRect2)    {
                 clipRect.setTopRight(QPoint(eventX, eventY));
                 clickMode = Inactive;
             }
             else if (clickMode == MoveClipRect3)    {
                 clipRect.setBottomLeft(QPoint(eventX, eventY));
                 clickMode = Inactive;
             }
             else if (clickMode == MoveClipRect4)    {
                 clipRect.setBottomRight(QPoint(eventX, eventY));
                 clickMode = Inactive;
             }
         }
     }
     else if (event->button() == Qt::RightButton)    {
         if (graphMode == Zoom)    {
             scale /= 2;
             centerPoint(eventX, eventY);
         }
         else if (clickMode == InsertPoint)    {
             if (DistBetweenPoints(clickPoint.x(), clickPoint.y(), eventX, eventY) < 5)    {
                                                   
                                                   
                 nodeDisplay[selected].append(NodeDisplay());  
                 if (p_selected != curve[selected]->getCount()-1)    {
                     for (i = curve[selected]->getCount(); i > p_selected+1; i--)    { 
                         nodeDisplay[selected][i] = nodeDisplay[selected][i-1];
                     }


                 }
                 if (p_selected >= 0)
                     nodeDisplay[selected][p_selected+1] =  nodeDisplay[selected][p_selected];
                 else
                      nodeDisplay[selected][p_selected+1] = NodeDisplay();
                                                   
                 curve[selected]->insertNode(p_selected, Curve::Node(eventX, eventY));  
                 setSelection(p_selected+1);
                 emit nodeAdded();
             }
             clickMode = Inactive;
         }
     }
     if (onFly)
         optimizeOnFly();
     update();
}

void GraphArea::mouseDoubleClickEvent(QMouseEvent* event)
{
     int eventX = event->x()/scale + viewOffset.x();
     int eventY = event->y()/scale + viewOffset.y();
     if (event->button() == Qt::RightButton)    {
         if (graphMode == Zoom)    {
             scale = 2; 
             centerPoint(eventX, eventY);
         }
     }
}

void GraphArea::contextMenuEvent(QContextMenuEvent* event)
{

}


void GraphArea::keyPressEvent(QKeyEvent* event)
{
    keyPress(event);
}

void GraphArea::keyPress(QKeyEvent* event)
{
    if (event->key() == Qt::Key_Delete)    {
         curve[selected]->deleteNode(p_selected);
         if (curve[selected]->getCount() <= 0)   {
             p_selected = -1;
         }
         else    {
             optimizeOnFly();
             setSelection(p_selected-1);
         }

         update();
     }
     if (event->key() == Qt::Key_Left || event->key() == Qt::Key_Down)    {//Plus without shift
         setSelection(p_selected-1);
         update();
     }
     if (event->key() == Qt::Key_Right || event->key() == Qt::Key_Up)    {
         setSelection(p_selected+1);
         update();
     }
     if (event->key() == Qt::Key_F2)    {
         optimize(Curve::Smart);
     }
}

void GraphArea::paintEvent(QPaintEvent* event)
{
     QFrame::paintEvent(event);
     QPainter painter(this);

     drawCurves(painter);

     painter.end();

}

void GraphArea::drawCurves(QPainter &painter, bool translate)
{
     
     painter.scale(scale, scale);
     if (translate)
         painter.translate(-viewOffset.x(), -viewOffset.y());
     painter.setRenderHint( QPainter::Antialiasing, true );
     painter.setBackground(Qt::red);

     QRectF rectangle;
     double dirx, diry;
     
     int p_count;
     int i, j, a, num_points;
     double u_save[6], v_save[6];

     painter.drawImage(QPointF(0, 0), *backgroundImage);
     
     for (a = 0; a < count; a++)    {
     p_count = curve[a]->count;
     
     NodeDisplay display;
     Curve::Node node;




     for (i = 0; i < p_count; i++)    {


         display = getActualNodeDisplay(a, i);

         painter.setPen(QPen(display.pieceColor, 1, Qt::SolidLine, Qt::RoundCap));
         painter.setBrush(QBrush(display.pieceColor, Qt::SolidPattern));

         if (curve[a]->getSubdivideCount() != 0 && (i != p_count-1 || curve[a]->isClosed()) && display.pointStyle != NoDisplay)    {

             QList<Curve::Node> nodes = curve[a]->getPhantNodes(i);
             for (int k = 0; k < nodes.size(); k++)    {
                 rectangle = QRect(nodes[k].x()-display.pointRadius/2, nodes[k].y()-display.pointRadius/2,  display.pointRadius, display.pointRadius);
                 painter.drawEllipse(rectangle);
             }
         }
         node = curve[a]->getNode(i);
         

            
         
         //Paint Curves
         if (display.showPiece && (i != p_count-1 || curve[a]->isClosed()))    {
        
             display = getActualNodeDisplay(a, i);
             painter.setPen(QPen(display.pieceColor, display.pieceWidth, Qt::SolidLine, Qt::RoundCap));

             if (curveMethod[a] != MethodCircles && curveMethod[a] != MethodDataBalls)    {
                 for (int k = 1; k < curve[a]->getPhantNodesPerPiece()+2; k++)    {
                     num_points = curve[a]->pointsOnCurve(i, cpointsx, cpointsy, k);
                     for (j = 0; j < num_points; j++)
                         cpoints[j] = QPointF(cpointsx[j], cpointsy[j]);
                     painter.drawPolyline(cpoints, num_points);
                 }
             }

             else if (curveMethod[a] == MethodCircles && (i % 2 == 0))    {
                 painter.setBrush(QBrush());
                 Curve::Node node2 = curve[a]->getNode(i+1);
                 QPoint centre = QPoint((node.x()+node2.x())/2, (node.y()+node2.y())/2);
                 int dist = DistBetweenPoints(node.x(), node.y(), node2.x(), node2.y())/2;
                 painter.drawEllipse(centre, dist, dist);
                 painter.setBrush(QBrush(display.pointColor, Qt::SolidPattern));
             }

         }



             painter.setPen(QPen(display.pointColor, 1, Qt::SolidLine, Qt::RoundCap));
             painter.setBrush(QBrush(display.pointColor, Qt::SolidPattern));

             if (display.pointStyle == DisplayPoint || display.pointStyle == DisplayPointWithArrow
                || display.pointStyle == DisplayPointWithTail || display.pointStyle == DisplayPointWithArrowAndTail
                || display.pointStyle == DisplayTick)    {
                 rectangle = QRect(node.x()-display.pointRadius, node.y()-display.pointRadius, display.pointRadius*2, display.pointRadius*2);
                 painter.drawEllipse(rectangle);
                 if (node.isDataBall())    {
                    painter.setBrush(QBrush());
                    double rad = node.dataBallRadius();
                    rectangle = QRect(node.x()-rad, node.y()-rad, rad*2, rad*2);
                    painter.drawEllipse(rectangle);
                    if (curveMethod[a] != MethodDataBalls)    {
                        double arrowTipX, arrowTipY;
                        PolarProjection(&arrowTipX, &arrowTipY, node.x(), node.y(), node.dataBallRadius(), node.corner()*degrees);
                        painter.drawLine(node.x(), node.y(), arrowTipX, arrowTipY);
                    }
                 }
             }
             if (display.pointStyle == DisplayPointWithArrow || display.pointStyle == DisplayPointWithArrowAndTail)    {
                 double arrowTipX, arrowTipY;
                 PolarProjection(&arrowTipX, &arrowTipY, node.x(), node.y(), display.pointRadius*ARROW_LENGTH, node.outDir()*degrees);
                 painter.drawLine(node.x(), node.y(), arrowTipX, arrowTipY);
                 PolarProjection(&dirx, &diry, arrowTipX, arrowTipY, display.pointRadius*2, (node.outDir()+150)*degrees);
                 painter.drawLine(arrowTipX, arrowTipY, dirx, diry);
                 PolarProjection(&dirx, &diry, arrowTipX, arrowTipY, display.pointRadius*2, (node.outDir()-150)*degrees);
                 painter.drawLine(arrowTipX, arrowTipY, dirx, diry);
                 
             }
  
             if (display.pointStyle == DisplayPointWithTail || display.pointStyle == DisplayPointWithArrowAndTail)    {
                 double arrowTipX, arrowTipY;
                 double dir = node.incDir()-180;
                 PolarProjection(&arrowTipX, &arrowTipY, node.x(), node.y(), display.pointRadius*ARROW_LENGTH, dir*degrees);
                 painter.drawLine(node.x(), node.y(), arrowTipX, arrowTipY);
                 PolarProjection(&dirx, &diry, arrowTipX, arrowTipY, display.pointRadius*2, (dir+30)*degrees);
                 painter.drawLine(arrowTipX, arrowTipY, dirx, diry);
                 PolarProjection(&dirx, &diry, arrowTipX, arrowTipY, display.pointRadius*2, (dir-30)*degrees);
                 painter.drawLine(arrowTipX, arrowTipY, dirx, diry);
                 
             }
             
             if (display.pointStyle == DisplayArrow || display.pointStyle == DisplayReverseArrow)    {
                 double arrowDir;
                 if (display.pointStyle == DisplayArrow)
                     arrowDir = node.outDir();
                 else
                     arrowDir = node.outDir() - 180;
                 PolarProjection(&dirx, &diry, node.x(), node.y(), display.pointRadius*1.5, (arrowDir+90)*degrees);
                 QPoint ArrowP1 = QPoint(dirx, diry);
                 PolarProjection(&dirx, &diry, node.x(), node.y(), display.pointRadius*1.5, arrowDir*degrees);
                 QPoint ArrowP2 = QPoint(dirx, diry);
                 PolarProjection(&dirx, &diry, node.x(), node.y(), display.pointRadius*1.5, (arrowDir-90)*degrees);
                 QPoint ArrowP3 = QPoint(dirx, diry);
                 
                 painter.drawPolygon(QPolygon() << ArrowP1 << ArrowP2 << ArrowP3);
             }
             if (display.pointStyle == DisplayTick)    {
                 double p1X, p1Y, p2X, p2Y;
                 PolarProjection(&p1X, &p1Y, node.x(), node.y(), display.pointRadius*2, (node.outDir()+90)*degrees);
                 PolarProjection(&p2X, &p2Y, node.x(), node.y(), display.pointRadius*2, (node.outDir()-90)*degrees);
                 painter.drawLine(p1X, p1Y, p2X, p2Y);
             }
             if (a == selected && i == p_selected)    {
                 painter.setPen(Qt::black);
                 painter.setBrush(QBrush());
                 rectangle.setRect(node.x() - display.pointRadius*1.5, node.y() - display.pointRadius*1.5, 
                   display.pointRadius*3, display.pointRadius*3);
                 painter.drawEllipse(rectangle);
            }
             
      
     }
     
     if (graphMode == ClipImageRect)    {
         painter.setPen(QPen(Qt::black, 1, Qt::SolidLine, Qt::RoundCap));
         painter.setBrush(QBrush());
         painter.drawRect(clipRect);
     }

     }
     if (graphMode == Transform)    {
          if (transformType->currentIndex() < 4 && transformType->currentIndex() != 0)    {
              painter.setPen(QPen(QColor("green"), 1, Qt::SolidLine, Qt::RoundCap));
              painter.setBrush(QBrush(QColor("green"), Qt::SolidPattern));
              rectangle = QRect(pivot.x()-5, pivot.y()-5, 10, 10);
              painter.drawEllipse(rectangle);
          }
          
          if (transformType->currentIndex() == 4)    {
              painter.setPen(QPen(QColor("cyan"), 1, Qt::SolidLine, Qt::RoundCap));
              painter.setBrush(QBrush(QColor("cyan"), Qt::SolidPattern));
              rectangle = QRect(reflectLinePos1.x()-5, reflectLinePos1.y()-5, 10, 10);
              painter.drawEllipse(rectangle);
              rectangle = QRect(reflectLinePos2.x()-5, reflectLinePos2.y()-5, 10, 10);
              painter.drawEllipse(rectangle);
              painter.drawLine(reflectLinePos1, reflectLinePos2);
          }
     }
}

NodeDisplay GraphArea::getActualNodeDisplay(int curveIndex, int pointIndex)
{
    NodeDisplay display = nodeDisplay[curveIndex][pointIndex];
    Curve::Node node = curve[curveIndex]->getNode(pointIndex);
    if (!useNodeDisplay)    {
        if (node.isRegular())
            display.pointStyle = DisplayPointWithTail;
        else    {
            if (node.hasIncClamp() && !node.hasOutClamp())
                display.pointStyle = DisplayPointWithTail;
            if (!node.hasIncClamp() && node.hasOutClamp())
                display.pointStyle = DisplayPointWithArrow;   
            if (node.hasIncClamp() && node.hasOutClamp())
                display.pointStyle = DisplayPointWithArrowAndTail;
        }
        
        if (display.pointStyle == NoDisplay)
            display.pointStyle = DisplayPoint;
    
        if (curve[curveIndex]->getForm(pointIndex) == 0 || curve[curveIndex]->getForm(pointIndex) == 1)
            display.pieceColor = QColor("red");
        else if (curve[curveIndex]->getForm(pointIndex) == 2)
            display.pieceColor = QColor("purple");
            
        if (curveIndex != selected)    {
            display.pieceColor = Qt::black;
            display.pointColor = Qt::black;
        }
    }
    return display;     
}

void GraphArea::setSelection(int n)
{
     if (n >= curve[selected]->getCount())
         n = 0;
     else if (n < 0)
         n = curve[selected]->getCount()-1;
     p_selected = n;

     emit selectionChanged();
}

void GraphArea::centerPoint(int x, int y)
{
     viewOffset = QPoint(x - size().width()/2/scale, y - size().height()/2/scale);
     update();
}

void GraphArea::centerNode(int index)
{
     Curve::Node node = getNode(index);
     centerPoint(node.x(), node.y());
     update();
}

void GraphArea::centerCurve(int index)
{
     Curve::Node node;
     QPointF average;
     int i;
     for (i = 0; i < curve[selected]->getCount(); i++)    {
         node = curve[selected]->getNode(i);
         average.setX(average.x() + node.x());
         average.setY(average.y() + node.y());
     }
     average.setX(average.x()/i);
     average.setY(average.y()/i);
     
     centerPoint(average.x(), average.y());
}

void GraphArea::changeCurve(int index)
{

     selected = index;
     p_selected =  curve[selected]->getCount()-1;
     update();

}

void GraphArea::removeCurve(int index)
{
     int i;
     for (i = index; i < count-1; i++)    {
         curve[i] = curve[i+1];
         nodeDisplay[i] = nodeDisplay[i+1];
         curveMethod[i] = curveMethod[i+1];
    }

    count--;
    selected = 0;
    p_selected = 0;
    update();
}

void GraphArea::addCurve(CurveMethod method, Curve* arg_curve)
{
     
     newCurve(count, method);

     if (arg_curve != NULL)
         curve[count]->copy(arg_curve);
     nodeDisplay.append(QList<NodeDisplay>());
     

     for (int i = 0; i < curve[count]->getCount(); i++)
         nodeDisplay[count].append(NodeDisplay());

     
     count++;
     selected = count-1;
     if (curve[count-1]->getCount() == 0)
         p_selected = -1;
     else
        p_selected = 0;
     update();
}

void GraphArea::addNode()
{

}

void GraphArea::setScale(double arg_scale)
{
     scale = arg_scale;
     update();
}

void GraphArea::setGraphMode(GraphMode mode)
{
    if (graphMode == Transform && mode != Transform)
        transformDialog->done(true);
    if (mode == Normal)
        setCursor(Qt::ArrowCursor);
    if (mode == Transform)    {
        setCursor(Qt::ArrowCursor);
        pivot = viewOffset + QPoint(size().width()/2, size().height()/2);
        reflectLinePos1 = viewOffset + QPoint(size().width()/4, size().height()/4);
        reflectLinePos2 = viewOffset + QPoint(size().width()*3/4, size().height()*3/4);
        transformDialog->show();

    }
    else if (mode == Zoom)
        setCursor(Qt::CrossCursor);
    else if (mode == ClipImageRect)
        setCursor(Qt::ArrowCursor);
    graphMode = mode;
    update();
}

GraphArea::GraphMode GraphArea::getGraphMode()
{
    return graphMode;
}

void GraphArea::open(QString fileName)
{
     char emptySpace[1000];
     int i, j;
     QRgb rgb;

     char buf[strlen(PREAMBLE)];
     char preamble[] = PREAMBLE;

     fstream file;
     file.open(fileName.toStdString().c_str(), ios::in | ios::binary);
     file.read((char*)buf, strlen(PREAMBLE));
     long int version;
     file.read((char*)&version, sizeof(long int));
     bool correct = true;
     for (int i = 0; i < strlen(PREAMBLE); i++)
         if (buf[i] != preamble[i])

             correct = false;
     bool old = !(correct);
     if (old)
         qDebug() << "File opened is old version.";
     else
         qDebug() << "File opened is new version.";
     file.read((char*)&count, sizeof(long int));

     
     nodeDisplay.clear();

   
    file.read((char*)&useNodeDisplay, sizeof(bool));
    file.read(emptySpace, 1000);
     for (i = 0; i < count; i++)    {
        file.read((char*)(&curveMethod[i]), sizeof(int));
        newCurve(i, curveMethod[i]);
        nodeDisplay.append(QList<NodeDisplay>());

        curve[i]->load(file, old);

        for (j = 0; j < curve[i]->count; j++)    { 

             nodeDisplay[i].append(NodeDisplay());

             file.read((char*)&nodeDisplay[i][j].pieceWidth, sizeof(double));
             file.read((char*)&rgb, sizeof(QRgb)); 
             nodeDisplay[i][j].pieceColor.setRgb(rgb);

             file.read((char*)&nodeDisplay[i][j].showPiece, sizeof(bool));
             file.read((char*)&nodeDisplay[i][j].pointRadius, sizeof(double));
             file.read((char*)&rgb, sizeof(QRgb));
             nodeDisplay[i][j].pointColor.setRgb(rgb);

             file.read((char*)&nodeDisplay[i][j].pointStyle, sizeof(PointStyle)); 
             file.read(emptySpace, 66);
             
        } 
        file.read(emptySpace, 35);
     }

     short int x, y, width, height;
     file.read((char*)&x, 2);
     file.read((char*)&y, 2);
     file.read((char*)&width, 2);
     file.read((char*)&height, 2);
     clipRect.setRect(x, y, width, height);

     centerCurve(selected);

     selected = 0;
     p_selected = 0;
     string str = fileName.toStdString();
     int p;
     p = str.rfind("/");

     update();
}

void GraphArea::save(QString fileName)
{
     char emptySpace[1000];
     int i, j;
     for (i = 0; i < 1000; i++)
         emptySpace[i] = 0;
     QRgb rgb;
     char preamble[] = PREAMBLE;

         
     fstream file;
     
     file.open(fileName.toStdString().c_str(), ios::out | ios::binary);
     file.write((char*)preamble, strlen(PREAMBLE));
     long int version = VERSION*1000;
     file.write((char*)&version, sizeof(long int));

     file.write((char*)&count, sizeof(long int));
     file.write((char*)&useNodeDisplay, sizeof(bool));
     file.write(emptySpace, 1000);
     for (i = 0; i < count; i++)    {
        file.write((char*)(&curveMethod[i]), sizeof(int));
        curve[i]->save(file);
        for (j = 0; j < curve[i]->count; j++)    {
            file.write((char*)&nodeDisplay[i][j].pieceWidth, sizeof(double));

            rgb = nodeDisplay[i][j].pieceColor.rgb();
            file.write((char*)&rgb, sizeof(QRgb));
            file.write((char*)&nodeDisplay[i][j].showPiece, sizeof(bool));
            file.write((char*)&nodeDisplay[i][j].pointRadius, sizeof(double));

            rgb = nodeDisplay[i][j].pointColor.rgb();
            file.write((char*)&rgb, sizeof(QRgb));
            file.write((char*)&nodeDisplay[i][j].pointStyle, sizeof(PointStyle));
            file.write(emptySpace, 66);
        }
        file.write(emptySpace, 35);
     }

     int x1, y1, width1, height1;
     clipRect.getRect(&x1, &y1, &width1, &height1);
     short int x2, y2, width2, height2;
     x2 = (short int)x1; y2 = (short int)y1; width2 = (short int)width1; height2 = (short int)height1;
     file.write((char*)&x2, 2);
     file.write((char*)&y2, 2);
     file.write((char*)&width2, 2);
     file.write((char*)&height2, 2);


     file.close();
}


void GraphArea::saveAsText()
{
     QString fileName = QFileDialog::getSaveFileName(this, "Save As Text", "", "Text Files (*.txt);;All Files(*)");
     if (fileName == NULL)
         return;
         
     fstream file;
     file.open(fileName.toStdString().c_str(), ios::out);
     file << "#Curve Ensemble Text File\n";
     curve[selected]->saveAsText(file);
     file.close();
}

void GraphArea::loadFromText()
{
     QString fileName = QFileDialog::getOpenFileName(this, "Load From Text", "", "Text Files (*.txt);;All Files(*)");
     if (fileName == NULL)
         return;
     fstream file;
     file.open(fileName.toStdString().c_str(), ios::in);
     curve[selected]->loadFromText(file);
     nodeDisplay[selected].clear();
     for (int i = 0; i < curve[selected]->getCount(); i++)    {
         nodeDisplay[selected].append(NodeDisplay());
     }
     file.close();
     setSelection(0);
     centerCurve(selected);
     update();
}

void GraphArea::saveToClipboard()
{   
     fstream file;
     file.open("clipboard657346.txt", ios::out);
     curve[selected]->saveAsText(file);
     file.close();
}

void GraphArea::loadFromClipboard()
{
     
     fstream file;
     file.open("clipboard657346.txt", ios::in);
     curve[selected]->loadFromText(file);
     file.close();
     nodeDisplay[selected].clear();
     for (int i = 0; i < curve[selected]->getCount(); i++)    {
         nodeDisplay[selected].append(NodeDisplay());
     }
     setSelection(0);
     update();
}

void GraphArea::clear()
{
     int i;
     nodeDisplay.clear();
     for (i = 0; i < count; i++)
         delete curve[i];
     newCurve(0);
     nodeDisplay.append(QList<NodeDisplay>());
     count = 1;
     selected = 0;
     p_selected = -1;
     centerPoint(0, 0);
     update();
     emit cleared();
}

void GraphArea::optimize(Curve::InitDir initDir = Curve::Smart)
{

     curve[selected]->optimize(initDir);

     emit directionChanged(getNode(p_selected).outDir());
     update();
}


void GraphArea::fineOptimize()
{
     curve[selected]->fineOptimize(onFly);
     emit directionChanged(getNode(p_selected).outDir());
     update();
}

void GraphArea::centerClipRect()
{

     clipRect.setCoords(viewOffset.x() + 10, viewOffset.y() + 10, viewOffset.x() - 10 + size().width(), viewOffset.x() - 10 + size().height());
     update();
}

void GraphArea::setCurveMethod(CurveMethod method, int i)
{
    if (i == -1)
        i = selected;

    if ((curveMethod[selected] == MethodC2ParCubic || curveMethod[selected] == MethodClassicalCubicSpline) && (method != MethodC2ParCubic && method != MethodClassicalCubicSpline))
        curve[i]->clearExceptPoints();

    curveMethod[selected] = method;
    if (method == MethodResElasticSpline)    {
        ResElasticSpline *temp = new ResElasticSpline;
        temp->copy(curve[i]);
        curve[i] = temp;
    }
    else if (method == MethodCubicQuasiElasticSpline)    {
        CubicQuasiElasticSpline *temp = new CubicQuasiElasticSpline;
        temp->copy(curve[i]);
        curve[i] = temp;
    }
    else if (method == MethodCubicSpline)    {
        CubicSpline *temp = new CubicSpline;
        temp->copy(curve[i]);
        curve[i] = temp;
    }
    else if (method == MethodElasticSpline)    {
        ElasticSpline *temp = new ElasticSpline;
        temp->copy(curve[i]);
        curve[i] = temp;
    }
    else if (method == MethodScottSpline1)    {
        ScottSpline1 *temp = new ScottSpline1;
        temp->copy(curve[i]);
        curve[i] = temp;
    }
    else if (method == MethodCondG2Cubic)    {
        CondG2Cubic *temp = new CondG2Cubic;
        temp->copy(curve[i]);
        curve[i] = temp;
    }
    else if (method == MethodJZSpline)    {
        JZSpline *temp = new JZSpline;
        temp->copy(curve[i]);
        curve[i] = temp;
    }
    else if (method == MethodC2ParCubic)    {
        C2ParCubic *temp = new C2ParCubic(false);
        temp->copy(curve[i]);
        curve[i] = temp;

    }
    else if (method == MethodClassicalCubicSpline)    {
        C2ParCubic *temp = new C2ParCubic(true);
        temp->copy(curve[i]);
        curve[i] = temp;
    }
    else if (method == MethodLines)    {
        Lines *temp = new Lines;
        temp->copy(curve[i]);
        curve[i] = temp;
    }
    else if (method == MethodDataBalls)    {
        Curve *temp = new Curve("DataBalls");
        temp->copy(curve[i]);
        curve[i] = temp;
    }
    else    {
        Curve *temp = new Curve;
        temp->copy(curve[i]);
        curve[i] = temp;
    }
}

void GraphArea::newCurve(int i, CurveMethod method)
{
    curveMethod[i] = method;

    if (method == MethodResElasticSpline)    {
       curve[i] = new ResElasticSpline;
    }
    else if (method == MethodCubicSpline)    {
       curve[i] = new CubicSpline;
    }
    else if (method == MethodCubicQuasiElasticSpline)    {
       curve[i] = new CubicQuasiElasticSpline;
    }
    else if (method == MethodElasticSpline)    {
       curve[i] = new ElasticSpline;
    }
    else if (method == MethodScottSpline1)    {
       curve[i] = new ScottSpline1;
    }
    else if (method == MethodCondG2Cubic)    {
       curve[i] = new CondG2Cubic;
    }
    else if (method == MethodJZSpline)    {
       curve[i] = new JZSpline;
    }
    else if (method == MethodC2ParCubic)    {
        curve[i] = new C2ParCubic(false);
    }
    else if (method == MethodClassicalCubicSpline)    {
        curve[i] = new C2ParCubic(true);
    }
    else if (method == MethodLines)    {
        curve[i] = new Lines;
    }
    else if (method == MethodDataBalls)    {
        curve[i] = new Curve("DataBalls");
    }
    else
        curve[i] = new Curve;

}

void GraphArea::endTransformMode()
{
    graphMode = Normal;
}

void GraphArea::optimizeOnFly()
{
    if (onFly)    {

        curve[selected]->optimize(Curve::Smart, true);

    }
}

void GraphArea::finishOptimizeOnFly()
{
    curve[selected]->copyFromTemp();
    emit directionChanged(getNode(p_selected).outDir());
    update();
}

void GraphArea::loadBackgroundImage()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Load Background");
    backgroundImage->load(fileName);
}

void GraphArea::anchorSubdivision()
{
     curve[selected]->anchorPhantomCurve();
     while (nodeDisplay[selected].size() < curve[selected]->getCount())
          nodeDisplay[selected].append(NodeDisplay());  
     
}

GraphArea::~GraphArea()
{
    optimizeThread.quit();
    optimizeThread.wait();
}
