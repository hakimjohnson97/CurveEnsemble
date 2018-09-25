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
#include <QTime>
#include <cmath>
#include <ctime>
#include <QFrame>
#include <QKeyEvent>
#include <QLineEdit>
#include <QCheckBox>
#include <fstream>
#include <QTabBar>
#include <QFileDialog>
#include <QLabel>
#include <QSlider>
#include <QScrollArea>
#include <QDebug>
#include "GraphUI.h"
using namespace std;

#define RANGE 1000000000


fstream file;



GraphUI::GraphUI (GraphArea *arg_graph, QWidget *parent) : QFrame(parent)
{
    graph = arg_graph;
                 
    setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    
    closed = new QCheckBox();
    
    posx = new QDoubleSpinBox();
    posx->setRange(-RANGE, RANGE);
    posy = new QDoubleSpinBox();
    posy->setRange(-RANGE, RANGE);
    
    regular = new QCheckBox();
    dataBall = new QCheckBox();
    
    incClamp = new QPushButton("Clamp");
    incClamp->setCheckable(true);
    incDir = new QDoubleSpinBox();
    incDir->setRange(-RANGE, RANGE);
    incDirLayout = new QHBoxLayout();
    incDirLayout->addWidget(incClamp);
    incDirLayout->addWidget(incDir);
    
    outClamp = new QPushButton("Clamp");
    outClamp->setCheckable(true);
    outDir = new QDoubleSpinBox();
    outDir->setRange(-RANGE, RANGE);
    outDirLayout = new QHBoxLayout();
    outDirLayout->addWidget(outClamp);
    outDirLayout->addWidget(outDir);
    
    corner = new QDoubleSpinBox();
    corner->setRange(-RANGE, RANGE);
    respectCorner = new QCheckBox("Respect");
    cornerLayout = new QHBoxLayout();
    cornerLayout->addWidget(corner);
    cornerLayout->addWidget(respectCorner);
    
    alphaMax = new QDoubleSpinBox();
    alphaMax->setRange(0, 360);

    lambda = new QDoubleSpinBox();
    lambda->setRange(0, 3);
    lambda->setSingleStep(0.05);

    lambda->setEnabled(false);

    curveMethod = new QComboBox();
    curveMethod->addItem("Restricted Elastic Spline");
    curveMethod->addItem("Quasi-Elastic Cubic");
    curveMethod->addItem("Elastic Spline");
    curveMethod->addItem("Cubic Spline");
    curveMethod->addItem("Yong-Cheng Cubic");
    curveMethod->addItem("Conditionally G2 Cubic");
    curveMethod->addItem("Jaklic-Zagar Cubic");
    curveMethod->insertSeparator(7);
    curveMethod->addItem("Lines");
    curveMethod->addItem("Circles");
    curveMethod->addItem("DataBalls");
    curveMethod->addItem("C2 Parametric Cubic Spline");
    curveMethod->addItem("Classical Cubic Spline");
    
    incDir->setEnabled(false);
    outDir->setEnabled(false);
    corner->setEnabled(true);
    
    connect(closed, SIGNAL(toggled(bool)), this, SLOT(handleChanges()));
    connect(posx, SIGNAL(editingFinished()), this, SLOT(handleChanges()));
    connect(posy, SIGNAL(editingFinished()), this, SLOT(handleChanges()));
    connect(incDir, SIGNAL(editingFinished()), this, SLOT(handleChanges()));
    connect(outDir, SIGNAL(editingFinished()), this, SLOT(handleChanges()));
    connect(corner, SIGNAL(editingFinished()), this, SLOT(handleChanges()));
    connect(respectCorner, SIGNAL(toggled(bool)), this, SLOT(handleChanges()));
    connect(regular, SIGNAL(toggled(bool)), this, SLOT(handleChanges()));
    connect(dataBall, SIGNAL(toggled(bool)), this, SLOT(handleChanges()));
    connect(incClamp, SIGNAL(toggled(bool)), this, SLOT(handleChanges()));
    connect(outClamp, SIGNAL(toggled(bool)), this, SLOT(handleChanges()));
    connect(alphaMax, SIGNAL(editingFinished()), this, SLOT(handleChanges()));
    connect(lambda, SIGNAL(editingFinished()), this, SLOT(handleChanges()));
    connect(curveMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(handleChanges()));

    QFormLayout *curveSpecific = new QFormLayout();
    curveSpecific->addRow("Closed:", closed);
    curveSpecific->addRow("Omega:", alphaMax);
    curveSpecific->addRow("Lambda:", lambda);
    curveSpecific->addRow("Curve Method:", curveMethod);

    nodeSpecific = new QFormLayout();
    nodeSpecific->addRow("Data Ball:", dataBall);
    nodeSpecific->addRow("Regularity:", regular);
    nodeSpecific->addRow("Incoming Direction:", incDirLayout);
    nodeSpecific->addRow("Outgoing Direction:", outDirLayout);
    nodeSpecific->addRow("Corner:", cornerLayout);
    nodeSpecific->addRow("Position X:", posx);
    nodeSpecific->addRow("Position Y:", posy);
   
    QVBoxLayout *layout = new QVBoxLayout;
    layout->addLayout(curveSpecific);
    layout->addSpacing(30);
    layout->addLayout(nodeSpecific);
    
    setLayout(layout);
    
    
    connect(graph, SIGNAL(selectionChanged()), this, SLOT(updateWidgets()));
    connect(graph, SIGNAL(directionChanged(double)), this, SLOT(setDirection(double)));
    connect(graph, SIGNAL(positionChanged(double, double)), this, SLOT(setPosition(double, double))); 

    
    userChange = true;
}

void GraphUI::updateWidgets()
{
    Curve::Node node = graph->getNode(graph->getSelectedPointIndex());
    
    userChange = false;

    //Update values
    dataBall->setChecked(node.isDataBall());
    closed->setChecked(graph->getSelectedCurve()->isClosed());
    incClamp->setChecked(node.hasIncClamp());
    outClamp->setChecked(node.hasOutClamp());
    outDir->setValue(node.outDir());
    corner->setValue(node.corner());
    respectCorner->setChecked(node.isCornerRespected());
    if (node.isDataBall())    {
        incDir->setValue(node.dataBallRadius());
        if (node.isTicked())
            regular->setChecked(true);
        else
            regular->setChecked(false);
    }
    else    {
        incDir->setValue(node.incDir());
        regular->setChecked(node.isRegular());
    }
    posx->setValue(node.x());
    posy->setValue(node.y());  
    userChange = true;


        
    incClamp->setEnabled(true);
    outClamp->setEnabled(true);
    incDir->setEnabled(true);
    outDir->setEnabled(true);
    corner->setEnabled(true);
    regular->setEnabled(true);
    dataBall->setEnabled(true);
    respectCorner->setEnabled(false); //Potentially bad
    
    if (node.isRegular())    {
        incClamp->setEnabled(false);
        outClamp->setEnabled(false);
        incDir->setEnabled(false);
        outDir->setEnabled(false);
    }
    else    {
        if (node.isDataBall())    {
            incDir->setEnabled(true);
            outDir->setEnabled(true);
            corner->setEnabled(true);
            incClamp->setEnabled(false);
            outClamp->setEnabled(false);
        }
        else     {
             corner->setEnabled(false);
             if (!node.hasIncClamp())    {
                 incDir->setEnabled(false);
             }
             if (!node.hasOutClamp())    {
                 outDir->setEnabled(false);
             }
             if (node.hasOutClamp() == true &&  node.hasIncClamp() == true)    {
                corner->setEnabled(true);
                respectCorner->setEnabled(true);
             }
        }
    }
    if (!graph->getSelectedCurve()->isClosed())    {
        if (graph->getSelectedPointIndex() == 0)    {
            regular->setEnabled(false);
            incClamp->setEnabled(false);
            incDir->setEnabled(false);
            dataBall->setEnabled(false);
        }
        else if (graph->getSelectedPointIndex() == graph->getSelectedCurve()->getCount()-1)    {
            regular->setEnabled(false);
            outClamp->setEnabled(false);
            outDir->setEnabled(false);
            dataBall->setEnabled(false);
        }
    }
    if (graph->getSelectedCurve()->alphaMax() != alphaMax->value())
        alphaMax->setValue(graph->getSelectedCurve()->alphaMax());
    if (graph->getSelectedCurve()->lambda() != lambda->value())
        lambda->setValue(graph->getSelectedCurve()->lambda());
    if (int(graph->getCurveMethod()) != curveMethod->currentIndex())    {
        curveMethod->setCurrentIndex(int(graph->getCurveMethod()));
        if (curveMethod->currentIndex() == 5 || curveMethod->currentIndex() == 11) //Cond G2 or C2 Cubic
            lambda->setEnabled(true);
        else
            lambda->setEnabled(false);
    }
        
    if (dataBall->isChecked() == true)    {
         QLabel* label;
         label = (QLabel*)nodeSpecific->labelForField(regular);
         label->setText("Marked:");
         label = (QLabel*)nodeSpecific->labelForField(incDirLayout);
         label->setText("Data Ball Radius:");  
         label = (QLabel*)nodeSpecific->labelForField(outDirLayout);
         label->setText("Piece Direction:");  
         label = (QLabel*)nodeSpecific->labelForField(cornerLayout);
         label->setText("Node Direction:"); 
     }
     else if (dataBall->isChecked() == false)    {
         QLabel* label;
         label = (QLabel*)nodeSpecific->labelForField(regular);
         label->setText("Regular:");
         label = (QLabel*)nodeSpecific->labelForField(incDirLayout);
         label->setText("Incoming Direction:");  
         label = (QLabel*)nodeSpecific->labelForField(outDirLayout);
         label->setText("Outgoing Direction:");  
         label = (QLabel*)nodeSpecific->labelForField(cornerLayout);
         label->setText("Corner:");  
     }

    

}

void GraphUI::handleChanges()
{
     if (userChange == false)
         return;
     
     Curve* curve = graph->getSelectedCurve();
     
     if (sender() == closed)    {

         curve->setClosed(closed->isChecked());
     }
     if (sender() == regular)    {
         Curve::Node node = graph->getNode(graph->getSelectedPointIndex());
         if (!node.isDataBall())
             node.setRegular(regular->isChecked());
         else    {
              node.toggleTick();

         }
         curve->setNode(graph->getSelectedPointIndex(), node);
     }
     if (sender() == dataBall)    {
         Curve::Node node = graph->getNode(graph->getSelectedPointIndex());
         node.setDataBall(dataBall->isChecked());
         curve->setNode(graph->getSelectedPointIndex(), node);
     }
          
     if (sender() == incClamp)    {
         Curve::Node node = graph->getNode(graph->getSelectedPointIndex());
         if (incClamp->isChecked())
             node.setIncClamp(node.incDir());
         else
             node.removeIncClamp();
         curve->setNode(graph->getSelectedPointIndex(), node);
             
     }
     if (sender() == outClamp)    {
         Curve::Node node = graph->getNode(graph->getSelectedPointIndex());
         if (outClamp->isChecked())
             node.setOutClamp(node.outDir());
         else
             node.removeOutClamp();
         curve->setNode(graph->getSelectedPointIndex(), node);
     }
     if (sender() == incDir)    {
           Curve::Node node = graph->getNode(graph->getSelectedPointIndex());
           if (node.isDataBall())
               node.setDataBall(true, incDir->value());
           else
               node.setIncClamp(incDir->value());
           curve->setNode(graph->getSelectedPointIndex(), node);
     }
     else if (sender() == outDir)    {
         Curve::Node node = graph->getNode(graph->getSelectedPointIndex());
         node.setOutClamp(outDir->value());
         curve->setNode(graph->getSelectedPointIndex(), node);
     }
     else if (sender() == corner)    {
          Curve::Node node = graph->getNode(graph->getSelectedPointIndex());
          node.setCorner(corner->value());
          curve->setNode(graph->getSelectedPointIndex(), node);
     }
     else if (sender() == posx || sender() == posy)    {
         Curve::Node node = graph->getNode(graph->getSelectedPointIndex());
         node.setPos(posx->value(), posy->value());
         curve->setNode(graph->getSelectedPointIndex(), node);
     } 
     else if (sender() == respectCorner)    {
         Curve::Node node = graph->getNode(graph->getSelectedPointIndex());
         node.respectCorner(respectCorner->isChecked());
         curve->setNode(graph->getSelectedPointIndex(), node);
     } 
     else if (sender() == alphaMax)
         graph->getSelectedCurve()->setAlphaMax(alphaMax->value());
     else if (sender() == lambda)
         graph->getSelectedCurve()->setLambda(lambda->value());
     else if (sender() == curveMethod)    {
             graph->setCurveMethod((GraphArea::CurveMethod)curveMethod->currentIndex());
             if (curveMethod->currentIndex() == 5 || curveMethod->currentIndex() == 11) //Cond G2 or C2 Cubic
                 lambda->setEnabled(true);
             else
                 lambda->setEnabled(false);
     }


    graph->optimizeOnFly();
    updateWidgets();
    graph->update();
}


void GraphUI::setDirection(double dir)
{
     updateWidgets();
}

void GraphUI::setPosition(double x, double y)
{
     updateWidgets();
}

void GraphUI::keyPressEvent(QKeyEvent* event)
{

}


