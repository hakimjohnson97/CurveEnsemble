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


#include "GraphPaintEditor.h"

#define CURVE_INDEX graph->getSelectedCurveIndex()
#define NODE_INDEX graph->getSelectedPointIndex()

GraphPaintEditor::GraphPaintEditor(GraphArea* arg_graph, QWidget *parent) : QFrame(parent)
{    
     graph = arg_graph;
     connect(graph, SIGNAL(selectionChanged()), this, SLOT(updateWidgets()));
     
     applyToAll = new QPushButton("Apply To All Nodes");
     connect(applyToAll, SIGNAL(clicked()), this, SLOT(setAllNodes()));
     
     useNodeDisplay = new QComboBox();
     useNodeDisplay->addItem("Use Recommended Display");
     useNodeDisplay->addItem("Use User Specified Display");
     connect(useNodeDisplay, SIGNAL(currentIndexChanged(int)), this, SLOT(setUseNodeDisplay()));
     
     backgroundColorButton = new QPushButton("Pick Color");
     connect(backgroundColorButton, SIGNAL(clicked()), this, SLOT(pickColor()));
     
     pieceWidth = new QDoubleSpinBox();
     connect(pieceWidth, SIGNAL(valueChanged(double)), this, SLOT(setPieceWidth()));
     pieceColorButton = new QPushButton("Pick Color");
     connect(pieceColorButton, SIGNAL(clicked()), this, SLOT(pickColor()));
     showPiece = new QCheckBox();
     connect(showPiece, SIGNAL(stateChanged(int)), this, SLOT(setShowPiece()));
     
     pointRadius = new QDoubleSpinBox();
     connect(pointRadius, SIGNAL(valueChanged(double)), this, SLOT(setPointRadius()));
     pointColorButton = new QPushButton("Pick Color");
     connect(pointColorButton, SIGNAL(clicked()), this, SLOT(pickColor()));
     pointStyle = new QComboBox();
     pointStyle->addItem("Display Point");
     pointStyle->addItem("Display Point With Arrow");
     pointStyle->addItem("Display Point With Tail");
     pointStyle->addItem("Display Point With Arrow and Tail");
     pointStyle->addItem("Display Arrow");
     pointStyle->addItem("Display Reverse Arrow");
     pointStyle->addItem("No Display");
     pointStyle->addItem("Display Point With Tick");
     connect(pointStyle, SIGNAL(currentIndexChanged(int)), this, SLOT(setPointStyle()));

     QFormLayout *layout = new QFormLayout (this);
     layout->addRow("", applyToAll);
     layout->addRow("", useNodeDisplay);
     layout->addRow("Background Color: ", backgroundColorButton);
     layout->addRow("Piece Width ", pieceWidth);
     layout->addRow("Piece Color ", pieceColorButton);
     layout->addRow("Show Piece ", showPiece);
     layout->addRow("Point Radius ", pointRadius);
     layout->addRow("Point Color: ", pointColorButton);
     layout->addRow("Point Style ", pointStyle);
     setLayout(layout);

}

void GraphPaintEditor::updateWidgets()
{
     NodeDisplay nodeDisplay = graph->getSelectedNodeDisplay();
     if (graph->isNodeDisplayUsed())
         useNodeDisplay->setCurrentIndex(1);
     else
         useNodeDisplay->setCurrentIndex(0);

     backgroundColor = graph->getBackgroundColor();

     pieceWidth->setValue(nodeDisplay.pieceWidth);
     pieceColor = nodeDisplay.pieceColor;
     showPiece->setChecked(nodeDisplay.showPiece);
     
     pointRadius->setValue(nodeDisplay.pointRadius);
     pointColor = nodeDisplay.pointColor;
     pointStyle->setCurrentIndex(nodeDisplay.pointStyle); 
}

void GraphPaintEditor::setAllNodes()
{
    int i;
    for (i = 0; i < graph->getSelectedCurve()->getCount(); i++)    {
        NodeDisplay nodeDisplay = graph->getNodeDisplay(i);
        nodeDisplay.pieceWidth = pieceWidth->value();
        nodeDisplay.showPiece = showPiece->isChecked();
        nodeDisplay.pointRadius = pointRadius->value();
        nodeDisplay.pointStyle = (PointStyle)pointStyle->currentIndex();
        nodeDisplay.pieceColor = pieceColor;
        nodeDisplay.pointColor = pointColor;
        graph->setNodeDisplay(i, nodeDisplay);
    }
    
}

void GraphPaintEditor::setUseNodeDisplay()
{
     if (useNodeDisplay->currentIndex() == 0)
         graph->setUseNodeDisplay(false);
     else
         graph->setUseNodeDisplay(true);
     graph->update();
}

void GraphPaintEditor::setPieceWidth()
{
    NodeDisplay nodeDisplay = graph->getSelectedNodeDisplay();
    nodeDisplay.pieceWidth = pieceWidth->value();
    graph->setSelectedNodeDisplay(nodeDisplay);
}

void GraphPaintEditor::setShowPiece()
{
    NodeDisplay nodeDisplay = graph->getSelectedNodeDisplay();
    nodeDisplay.showPiece = showPiece->isChecked();
    graph->setSelectedNodeDisplay(nodeDisplay);
}

void GraphPaintEditor::setPointRadius()
{
    NodeDisplay nodeDisplay = graph->getSelectedNodeDisplay();
    nodeDisplay.pointRadius = pointRadius->value();
    graph->setSelectedNodeDisplay(nodeDisplay);
}

void GraphPaintEditor::setPointStyle()
{
    NodeDisplay nodeDisplay = graph->getSelectedNodeDisplay();
    nodeDisplay.pointStyle = (PointStyle)pointStyle->currentIndex();
    graph->setSelectedNodeDisplay(nodeDisplay);
}

void GraphPaintEditor::pickColor()
{
     NodeDisplay nodeDisplay = graph->getSelectedNodeDisplay();
     QColor temp;
     if (sender() == backgroundColorButton)     {
         temp = QColorDialog::getColor(backgroundColor, this);
         if (temp.isValid())    {
             backgroundColor = temp;
             graph->setBackgroundColor(backgroundColor);
         }
     }
     if (sender() == pieceColorButton)    {
         temp = QColorDialog::getColor(pieceColor, this);
         if (temp.isValid())    {
             pieceColor = temp;
             nodeDisplay.pieceColor = pieceColor;
         }
     }
     if (sender() == pointColorButton)      {
         temp = QColorDialog::getColor(pointColor, this);
         if (temp.isValid())    {
             pointColor = temp;
             nodeDisplay.pointColor = pointColor;
         }
     }
     graph->setSelectedNodeDisplay(nodeDisplay);
}

void GraphPaintEditor::keyPressEvent(QKeyEvent* event)
{
     if (event->key() == Qt::Key_Left)    {
         graph->setSelection(graph->getSelectedPointIndex()-1);
         graph->update();
     }
     else if (event->key() == Qt::Key_Right)    {
         graph->setSelection(graph->getSelectedPointIndex()+1);
         graph->update();
     }
     else if (event->key() == Qt::Key_F2)    {
         graph->optimize(Curve::Smart);
     }
}

