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


#ifndef WINDOW_H
#define WINDOW_H

#include <QMainWindow>
#include <QWidget>
#include <QPushButton>
#include <QScrollArea>
#include <QComboBox>
#include <QSpinBox>
#include <QCheckBox>
#include <QActionGroup>
#include <QCloseEvent>
#include <QSignalMapper>
#include <QPaintEvent>
#include "GraphArea.h"
#include "GraphUI.h"
#include "GraphPaintEditor.h"

/*
  Window is the main window of the GUI and is the parent of all the other objects in the project.
  It has the toolbar with File, View and Tools which each have many options; it allows the user to change
  the window mode to Transform or Zoom; it also has the Optimize button which optimizes the curve.
  Window holds instances of GraphArea, GraphUI and GraphPaintEditor which form the core of editing and displaying the curves

*/

class Window : public QMainWindow
{
  Q_OBJECT
  public:
      Window();
      void createPaintOptions();
      void createImageDialog();
      void paintEvent(QPaintEvent* event);
  protected:
      void closeEvent(QCloseEvent*);
      bool eventFilter(QObject*, QEvent*);
  public slots:
      void addTab();
      void tabChanged(int index);
      void removeTab(int index);
      void pickColor();
      void centerPoint();
      void setGraphMode(QAction* action);
      void setStatusBarMouseCoords(int, int);

      void newProject();
      void saveProject();
      void saveAsProject();
      void openProject();
      void saveToImage();
      void saveToImage2();
      void saveImageReject();
      void loadBackground();
      void transform();
      void viewCurveEditor();
      void viewPaintEditor();
      void updateCurveTab();
      void optimizeCurve();
      void fineOptimizeCurve();
      void endTransformMode();
      void subdivideCurve();
      void deleteSubdivision();
      void anchorSubdivision();
      void onFlyChanged(bool);
      void convertCirclesToDataBalls();
      void createShortestPath();
      void createSmoothestCurve();
      void continueStep();
      void refreshBE();
      void refreshCurveLength();
  private:
      GraphArea* graph;
      GraphUI* graphUI;
      GraphPaintEditor* graphPaintEditor;

      QString saveFile;
      QTabBar* curveTab;
      QPushButton* addCurveTab;
      QDockWidget* curveEditorDock;
      QDockWidget* paintEditorDock;
      QPushButton *optimize, *fineOptimize, *subdivide, *step;
      QComboBox *initDir;
      QCheckBox *onFly;
      QLabel* mouseCoords, *bendingEnergy, *curveLength;

      //Actions
      QActionGroup *modes;
      QAction *normalMode, *transformMode, *zoomMode;

      //Paint Options
      QDialog *paintOptions;
      QComboBox *applyToCurve;
      QComboBox *applyToPoint;

      QCheckBox *closed;
      QPushButton *backgroundColorButton;
      QColor backgroundColor;

      QSpinBox *pieceWidth;
      QPushButton *pieceColorButton;
      QColor pieceColor;
      QCheckBox *showPiece;

      QSpinBox *pointRadius;
      QPushButton *pointColorButton;
      QColor pointColor;
      QComboBox *pointStyle;

      QPushButton *done, *apply;

      //Save Image dialog
      QDialog *imageDialog;
      QLabel *saveImageLabel;
      QPushButton *centerClipRectButton;
      QDoubleSpinBox *scaleWidget, *steph;
      QPushButton *saveImageButton;
      QString imageFormat;

      Curve* dataBallCurve;

      void checkIfCurrentTabIsEmpty(); //Fix bug
};

#endif
