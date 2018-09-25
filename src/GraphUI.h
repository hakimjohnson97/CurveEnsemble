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


#ifndef GRAPHUI_H
#define GRAPHUI_H

#include <QWidget>
#include <QFrame>
#include <QLineEdit>
#include <QCheckBox>
#include <QTabBar>
#include <QSlider>
#include <QPushButton>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QButtonGroup>
#include <QRadioButton>
#include <QComboBox>
#include <QFormLayout>
#include <QHBoxLayout>
#include "GraphArea.h"

/*
  GraphUI is a separate window that shows non-display oriented properties of the currently
  selected node (e.g coordinates) and allows the user to change them.
  It also has options that affect the curve as a whole such as the curve's closed properties.
*/
class GraphUI : public QFrame
{
  Q_OBJECT
  public:
      GraphUI (GraphArea *graph, QWidget *parent = 0);
  signals:
      void directionChanged(double);
      void positionChanged(double, double);
      void stateChanged(int);
      void cornerChanged(double);
      void newTab();
      void tabChanged(int index);
      void tabRemoved(int index);
      void scaleChanged(double);
      void findPoint(int, int);
  protected:
      void keyPressEvent(QKeyEvent*);

  public slots:
      void updateWidgets();
      void setDirection(double);
      void setPosition(double, double);
      void handleChanges();

  private:
      GraphArea *graph;
      QCheckBox *closed;
      QDoubleSpinBox *posx, *posy;
      QCheckBox *regular, *dataBall;
      QPushButton *incClamp, *outClamp;
      QDoubleSpinBox *incDir, *outDir, *corner;
      QCheckBox *respectCorner;

      QCheckBox *clamped, *free;
      QDoubleSpinBox *alphaMax, *lambda;
      QComboBox *curveMethod;

      QFormLayout *nodeSpecific;
      QHBoxLayout *incDirLayout, *outDirLayout, *cornerLayout;




      bool userChange;
};

#endif
