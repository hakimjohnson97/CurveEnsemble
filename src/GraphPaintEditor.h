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


#ifndef GRAPHPAINTEDITOR_H
#define GRAPHPAINTEDITOR_H

#include <QWidget>
#include <QFrame>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QPushButton>
#include <QFormLayout>
#include <QColorDialog>
#include <QKeyEvent>
#include "GraphArea.h"



/*
  GraphPaintEditor is a separate window that shows the display-oriented properties of the
  currently selected node and piece and allows the user to change them.
*/
class GraphPaintEditor : public QFrame
{
  Q_OBJECT
  public:
      GraphPaintEditor (GraphArea* arg_graph, QWidget *parent = 0);
  public slots:
      void updateWidgets();

      void setAllNodes();

      void setUseNodeDisplay();
      void setPieceWidth();
      void setShowPiece();
      void setPointRadius();
      void setPointStyle();

      void pickColor();
  protected:
      void keyPressEvent(QKeyEvent*);

  private:
      GraphArea *graph;

      QPushButton *applyToAll;

      QComboBox *useNodeDisplay;
      QPushButton *backgroundColorButton;
      QColor backgroundColor;

      QDoubleSpinBox *pieceWidth;
      QPushButton *pieceColorButton;
      QColor pieceColor;

      QCheckBox *showPiece;

      QDoubleSpinBox *pointRadius;
      QPushButton *pointColorButton;
      QColor pointColor;
      QComboBox *pointStyle;
};

#endif
