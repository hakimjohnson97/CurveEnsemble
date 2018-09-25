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
#include <QVBoxLayout>
#include <QTime>
#include <cmath>
#include <QTabBar>
#include <QToolBar>
#include <QMainWindow>
#include <QStatusBar>
#include <QMenuBar>
#include <QScrollArea>
#include <fstream>
#include <QDockWidget>
#include <QDialog>
#include <QComboBox>
#include <QSpinBox>
#include <QFormLayout>
#include <QImageWriter>
#include <QFileDialog>
#include <QColorDialog>
#include <QInputDialog>
#include <QCheckBox>
#include <QMouseEvent>
#include <QDebug>
#include <stdio.h>
#include <time.h>
#include "Window.h"
using namespace std;

fstream file6;

#define X 0
#define Y 1
#define DIR 2
#define CORNER 3
#define STATE 4
#define BE 5

#define TITLE "Curve Ensemble - "

QString nameOfFile(QString loc)
{
     string str = loc.toStdString();
     int p;
     p = str.rfind("/");
     return QString(str.substr(p+1, str.length() - p -1).c_str());
}

Window::Window()
{ 
                
    setContextMenuPolicy(Qt::NoContextMenu);            
    QWidget *centralWidget = new QWidget(this);
    graph = new GraphArea(); 
    graphUI = new GraphUI(graph); 
    graphPaintEditor = new GraphPaintEditor(graph); 
    
    setWindowTitle(TITLE + QString("Untitled"));
    
    
    optimize = new QPushButton ("Optimize");
    QLabel *initDirLabel = new QLabel("with starting directions chosen ");
    initDir = new QComboBox();
    onFly = new QCheckBox("on the fly");
    onFly->setChecked(true);
    initDir->addItem("Smartly");
    initDir->addItem("Randomly");
    initDir->addItem("AsIs");

    QHBoxLayout *optimizeLayout = new QHBoxLayout;
    optimizeLayout->addWidget(optimize);
    optimizeLayout->addWidget(initDirLabel);
    optimizeLayout->addWidget(initDir);
    optimizeLayout->addWidget(onFly);
    

    fineOptimize = new QPushButton ("Fine Optimize", this);
    connect(optimize, SIGNAL(clicked()), this, SLOT(optimizeCurve()));
    connect(onFly, SIGNAL(toggled(bool)), this, SLOT(onFlyChanged(bool)));
    connect(fineOptimize, SIGNAL(clicked()), this, SLOT(fineOptimizeCurve()));


    step = new QPushButton ("Run Command", this);
    connect(step, SIGNAL(clicked()), this, SLOT(continueStep()));

    steph = new QDoubleSpinBox(this);
    steph->setRange(-100000000000, 10000000000);
    steph->setDecimals(6);
    steph->setValue(25);

    QHBoxLayout *fineOptimizeLayout = new QHBoxLayout;
    fineOptimizeLayout->addWidget(fineOptimize);
    fineOptimizeLayout->addWidget(step);
    fineOptimizeLayout->addWidget(steph);
    
    connect(graphUI, SIGNAL(findPoint(int, int)), this, SLOT(showArea(int, int)));
    
    curveTab = new QTabBar;
    curveTab->addTab("1");
    curveTab->setTabsClosable(true);
    connect(curveTab, SIGNAL(currentChanged(int)), this, SLOT(tabChanged(int)));
    connect(curveTab, SIGNAL(tabCloseRequested(int)), this, SLOT(removeTab(int)));
    
    addCurveTab = new QPushButton("+");
    addCurveTab->setMaximumSize(30, 30);
    QHBoxLayout* curveTabLayout = new QHBoxLayout;
    curveTabLayout->addWidget(curveTab);
    curveTabLayout->addWidget(addCurveTab);
    connect(addCurveTab, SIGNAL(clicked()), this, SLOT(addTab()));
    

    setDockOptions(QMainWindow::AnimatedDocks | QMainWindow::ForceTabbedDocks);
     curveEditorDock = new QDockWidget("Curve Editor", this);
     curveEditorDock->setAllowedAreas(Qt::BottomDockWidgetArea);
     curveEditorDock->setWidget(graphUI);
     addDockWidget(Qt::NoDockWidgetArea, curveEditorDock);
     curveEditorDock->setFloating(true);
     
     paintEditorDock = new QDockWidget("Paint Editor", this);
     paintEditorDock->setAllowedAreas(Qt::BottomDockWidgetArea);
     paintEditorDock->setWidget(graphPaintEditor);
     addDockWidget(Qt::NoDockWidgetArea, paintEditorDock); 
     paintEditorDock->setVisible(false);
     paintEditorDock->setFloating(true);
    
    
    
    
    QVBoxLayout *layout = new QVBoxLayout;
    layout->addLayout(curveTabLayout);
    layout->addWidget(graph);
    layout->addLayout(optimizeLayout);
    layout->addLayout(fineOptimizeLayout);
    centralWidget->setLayout(layout);
    setCentralWidget(centralWidget);
    
    QMenuBar *menubar = new QMenuBar();
    
    QMenu *file = menubar->addMenu("File");
    QAction *clear = file->addAction("New");
    connect(clear, SIGNAL(triggered()), this, SLOT(newProject()));
    QAction *open = file->addAction("Open");
    connect((QObject*)open, SIGNAL(triggered()), this, SLOT(openProject()));
    QAction *save = file->addAction("Save");
    connect((QObject*)save, SIGNAL(triggered()), this, SLOT(saveProject()));
    QAction *saveAs = file->addAction("SaveAs");
    connect((QObject*)saveAs, SIGNAL(triggered()), this, SLOT(saveAsProject()));
    QAction *loadBackground = file->addAction("Load Background");
    connect((QObject*)loadBackground, SIGNAL(triggered()), this, SLOT(loadBackground()));
    
    QMenu *exportAct = file->addMenu("Export");
    
    QMenu *exportAsImage = exportAct->addMenu("as Image");
    QList<QAction *> imageFormatActs;
    foreach (QByteArray format, QImageWriter::supportedImageFormats()) {
         QString text = tr("%1...").arg(QString(format).toUpper());

         QAction *action = new QAction(text, this);
         action->setData(format);
         connect(action, SIGNAL(triggered()), this, SLOT(saveToImage()));
         imageFormatActs.append(action);
     }
     foreach (QAction *action, imageFormatActs)
         exportAsImage->addAction(action);
  
   QAction *exportAsText = exportAct->addAction("as Text");
   connect(exportAsText, SIGNAL(triggered()), graph, SLOT(saveAsText()));
   QAction *exportToClipboard = exportAct->addAction("to Clipboard");
   connect(exportToClipboard, SIGNAL(triggered()), graph, SLOT(saveToClipboard()));
   
   QMenu *importMenu = file->addMenu("Import");
   
   QAction *importAsText = importMenu->addAction("from Text");
   connect(importAsText, SIGNAL(triggered()), graph, SLOT(loadFromText()));
   QAction *importToClipboard = importMenu->addAction("from Clipboard");
   connect(importToClipboard, SIGNAL(triggered()), graph, SLOT(loadFromClipboard()));
   
   QAction *quit = file->addAction("Quit");
   connect(quit, SIGNAL(triggered()), this, SLOT(close()));
  
   QMenu *view = menubar->addMenu("View");
   QAction *viewCurveEditor = view->addAction("Curve Editor");

   connect((QObject*)viewCurveEditor, SIGNAL(triggered()), this, SLOT(viewCurveEditor()));
   QAction *viewPaintEditor = view->addAction("Paint Editor");

   connect((QObject*)viewPaintEditor, SIGNAL(triggered()), this, SLOT(viewPaintEditor()));
   
    
    QMenu* tools = menubar->addMenu("Tools");
    QAction* convertToDataBall = tools->addAction("Convert Circles To Data Balls");
    connect((QObject*)convertToDataBall, SIGNAL(triggered()), this, SLOT(convertCirclesToDataBalls()));
    QAction* createShortestPathAction = tools->addAction("Shortest Path through Data Balls");
    connect((QObject*)createShortestPathAction, SIGNAL(triggered()), this, SLOT(createShortestPath()));
    QAction* createSmoothestCurveAction = tools->addAction("Smoothest Curve through Data Balls");
    connect((QObject*)createSmoothestCurveAction, SIGNAL(triggered()), this, SLOT(createSmoothestCurve()));
    tools->addSeparator();
    QAction* subdivideAction = tools->addAction("Subdivide Curve");
    connect((QObject*)subdivideAction, SIGNAL(triggered()), this, SLOT(subdivideCurve()));
    QAction* deleteSubdivisionAction = tools->addAction("Delete Subdivision Nodes");
    connect((QObject*)deleteSubdivisionAction, SIGNAL(triggered()), this, SLOT(deleteSubdivision()));
    QAction* anchorSubdivisionAction = tools->addAction("Anchor Subdivision Nodes");
    connect((QObject*)anchorSubdivisionAction, SIGNAL(triggered()), this, SLOT(anchorSubdivision()));
    
    setMenuBar(menubar);
    
    QToolBar* toolbar = addToolBar("main toolbar");
    normalMode = toolbar->addAction("Normal");
    normalMode->setCheckable(true);
    normalMode->setChecked(true);
    transformMode = toolbar->addAction("Transform");
    transformMode->setCheckable(true);

    zoomMode = toolbar->addAction("Zoom");
    zoomMode->setCheckable(true);
    modes = new QActionGroup(this);
    modes->addAction(normalMode);
    modes->addAction(transformMode);
    modes->addAction(zoomMode);
    connect(modes, SIGNAL(triggered(QAction*)), this, SLOT(setGraphMode(QAction*)));
    toolbar->addSeparator();
    
    QAction* centerPoint = toolbar->addAction("Center Point");
    connect((QObject*)centerPoint, SIGNAL(triggered()), this, SLOT(centerPoint()));
    
    connect(graph, SIGNAL(cleared()), this, SLOT(updateCurveTab()));  
 
   
    
    mouseCoords = new QLabel(" 0, 0");
    statusBar()->addWidget(mouseCoords);
    connect(graph, SIGNAL(mouseCoordsChanged(int, int)), this, SLOT(setStatusBarMouseCoords(int, int)));
    bendingEnergy = new QLabel("Energy: 0");
    statusBar()->addWidget(bendingEnergy);
    curveLength = new QLabel("Length: 0");
    statusBar()->addWidget(curveLength);
    
    createPaintOptions();
    createImageDialog();

    connect(graph, SIGNAL(transformWindowClosed()), this, SLOT(endTransformMode()));
    
    qApp->installEventFilter(this);
    
    
}

void Window::createImageDialog()
{
    imageDialog = new QDialog(this);
    connect(imageDialog, SIGNAL(rejected()), this, SLOT(saveImageReject())); 
    
    saveImageLabel = new QLabel("Change the bounding box shown on the graph");
    centerClipRectButton = new QPushButton("Center");
    connect(centerClipRectButton, SIGNAL(clicked()), graph, SLOT(centerClipRect()));
    scaleWidget = new QDoubleSpinBox();
    saveImageButton = new QPushButton("Save");
    connect(saveImageButton, SIGNAL(clicked()), this, SLOT(saveToImage2()));
    QFormLayout *layout = new QFormLayout();
    layout->addRow("", saveImageLabel);
    layout->addRow("Scale Factor:", scaleWidget);
    layout->addRow("Bounding Box", centerClipRectButton);
    layout->addRow("", saveImageButton);
    imageDialog->setLayout(layout);
}

void Window::createPaintOptions()
{
     paintOptions = new QDialog(this);
     applyToCurve = new QComboBox();
     applyToPoint = new QComboBox();
     
     backgroundColorButton = new QPushButton("Pick Color");
     connect(backgroundColorButton, SIGNAL(clicked()), this, SLOT(pickColor()));
     closed = new QCheckBox();
     
     pieceWidth = new QSpinBox();
     pieceColorButton = new QPushButton("Pick Color");
     connect(pieceColorButton, SIGNAL(clicked()), this, SLOT(pickColor()));
     showPiece = new QCheckBox();
     
     pointRadius = new QSpinBox();
     pointColorButton = new QPushButton("Pick Color");
     connect(pointColorButton, SIGNAL(clicked()), this, SLOT(pickColor()));
     pointStyle = new QComboBox();
     pointStyle->addItem("Display Point");
     pointStyle->addItem("Display Arrow");
     pointStyle->addItem("Display Point With Arrow");
     pointStyle->addItem("No Display");
     pointStyle->addItem("Display Reverse Arrow");
             
     done = new QPushButton("Done");
     connect(done, SIGNAL(clicked()), paintOptions, SLOT(accept()));
     apply = new QPushButton("Apply");


     QFormLayout *layout = new QFormLayout (paintOptions);
     layout->addRow("Apply to Curve: ", applyToCurve);
     layout->addRow("Apply to Point: ", applyToPoint);
     layout->addRow("Closed: ", closed);
     layout->addRow("Background Color: ", backgroundColorButton);
     layout->addRow("Piece Width ", pieceWidth);
     layout->addRow("Piece Color ", pieceColorButton);
     layout->addRow("Show Piece ", showPiece);
     layout->addRow("Point Radius ", pointRadius);
     layout->addRow("Point Color: ", pointColorButton);
     layout->addRow("Point Style ", pointStyle);
     layout->addRow("", done);
     paintOptions->setLayout(layout);
     
     
}

void Window::setStatusBarMouseCoords(int x, int y)
{
     mouseCoords->setText(" " + QString::number(x) + " , " + QString::number(y) );
}

void Window::checkIfCurrentTabIsEmpty()
{
    if (graph->getSelectedCurve()->getCount() == 0)  {//Fix bug regarding p-selected = -1
        curveTab->blockSignals(true);
        graph->removeCurve(graph->getSelectedCurveIndex());
        curveTab->removeTab(graph->getSelectedCurveIndex());
        curveTab->blockSignals(false);
     }
}
void Window::addTab()
{
     checkIfCurrentTabIsEmpty();
     curveTab->addTab(QString::number(curveTab->count()+1));
     curveTab->blockSignals(true);
     curveTab->setCurrentIndex(curveTab->count()-1);
     curveTab->blockSignals(false);
     graph->addCurve();
     graphUI->updateWidgets();
}

void Window::tabChanged(int index)
{
     checkIfCurrentTabIsEmpty();
     graph->changeCurve(index);
     graphUI->updateWidgets();
     graphPaintEditor->updateWidgets();
}

void Window::removeTab(int index)
{
    if (index != graph->getSelectedCurveIndex())
        checkIfCurrentTabIsEmpty();
    if (curveTab->count() == 1)
        graph->clear();
    else    {
        graph->removeCurve(index);
        curveTab->removeTab(index);
        curveTab->setCurrentIndex(index-1);
    }
}

void Window::updateCurveTab()
{
     int i;
     curveTab->blockSignals(true);
    while (curveTab->count() > 0)
         curveTab->removeTab(0);
     for (i = 0; i < graph->getCurveCount(); i++)
         curveTab->addTab(QString::number(curveTab->count()+1));
         curveTab->blockSignals(false);
}

void Window::setGraphMode(QAction* action)
{
    if (action == normalMode)
        graph->setGraphMode(GraphArea::Normal);
    else if (action == transformMode)
        graph->setGraphMode(GraphArea::Transform);
    else if (action == zoomMode)
        graph->setGraphMode(GraphArea::Zoom);
}

void Window::transform()
{
     QDialog transformDialog;
     QDoubleSpinBox c1x, c1y;
     QDoubleSpinBox c2x, c2y;
     QCheckBox reflection;
     c1x.setRange(-1000, 1000);
     c1y.setRange(-1000, 1000);
     c2x.setRange(-1000, 1000);
     c2y.setRange(-1000, 1000);
     c1x.setValue(1);
     c1y.setValue(0);
     
     QHBoxLayout c1Layout, c2Layout;
     c1Layout.addWidget(&c1x);
     c1Layout.addWidget(&c1y);
     c2Layout.addWidget(&c2x);
     c2Layout.addWidget(&c2y);
     
     QPushButton done("Done");

     connect(&done, SIGNAL(clicked()), &transformDialog, SLOT(accept()));

     
     QFormLayout *layout = new QFormLayout();
     layout->addRow("C1:", &c1Layout);
     layout->addRow("C2:", &c2Layout);
     layout->addRow("With Reflection:", &reflection);
     layout->addRow("", &done);

     transformDialog.setLayout(layout); 
     
     if (transformDialog.exec())   {
         graph->getSelectedCurve()->transform(c1x.value(), c1y.value(), c2x.value(), c2y.value(), reflection.isChecked());
         if (reflection.isChecked())
             graph->optimizeOnFly();
     }
     
}

void Window::centerPoint()
{ 
     graph->centerNode(graph->getSelectedPointIndex());
}

void Window::pickColor()
{
     if (sender() == backgroundColorButton)
         backgroundColor = QColorDialog::getColor(backgroundColor, this);
     if (sender() == pieceColorButton)
         pieceColor = QColorDialog::getColor(pieceColor, this);
     if (sender() == pointColorButton)
         pointColor = QColorDialog::getColor(pointColor, this);
     
}

void Window::loadPaintOptions()
{ 

}

void Window::savePaintOptions(int curveIndex, int pointIndex)
{ 

}

void Window::setPaintOptions()
{

} 

void Window::newProject()
{
    setWindowTitle(TITLE + QString("Untitled"));
    graph->clear();
    saveFile = QString();
    graphUI->updateWidgets();
}

void Window::saveProject()
{
     if (saveFile == NULL)    {
         saveFile = QFileDialog::getSaveFileName(this, "Save Project", "", "Curve Project (*.cvp);;All Files(*)");
         if (saveFile == NULL)
             return;
     }
     setWindowTitle(TITLE + nameOfFile(saveFile));
     graph->save(saveFile);
}

void Window::saveAsProject()
{
     QString temp = QFileDialog::getSaveFileName(this, "Save Project", "", "Curve Project (*.cvp);;All Files(*)");
     if (temp == NULL)
         return;
     else
         saveFile = temp;
     setWindowTitle(TITLE + nameOfFile(saveFile));
     graph->save(saveFile);
}

void Window::openProject()
{
     QString fileName = QFileDialog::getOpenFileName(this, "Open Project", "", "Curve Project (*.cvp);;All Files(*)");
     if (fileName == NULL)
         return;
     saveFile = fileName;
     setWindowTitle(TITLE + nameOfFile(saveFile));
     graph->open(saveFile);
     
    updateCurveTab();
    graphUI->updateWidgets();
    graphPaintEditor->updateWidgets(); 
}

void Window::saveToImage()
{
    QAction *action = qobject_cast<QAction *>(sender());
    imageFormat = action->data().toString();
    graph->setGraphMode(GraphArea::ClipImageRect);
    scaleWidget->setValue(1);
    imageDialog->show();
    
}

void Window::saveToImage2()
{

    double scale = scaleWidget->value();
     

    QString fileName = QFileDialog::getSaveFileName(this, "Export As " + imageFormat, "", "(*." + imageFormat + ")");
    if (fileName != NULL)
        imageDialog->accept(); 
    else
        return;
    
    QRect clipRect = graph->getClipRect();

    clipRect.setCoords(clipRect.x()*scale, clipRect.y()*scale,
     clipRect.x()*scale+clipRect.width()*scale, clipRect.y()*scale+clipRect.height()*scale);
    QImage image(clipRect.size(), QImage::Format_RGB32);
    image.fill(graph->getBackgroundColor());
    QPainter painter(&image);
    painter.translate(-clipRect.x(), -clipRect.y());
    graph->setScale(graph->getScale()*scale);
    
    graph->setGraphMode(GraphArea::Normal);
    
    graph->drawCurves(painter, false);
    graph->setScale(graph->getScale()/scale);
    image.save(fileName.toStdString().c_str(), imageFormat.toStdString().c_str());
}

void Window::saveImageReject()
{
     graph->setGraphMode(GraphArea::Normal);
}

void Window::loadBackground()
{
    graph->loadBackgroundImage();
}

void Window::viewCurveEditor()
{
         curveEditorDock->setVisible(true);
}

void Window::viewPaintEditor()
{
         paintEditorDock->setVisible(true);
}

void Window::optimizeCurve()
{
    if (!(graph->getCurveMethod() == GraphArea::MethodDataBalls))
        graph->optimize((Curve::InitDir)initDir->currentIndex());

}

void Window::fineOptimizeCurve()
{
     graph->fineOptimize();
}

void Window::closeEvent(QCloseEvent* event)
{
    QDialog maybeSaveDialog;
    QLabel maybeSaveLabel("Do you want to save before exiting?");
    QPushButton maybeSaveYes("Yes");

    QPushButton maybeSaveNo("No");
    QPushButton maybeSaveCancel("Cancel");
    
    QSignalMapper signalMapper;
    connect(&maybeSaveNo, SIGNAL(clicked()), &signalMapper, SLOT(map()));
    signalMapper.setMapping(&maybeSaveNo, 2);
    connect(&maybeSaveYes, SIGNAL(clicked()), &signalMapper, SLOT(map()));
    signalMapper.setMapping(&maybeSaveYes, 1);
    connect(&maybeSaveCancel, SIGNAL(clicked()), &signalMapper, SLOT(map()));
    signalMapper.setMapping(&maybeSaveCancel, 0);
    connect(&signalMapper, SIGNAL(mapped(int)), &maybeSaveDialog, SLOT(done(int)));
    
    
    QHBoxLayout *inputLayout = new QHBoxLayout;
    inputLayout->addWidget(&maybeSaveYes);
    inputLayout->addWidget(&maybeSaveNo);
    inputLayout->addWidget(&maybeSaveCancel);
    
    QVBoxLayout *maybeSaveLayout = new QVBoxLayout();
    maybeSaveLayout->addWidget(&maybeSaveLabel);
    maybeSaveLayout->addLayout(inputLayout);
    
    maybeSaveDialog.setLayout(maybeSaveLayout);
    int done = maybeSaveDialog.exec();
    if (done == 2)
        event->accept(); 
    else if (done == 1)
        saveProject();
    else if (done == 0)
        event->ignore();
}

bool Window::eventFilter(QObject *object, QEvent *event)
{
        
     if (event->type() == QEvent::KeyPress) {
         QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
         if (keyEvent->key() == Qt::Key_Up || keyEvent->key() == Qt::Key_Down) {
             graph->keyPress(keyEvent);
             return true;
         }
         else
             return false;
     }  
     return false;
} 

void Window::endTransformMode()
{
    normalMode->setChecked(true);
}

void Window::subdivideCurve()
{
    graph->getSelectedCurve()->subdivide();
    graph->update();
    graphUI->updateWidgets();
}

void Window::deleteSubdivision()
{
    graph->getSelectedCurve()->endSubdivide();
    graph->update();
    graphUI->updateWidgets();
}

void Window::anchorSubdivision()
{
    graph->anchorSubdivision();
    graph->update();
    graphUI->updateWidgets();
}

void Window::onFlyChanged(bool fly)
{
    graph->setOnFlyOptimization(fly);
}

void Window::convertCirclesToDataBalls()
{
    graph->setCurveMethod(GraphArea::MethodDataBalls);
    graph->getSelectedCurve()->convertCirclesToDataBalls();

    graph->setSelection(0);
    graphUI->updateWidgets();
    graph->update();


}
void Window::createShortestPath()
{

      Curve* curve1 = graph->getSelectedCurve();
      Lines curve2;

      checkIfCurrentTabIsEmpty();
      curveTab->addTab(QString::number(curveTab->count()+1));
      curveTab->blockSignals(true);
      curveTab->setCurrentIndex(curveTab->count()-1);
      curveTab->blockSignals(false);

      double r = curve2.createDataBallCurve(curve1);
      steph->blockSignals(true);
      steph->setValue(r);
      steph->blockSignals(false);

      onFly->setChecked(false);
      onFlyChanged(false);

      graph->addCurve(GraphArea::MethodLines, (Curve*)&curve2);

      graph->update();
      graphUI->updateWidgets();

 
}

void Window::createSmoothestCurve()
{
      Curve* curve1 = graph->getSelectedCurve();
      ResElasticSpline curve2;

      checkIfCurrentTabIsEmpty();
      curveTab->addTab(QString::number(curveTab->count()+1));
      curveTab->blockSignals(true);
      curveTab->setCurrentIndex(curveTab->count()-1);
      curveTab->blockSignals(false);

      double r = curve2.createDataBallCurve(curve1);
      steph->blockSignals(true);
      steph->setValue(r);
      steph->blockSignals(false);

      onFly->setChecked(false);
      onFlyChanged(false);

      graph->addCurve(GraphArea::MethodResElasticSpline, (Curve*)&curve2);

      graph->update();
      graphUI->updateWidgets();
}

void Window::continueStep()
{
    double c = 0;

    QString debugText = QInputDialog::getText(NULL, "Run Command", "Enter your command (an integer): ");
    if (debugText == "")
        return;
    c = double(debugText.toInt());

    double r = graph->getSelectedCurve()->runCommand(c, steph->value());
    steph->blockSignals(true);
    steph->setValue(r);
    steph->blockSignals(false);

    graph->update();
    graphUI->updateWidgets();
    update();
}

void Window::paintEvent(QPaintEvent* event)
{
     refreshBE();
     refreshCurveLength();
     QMainWindow::paintEvent(event);
}

void Window::refreshBE()
{
    bendingEnergy->setText("Energy: " + QString::number(graph->getSelectedCurve()->getCurveBE()));
}

void Window::refreshCurveLength()
{
    curveLength->setText("Length: " + QString::number(graph->getSelectedCurve()->getCurveLength()));
}






/* Old method
void Window::createSmoothestCurve()
{
      Curve* curve1 = graph->getSelectedCurve();
      Curve curve2;

      dataBallCurve = graph->getSelectedCurve();

      checkIfCurrentTabIsEmpty();
      curveTab->addTab(QString::number(curveTab->count()+1));
      curveTab->blockSignals(true);
      curveTab->setCurrentIndex(curveTab->count()-1);
      curveTab->blockSignals(false);

      curve1->createSmoothestCurve(&curve2);

      graph->addCurve(GraphArea::MethodResElasticSpline, curve2);
      graph->getSelectedCurve()->optimize();
      graph->getSelectedCurve()->fineOptimize();
      graph->update();
      graphUI->updateWidgets();
}*/


/* Daddys' alterations
void Window::continueStep()
{
    Curve* curve1 = graph->getSelectedCurve();
    Curve* curve2 = graph->getCurve(1);
    double c = 0;

    QString debugText = QInputDialog::getText(NULL, "Input debug text", "Enter your command (an integer): ");

    c = double(debugText.toInt());

    curve1->createSmoothestCurveStep(curve2, c);

    graph->update();
    graphUI->updateWidgets();
    update();
}*/

/* Old method
void Window::continueStep()
{
    Curve* curve1 = graph->getSelectedCurve();
    Curve* curve2 = dataBallCurve;
    
    int num = 1;

    QString debugText = QInputDialog::getText(NULL, "Input debug text", "How many times would you like to run this?");

    num = debugText.toInt();
    qDebug() << "Running " << num << " times.";

    int t = clock();
    for (int i = 0; i < num; i++)    {
        curve1->createSmoothestCurveStep(dataBallCurve, steph->value());
        graph->update();
        graphUI->updateWidgets();
        update();
    }
    int timeTaken = clock() - t;
    qDebug() << "Finished. Took " << (double(timeTaken)/1000000.0) << " seconds";


}*/
