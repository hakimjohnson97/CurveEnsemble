#-------------------------------------------------
#
# Project created by QtCreator 2013-03-16T20:05:02
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Curve-Ensemble
TEMPLATE = app


SOURCES += main.cpp\
    GraphUI.cpp \
    GraphPaintEditor.cpp \
    GraphArea.cpp \
    Window.cpp \
    RectWidget.cpp \
    Specific_Curves/Curve.cpp \
    Globals.cpp \
    Specific_Curves/ScottSpline1.cpp \
    Specific_Curves/ResElasticSpline.cpp \
    Specific_Curves/ElasticSpline.cpp \
    Specific_Curves/CubicSpline.cpp \
    Specific_Curves/CubicQuasiElasticSpline.cpp \
    Specific_Curves/Curve_Libs/common_curve_lib.cpp \
    Specific_Curves/CondG2Cubic.cpp \
    Specific_Curves/JZSpline.cpp \
    OptimizeWorker.cpp \
    Specific_Curves/C2ParCubic.cpp \
    Specific_Curves/Lines.cpp



HEADERS  += \
    GraphUI.h \
    GraphPaintEditor.h \
    GraphArea.h \
    Window.h \
    S_curve_lib.h \
    RectWidget.h \
    Specific_Curves/Curve.h \
    Specific_Curves/ScottSpline1.h \
    Specific_Curves/Curve_Libs/scott_spline1_lib.h \
    Specific_Curves/Curve_Libs/restricted_elastic_spline_lib.h \
    Specific_Curves/ResElasticSpline.h \
    Specific_Curves/ElasticSpline.h \
    Specific_Curves/Curve_Libs/elastic_spline_lib.h \
    Specific_Curves/CubicSpline.h \
    Specific_Curves/CubicQuasiElasticSpline.h \
    Specific_Curves/Curve_Libs/cubic_spline_lib.h \
    Specific_Curves/Curve_Libs/cubic_quasi_elastic_spline_lib.h \
    Specific_Curves/Curve_Libs/common_curve_lib.h \
    Specific_Curves/CondG2Cubic.h \
    Specific_Curves/Curve_Libs/cond_G2_cubic_lib.h \
    Specific_Curves/JZSpline.h \
    Specific_Curves/Curve_Libs/JZ_spline_lib.h \
    OptimizeWorker.h \
    Specific_Curves/Curve_Libs/c2_par_cubic_lib.h \
    Specific_Curves/C2ParCubic.h \
    Specific_Curves/Lines.h \
    Specific_Curves/Curve_Libs/lines_lib.h 

CONFIG += release

FORMS    +=
