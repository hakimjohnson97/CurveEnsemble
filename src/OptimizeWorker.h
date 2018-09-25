#ifndef OPTIMIZEWORKER_H
#define OPTIMIZEWORKER_H

#include <QObject>
#include <QMutex>
#include "Specific_Curves/Curve.h"

extern QMutex mutex;

/*
  This class is used in achieving optimizing on the fly. Optimizing a curve can take a long time,
  so we optimize the curve on a separate thread to avoid the program from lagging. This class is used
  in conjunction with the QThread class.
*/
class OptimizeWorker : public QObject
{
    Q_OBJECT
public:
    explicit OptimizeWorker(QObject *parent = 0);
    void setCurve(Curve* curve_arg) {curve = curve_arg;}
signals:
    void optimizeFinished();
public slots:
    void optimizeOnFly();
private:
    Curve* curve;
};

#endif // OPTIMIZEWORKER_H
