#include "OptimizeWorker.h"

OptimizeWorker::OptimizeWorker(QObject *parent) :
    QObject(parent)
{
}

void OptimizeWorker::optimizeOnFly()
{
    mutex.lock();
    curve->optimizeToTemp();
    mutex.unlock();
    emit optimizeFinished();
}
