#ifndef KERNEL_H
#define KERNEL_H

#include "Point3D.h"

class Kernel
{
public:
    enum Type{K6, K8};

    Kernel(Type type, double epsilon);

    void changeType(Type newType);

    double value(const Point3D &center, const Point3D &point);

private:

    void computeA();

    Type type_;
    double epsilon_, A_;
};

#endif
