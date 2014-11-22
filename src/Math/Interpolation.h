#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "Point3D.h"

class Interpolation
{

private:

public:

    static double linear(const Point3D& point,
                         const double* values, const Point3D* points,
                         int nPoints);

};

#endif
