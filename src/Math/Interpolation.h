#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>

#include "Matrix.h"
#include "Point3D.h"

class Interpolation
{
private:

    static Matrix M_, b_;

public:

    static double linear(const Point3D* points,
                         const double* values,
                         int nPoints,
                         const Point3D& interpolationPoint);

    static double trilinear(const Point3D* points,
                            const double* values,
                            int nPoints,
                            const Point3D& interpolationPoint);

    static double quadratic(const Point3D* points,
                            const double* values,
                            int nPoints,
                            const Point3D& interpolationPoint);

    static Matrix computeTrilinearCoeffs(const Point3D* points,
                                         int nPoints,
                                         const Point3D& interpolationPoint);
};

#endif
