#include "Interpolation.h"

Matrix Interpolation::M_;
Matrix Interpolation::b_;

double Interpolation::linear(const Point3D* points,
                             const double* values,
                             int nPoints,
                             const Point3D &interpolationPoint)
{
    int i;

    if(nPoints < 4)
        throw "in \"Interpolation::linear\", linear interpolation requires at least 4 points.";

    M_.reallocate(nPoints, 4);
    b_.reallocate(nPoints, 1);

    for(i = 0; i < nPoints; ++i)
    {
        M_(i, 0) = 1.;
        M_(i, 1) = points[i].x;
        M_(i, 2) = points[i].y;
        M_(i, 3) = points[i].z;

        b_(i, 0) = values[i];
    }

    M_.solveLeastSquares(b_);

    return b_(0, 0)
            + b_(1, 0)*interpolationPoint.x
            + b_(2, 0)*interpolationPoint.y
            + b_(3, 0)*interpolationPoint.z;
}

double Interpolation::quadratic(const Point3D* points,
                                const double* values,
                                int nPoints,
                                const Point3D &interpolationPoint)
{
    int i;
    double x, y, z;

    if(nPoints < 10)
        throw "in \"Interpolation::quadratic\", quadratic interpolation requires at least 10 points.";

    M_.reallocate(nPoints, 10);
    b_.reallocate(nPoints, 1);

    for(i = 0; i < nPoints; ++i)
    {
        x = points[i].x;
        y = points[i].y;
        z = points[i].z;

        M_(i, 0) = 1.;
        M_(i, 1) = x;
        M_(i, 2) = y;
        M_(i, 3) = z;
        M_(i, 4) = x*x;
        M_(i, 5) = y*y;
        M_(i, 6) = z*z;
        M_(i, 7) = x*y;
        M_(i, 8) = y*z;
        M_(i, 9) = x*z;

        b_(i, 0) = values[i];
    }

    M_.solveLeastSquares(b_);

    x = interpolationPoint.x;
    y = interpolationPoint.y;
    z = interpolationPoint.z;

    return b_(0, 0)
            + b_(1, 0)*x
            + b_(2, 0)*y
            + b_(3, 0)*z
            + b_(4, 0)*x*x
            + b_(5, 0)*y*y
            + b_(6, 0)*z*z
            + b_(7, 0)*x*y
            + b_(8, 0)*y*z
            + b_(9, 0)*x*z;
}
