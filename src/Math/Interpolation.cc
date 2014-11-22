#include "Interpolation.h"

double Interpolation::linear(const Point3D& point,
                             const double* values, const Point3D* points,
                             int nPoints)
{

    int i;
    double distanceSum(0.), interpolatedValue(0.);

    for(i = 0; i < nPoints; ++i)
    {

        distanceSum += (points[i] - point).mag();

    }

    for(i = 0; i < nPoints; ++i)
    {

        interpolatedValue += values[i]*(points[i] - point).mag()/distanceSum;

    }

    return interpolatedValue;

}
