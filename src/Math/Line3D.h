#ifndef LINE3D_H
#define LINE3D_H

#include "Point3D.h"

class Line3D
{

private:

    double length_;

    Point3D point1_, point2_;

public:

    Line3D(const Point3D& point1 = Point3D(),
           const Point3D& point2 = Point3D());

    double length();
    Line3D& rotate(double radians);
    Line3D& rotate(double radians, const Point3D& origin);
};

#endif
