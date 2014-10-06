#include "Line3D.h"

Line3D::Line3D(const Point3D& point1,
               const Point3D& point2)
    :
      point1_(point1),
      point2_(point2)
{
    length_ = (point2_ - point1_).mag();
}

double Line3D::length()
{
    return length_;
}
