#ifndef SPHERE_H
#define SPHERE_H

#include "Geometry.h"
#include "Point3D.h"

class Sphere : public Geometry
{
private:

public:

    double radius;
    Point3D center;

    Sphere(double radius = 0., Point3D center = Point3D(0., 0., 0.));

    double surfaceArea();
    double volume();
    double diameter(){ return 2.*radius; }

    Point3D nearestIntersect(const Point3D& point);
    std::pair<Point3D, Point3D> lineIntersect(const Point3D& pt1, const Point3D& pt2);

    bool isInside(const Point3D& point);

    Sphere& operator+=(const Vector3D& translationVec);
    Sphere& operator-=(const Vector3D& translationVec);
    Sphere& rotate(const double radians);
    Sphere& scale(const double a);
};

#endif
