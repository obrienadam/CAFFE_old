#include <math.h>

#include "Sphere.h"

Sphere::Sphere(double radius, Point3D center)
    :
      radius(radius),
      center(center)
{

}

double Sphere::surfaceArea()
{
    return 4.*M_PI*radius*radius;
}

double Sphere::volume()
{
    return 4.*M_PI*radius*radius*radius/3.;
}

Point3D Sphere::nearestIntersect(const Point3D &point)
{
    Vector3D l;

    l = (point - center).unitVector();

    return center + radius*l;
}

std::pair<Point3D, Point3D> Sphere::lineIntersect(const Point3D &pt1, const Point3D &pt2)
{
    Vector3D l;
    double a, b, c, d1, d2;

    l = (pt2 - pt1).unitVector();

    a = dot(l, l);
    b = 2.*dot(l, pt1 - center);
    c = dot(pt1 - center, pt1 - center) - radius*radius;

    d1 = (-b + sqrt(b*b - 4.*a*c))/(2.*a);
    d2 = (-b - sqrt(b*b - 4.*a*c))/(2.*a);

    return std::pair<Point3D, Point3D>(pt1 + l*d1, pt1 + l*d2);
}

bool Sphere::isInside(const Point3D &point)
{
    if((point - center).mag() <= radius)
        return true;

    return false;
}

Sphere& Sphere::operator+=(const Vector3D& translationVec)
{
    center += translationVec;
    return *this;
}

Sphere& Sphere::operator-=(const Vector3D& translationVec)
{
    center -= translationVec;
    return *this;
}

Sphere& Sphere::rotate(const double radians)
{
    return *this;
}

Sphere& Sphere::scale(const double a)
{
    radius *= a;
    return *this;
}
