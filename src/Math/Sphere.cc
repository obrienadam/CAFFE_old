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

Point3D Sphere::nearestIntersect(const Point3D &origin, const Point3D &point)
{
    Vector3D l, oc;
    double a, b, d1, d2;

    l = (point - origin).unitVector();
    oc = origin - center;
    a = -dot(l, oc);
    b = sqrt(a*a - dot(oc, oc) + radius*radius);

    d1 = a + b;
    d2 = a - b;

    if(fabs(d1) < fabs(d2))
        return origin + d1*l;
    else
        return origin + d2*l;
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
