#include <math.h>

#include "Vector3D.h"

Vector3D::Vector3D(double x, double y, double z) :
  x(x),
  y(y),
  z(z)
{

}

double& Vector3D::operator()(int i)
{
  if(i < 0 || i > 2)
    throw "Attempted to acces element outside the bounds of Vector3D.";
  
  switch(i)
    {
    case 0:

      return x;

    case 1:

      return y;

    default:

      return z;

    };
}

double Vector3D::mag()
{
  return sqrt(x*x + y*y + z*z);
}

Vector3D Vector3D::unitVector()
{
  double invMag(1./mag());

  return Vector3D(invMag*x, invMag*y, invMag*z);
}

Vector3D& Vector3D::operator+=(const Vector3D& rhs)
{
  x += rhs.x;
  y += rhs.y;
  z += rhs.z;

  return *this;
}

Vector3D& Vector3D::operator-=(const Vector3D& rhs)
{
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;

  return *this;
}

Vector3D& Vector3D::operator*=(double rhs)
{
  x *= rhs;
  y *= rhs;
  z *= rhs;

  return *this;
}

Vector3D& Vector3D::operator/=(double rhs)
{
  rhs = 1./rhs;

  x *= rhs;
  y *= rhs;
  z *= rhs;

  return *this;
}

Vector3D operator+(Vector3D lhs, const Vector3D& rhs)
{
  return lhs += rhs;
}

Vector3D operator-(Vector3D lhs, const Vector3D& rhs)
{
  return lhs -= rhs;
}

Vector3D operator*(Vector3D lhs, double rhs)
{
  return lhs *= rhs;
}

Vector3D operator*(double lhs, Vector3D rhs)
{
  return rhs *= lhs;
}

Vector3D operator/(Vector3D lhs, double rhs)
{
  return lhs /= rhs;
}

double dotProduct(const Vector3D& u, const Vector3D& v)
{
  return u.x*v.x + u.y*v.y + u.z*v.z;
}

Vector3D crossProduct(const Vector3D& u, const Vector3D& v)
{
  return Vector3D(u.y*v.z - u.z*v.y,
		  u.z*v.x - u.x*v.z,
		  u.x*v.y - u.y*v.x);
}

Vector3D relativeVector(const Vector3D& u, const Vector3D& v)
{

  return v - u;

}
