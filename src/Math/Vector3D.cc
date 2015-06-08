#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include "Vector3D.h"

Vector3D::Vector3D(double x, double y, double z) :
  x(x),
  y(y),
  z(z)
{

}

Vector3D::Vector3D(const Vector3D &other)
    :
      x(other.x),
      y(other.y),
      z(other.z)
{

}

double& Vector3D::operator()(int i)
{ 
    switch(i)
    {
    case 0:
        return x;
    case 1:
        return y;
    case 2:
        return z;
    };

    throw "Attempted to access an element outside the range of Vector3D.";
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

void Vector3D::print()
{
    std::cout << "(" << x << ", " << y << ", " << z << ")\n";
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

Vector3D operator-(const Vector3D& rhs)
{
    return Vector3D(-rhs.x, -rhs.y, -rhs.z);
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

Vector3D sqr(const Vector3D &u)
{
    return Vector3D(u.x*u.x, u.y*u.y, u.z*u.z);
}

Vector3D sqrt(const Vector3D& u)
{
    return Vector3D(sqrt(u.x), sqrt(u.y), sqrt(u.z));
}

std::ostream& operator<<(std::ostream& os, const Vector3D& vec)
{
    os << vec.x << "," << vec.y << "," << vec.z;

    return os;
}

double dot(const Vector3D& u, const Vector3D& v)
{
  return u.x*v.x + u.y*v.y + u.z*v.z;
}

Tensor3D tensor(const Vector3D &u, const Vector3D &v)
{
    return Tensor3D(u.x*v.x, u.x*v.y, u.x*v.z,
                    u.y*v.x, u.y*v.y, u.y*v.z,
                    u.z*v.x, u.z*v.y, u.z*v.z);
}

Vector3D cross(const Vector3D& u, const Vector3D& v)
{
  return Vector3D(u.y*v.z - u.z*v.y,
		  u.z*v.x - u.x*v.z,
		  u.x*v.y - u.y*v.x);
}

Vector3D relativeVector(const Vector3D& u, const Vector3D& v)
{
  return v - u;
}

namespace std
{
Vector3D max(Vector3D u, const Vector3D &v)
{
    u.x = max(u.x, v.x);
    u.y = max(u.y, v.y);
    u.z = max(u.z, v.z);

    return u;
}

Vector3D min(Vector3D u, const Vector3D &v)
{
    u.x = min(u.x, v.x);
    u.y = min(u.y, v.y);
    u.z = min(u.z, v.z);

    return u;
}

Vector3D stov(string vecStr)
{
    Vector3D vec;
    string nextElement, delim = ", ";
    vecStr = vecStr.substr(vecStr.find_first_of("(") + 1, vecStr.find_first_of(")") - 1);

    //- Extract x-component
    vecStr = vecStr.substr(vecStr.find_first_not_of(delim), vecStr.length());
    nextElement = vecStr.substr(0, vecStr.find_first_of(delim));
    vec.x = stod(nextElement);

    //- Extract y-component
    vecStr = vecStr.substr(vecStr.find_first_of(delim), vecStr.length());
    vecStr = vecStr.substr(vecStr.find_first_not_of(delim), vecStr.length());
    nextElement = vecStr.substr(0, vecStr.find_first_of(delim));
    vec.y = stod(nextElement);

    //- Extract z-component
    vecStr = vecStr.substr(vecStr.find_first_of(delim), vecStr.length());
    vecStr = vecStr.substr(vecStr.find_first_not_of(delim), vecStr.length());
    nextElement = vecStr.substr(0, vecStr.find_first_of(delim));
    vec.z = stod(nextElement);

    return vec;
}
}
