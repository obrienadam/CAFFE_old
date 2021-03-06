#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>
#include <string>
#include <sstream>

#include "Tensor3D.h"

class Tensor3D;

class Vector3D
{
 public:

  Vector3D(double x = 0., double y = 0., double z = 0.);
  Vector3D(const Vector3D& other);

  double& operator()(int i);
  const double& operator()(int i) const;

  double mag() const;
  Vector3D unitVector() const;

  //- Print to screen
  void print() const;

  //- Member operators
  Vector3D& operator+=(const Vector3D& rhs);
  Vector3D& operator-=(const Vector3D& rhs);
  Vector3D& operator*=(double rhs);
  Vector3D& operator/=(double rhs);

  double x, y, z;
};

//- Vector3D related operators
Vector3D operator+(Vector3D lhs, const Vector3D& rhs);
Vector3D operator-(Vector3D lhs, const Vector3D& rhs);
Vector3D operator-(const Vector3D& rhs);
Vector3D operator*(Vector3D lhs, double rhs);
Vector3D operator*(double lhs, Vector3D rhs);
Vector3D operator/(Vector3D lhs, double rhs);

bool operator ==(const Vector3D &lhs, const Vector3D &rhs);
bool operator !=(const Vector3D &lhs, const Vector3D &rhs);

Vector3D sqr(const Vector3D& u);
Vector3D sqrt(const Vector3D &u);

std::ostream& operator<<(std::ostream& os, const Vector3D& vec);

//- Vector3D related functions
double dot(const Vector3D& u, const Vector3D& v);
Tensor3D tensor(const Vector3D& u, const Vector3D& v);
Vector3D cross(const Vector3D& u, const Vector3D& v);
Vector3D relativeVector(const Vector3D& u, const Vector3D& v);

//- function with the std interface
namespace std
{
Vector3D max(Vector3D u, const Vector3D& v);
Vector3D min(Vector3D u, const Vector3D& v);
Vector3D stov(string vecStr);
string to_string(const Vector3D& u);
}

#endif
