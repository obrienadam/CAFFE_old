#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>
#include <string>
#include <sstream>

#include "Tensor3D.h"

class Tensor3D;

class Vector3D
{
 private:

 public:

  Vector3D(double x = 0., double y = 0., double z = 0.);
  Vector3D(std::string vectorStr);

  double x, y, z;

  double& operator()(int i);

  double mag();
  Vector3D unitVector();

  //- Special initializations

  void initialize(std::string vectorStr);

  //- Print to screen

  void print();

  //- Member operators

  Vector3D& operator+=(const Vector3D& rhs);
  Vector3D& operator-=(const Vector3D& rhs);
  Vector3D& operator*=(double rhs);
  Vector3D& operator/=(double rhs);

};

//- Vector3D related operators

Vector3D operator+(Vector3D lhs, const Vector3D& rhs);
Vector3D operator-(Vector3D lhs, const Vector3D& rhs);
Vector3D operator-(const Vector3D& rhs);
Vector3D operator*(Vector3D lhs, double rhs);
Vector3D operator*(double lhs, Vector3D rhs);
Vector3D operator/(Vector3D lhs, double rhs);
std::ostream& operator<<(std::ostream& os, const Vector3D& vec);

//- Vector3D related functions

double dot(const Vector3D& u, const Vector3D& v);
Tensor3D tensor(const Vector3D& u, const Vector3D& v);
Vector3D cross(const Vector3D& u, const Vector3D& v);
Vector3D relativeVector(const Vector3D& u, const Vector3D& v);

#endif
