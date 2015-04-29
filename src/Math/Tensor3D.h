#ifndef TENSOR_3D_H
#define TENSOR_3D_H

#include <iostream>

#include "Vector3D.h"

class Vector3D;

class Tensor3D
{
 private:

 public:

  Tensor3D(double xx = 0.,
	   double xy = 0.,
	   double xz = 0.,
	   double yx = 0.,
	   double yy = 0.,
	   double yz = 0.,
	   double zx = 0.,
	   double zy = 0.,
	   double zz = 0.);

  double xx, xy, xz, yx, yy, yz, zx, zy, zz;

  double& operator()(int i, int j);
  double& operator()(int i);
  Vector3D row(int rowNo);

  double mag();
  Tensor3D& transpose();

  Tensor3D& operator+=(const Tensor3D& rhs);
  Tensor3D& operator-=(const Tensor3D& rhs);
  Tensor3D& operator*=(const double& rhs);
  Tensor3D& operator/=(const double& rhs);

  void print();
};

Vector3D operator*(const Tensor3D& lhs, const Vector3D& rhs);
Tensor3D operator*(const Tensor3D& lhs, const Tensor3D& rhs);
Tensor3D operator+(Tensor3D lhs, const Tensor3D& rhs);
Tensor3D operator-(Tensor3D lhs, const Tensor3D& rhs);

Tensor3D transpose(Tensor3D tensor);
Vector3D dot(const Tensor3D& tensor, const Vector3D& vec);

std::ostream& operator<<(std::ostream& os, const Tensor3D& tensor);

#endif
