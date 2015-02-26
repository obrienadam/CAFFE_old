#ifndef TENSOR_3D_H
#define TENSOR_3D_H

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

  double operator()(int i, int j);

  double mag();

  Tensor3D& operator+=(const Tensor3D& rhs);
  Tensor3D& operator-=(const Tensor3D& rhs);
  Tensor3D& operator*=(const double& rhs);
  Tensor3D& operator/=(const double& rhs);
};

Vector3D operator*(const Tensor3D& lhs, const Vector3D& rhs);

#endif
