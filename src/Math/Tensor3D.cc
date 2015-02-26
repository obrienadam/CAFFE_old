#include "Tensor3D.h"

Tensor3D::Tensor3D(double xx,
                   double xy,
                   double xz,
                   double yx,
                   double yy,
                   double yz,
                   double zx,
                   double zy,
                   double zz)
    :
      xx(xx),
      xy(xy),
      xz(xz),
      yx(yx),
      yy(yy),
      yz(yz),
      zx(zx),
      zy(zy),
      zz(zz)
{

}

Vector3D operator*(const Tensor3D& lhs, const Vector3D& rhs)
{
    return Vector3D(lhs.xx*rhs.x + lhs.xy*rhs.y + lhs.xz*rhs.z,
                    lhs.yx*rhs.x + lhs.yy*rhs.y + lhs.yz*rhs.z,
                    lhs.zx*rhs.x + lhs.zy*rhs.z + lhs.zz*rhs.z);
}
