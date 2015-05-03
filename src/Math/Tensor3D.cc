#include <algorithm>

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

double& Tensor3D::operator()(int i, int j)
{
    switch(i)
    {
    case 0:

        switch(j)
        {
        case 0:
            return xx;
        case 1:
            return xy;
        case 2:
            return xz;
        };

        break;
    case 1:

        switch(j)
        {
        case 0:
            return yx;
        case 1:
            return yy;
        case 2:
            return yz;
        };

        break;
    case 2:

        switch(j)
        {
        case 0:
            return zx;
        case 1:
            return zy;
        case 2:
            return zz;
        };
    };

    throw("Attempted to access element outside the bounds of Tensor3D.");
}

double& Tensor3D::operator()(int i)
{
    switch(i)
    {
    case 0:
        return xx;
    case 1:
        return xy;
    case 2:
        return xz;
    case 3:
        return yx;
    case 4:
        return yy;
    case 5:
        return yz;
    case 6:
        return zx;
    case 7:
        return zy;
    case 8:
        return zz;
    };

    throw("Attempted to access element outside the bounds of Tensor3D.");
}

Vector3D Tensor3D::row(int rowNo)
{
    switch(rowNo)
    {
    case 0:
        return Vector3D(xx, xy, xz);
    case 1:
        return Vector3D(yx, yy, yz);
    case 2:
        return Vector3D(zx, zy, zz);
    };

    throw("Attempted to access row outside the bounds of Tensor3D.");
}

Tensor3D& Tensor3D::transpose()
{
    std::swap(xy, yx);
    std::swap(xz, zx);
    std::swap(yz, zy);

    return *this;
}

Tensor3D& Tensor3D::operator+=(const Tensor3D& rhs)
{
    xx += rhs.xx;
    xy += rhs.xy;
    xz += rhs.xz;
    yx += rhs.yx;
    yy += rhs.yy;
    yz += rhs.yz;
    zx += rhs.zx;
    zy += rhs.zy;
    zz += rhs.zz;

    return *this;
}

Tensor3D& Tensor3D::operator-=(const Tensor3D& rhs)
{
    xx -= rhs.xx;
    xy -= rhs.xy;
    xz -= rhs.xz;
    yx -= rhs.yx;
    yy -= rhs.yy;
    yz -= rhs.yz;
    zx -= rhs.zx;
    zy -= rhs.zy;
    zz -= rhs.zz;

    return *this;
}

void Tensor3D::print()
{
    using namespace std;

    for(int j = 0; j < 3; ++j)
    {
        for(int i = 0; i < 3; ++i)
        {
            cout << operator()(i, j) << ", ";
        }

        cout << endl;
    }
}

Vector3D operator*(const Tensor3D& lhs, const Vector3D& rhs)
{
    return Vector3D(lhs.xx*rhs.x + lhs.xy*rhs.y + lhs.xz*rhs.z,
                    lhs.yx*rhs.x + lhs.yy*rhs.y + lhs.yz*rhs.z,
                    lhs.zx*rhs.x + lhs.zy*rhs.y + lhs.zz*rhs.z);
}

Tensor3D operator*(const Tensor3D& lhs, const Tensor3D& rhs)
{
    return Tensor3D(lhs.xx*rhs.xx + lhs.xy*rhs.yx + lhs.xz*rhs.zx, lhs.yx*rhs.xx + lhs.yy*rhs.yx + lhs.yz*rhs.zx, lhs.zx*rhs.xx + lhs.zy*rhs.yx + lhs.zz*rhs.zx,
                    lhs.xx*rhs.xy + lhs.xy*rhs.yy + lhs.xz*rhs.zy, lhs.yx*rhs.xy + lhs.yy*rhs.yy + lhs.yz*rhs.zy, lhs.zx*rhs.xy + lhs.zy*rhs.yy + lhs.zz*rhs.zy,
                    lhs.xx*rhs.xz + lhs.xy*rhs.yz + lhs.xz*rhs.zz, lhs.yx*rhs.xz + lhs.yy*rhs.yz + lhs.yz*rhs.zz, lhs.zx*rhs.xz + lhs.zy*rhs.yz + lhs.zz*rhs.zz);
}

Tensor3D operator+(Tensor3D lhs, const Tensor3D& rhs)
{
    lhs += rhs;
    return lhs;
}

Tensor3D operator-(Tensor3D lhs, const Tensor3D& rhs)
{
    lhs -= rhs;
    return lhs;
}

Tensor3D transpose(Tensor3D tensor)
{
    tensor.transpose();
    return tensor;
}

Vector3D dot(const Tensor3D &tensor, const Vector3D &vec)
{
    return Vector3D(tensor.xx*vec.x + tensor.xy*vec.y + tensor.xz*vec.z,
                    tensor.yx*vec.x + tensor.yy*vec.y + tensor.yz*vec.z,
                    tensor.zx*vec.x + tensor.zy*vec.y + tensor.zz*vec.z);
}

double trace(const Tensor3D &tensor)
{
    return tensor.xx + tensor.yy + tensor.zz;
}

std::ostream& operator<<(std::ostream& os, const Tensor3D& tensor)
{
    using namespace std;

    os << "(" << tensor.xx << ", " << tensor.xy << ", " << tensor.xz << ")" << endl
       << "(" << tensor.yx << ", " << tensor.yy << ", " << tensor.yz << ")" << endl
       << "(" << tensor.zx << ", " << tensor.zy << ", " << tensor.zz << ")" << endl;

    return os;
}
