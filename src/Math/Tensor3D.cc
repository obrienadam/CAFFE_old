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

std::ostream& operator<<(std::ostream& os, const Tensor3D& tensor)
{
    using namespace std;

    os << "(" << tensor.xx << ", " << tensor.xy << ", " << tensor.xz << ")" << endl
       << "(" << tensor.yx << ", " << tensor.yy << ", " << tensor.yz << ")" << endl
       << "(" << tensor.zx << ", " << tensor.zy << ", " << tensor.zz << ")" << endl;

    return os;
}
