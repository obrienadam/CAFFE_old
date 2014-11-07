#include "AdvectionDiffusion.h"

AdvectionDiffusion::AdvectionDiffusion(double phi,
                                       double alpha,
                                       Vector3D v)
    :
      phi(phi),
      alpha(alpha),
      v(v)
{



}

AdvectionDiffusion& AdvectionDiffusion::operator+=(const AdvectionDiffusion& rhs)
{

    phi += rhs.phi;

    return *this;

}

AdvectionDiffusion& AdvectionDiffusion::operator-=(const AdvectionDiffusion& rhs)
{

    phi -= rhs.phi;

    return *this;

}

AdvectionDiffusion& AdvectionDiffusion::operator*=(double rhs)
{

    phi *= rhs;

    return *this;

}

AdvectionDiffusion& AdvectionDiffusion::operator/=(double rhs)
{

    phi /= rhs;

    return *this;

}

//- Functions

AdvectionDiffusion operator+(AdvectionDiffusion lhs, const AdvectionDiffusion& rhs)
{

    lhs += rhs;

    return lhs;

}

AdvectionDiffusion operator-(AdvectionDiffusion lhs, const AdvectionDiffusion& rhs)
{

    lhs -= rhs;

    return lhs;

}

AdvectionDiffusion operator*(AdvectionDiffusion lhs, double rhs)
{

    lhs *= rhs;

    return lhs;

}

AdvectionDiffusion operator*(double lhs, AdvectionDiffusion rhs)
{

    rhs *= lhs;

    return rhs;

}

AdvectionDiffusion operator/(AdvectionDiffusion lhs, double rhs)
{

    lhs /= rhs;

    return lhs;

}
