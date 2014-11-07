#include "AdvectionDiffusion.h"
#include "FiniteDifference.h"

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

AdvectionDiffusion evaluateGradient(const AdvectionDiffusion& phiMinus1, double xMinus1,
				    const AdvectionDiffusion& phi, double x,
				    const AdvectionDiffusion& phiPlus1, double xPlus1)
{

    return AdvectionDiffusion(fd2ndDeriv2ndOrderCentral(phiMinus1.phi,
							phi.phi,
							phiPlus1.phi,
							x - xMinus1,
							xPlus1 - x));

}
