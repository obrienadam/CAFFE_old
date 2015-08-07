#include <math.h>

#include "Kernel.h"

Kernel::Kernel(Type type, double epsilon)
    :
      type_(type),
      epsilon_(epsilon)
{
    computeA();
}

void Kernel::changeType(Type newType)
{
    type_ = newType;
    computeA();
}

double Kernel::value(const Point3D &center, const Point3D &point)
{
    double radius = (center - point).mag();

    if(radius <= epsilon_)
    {
        switch(type_)
        {
        case K6:
            return A_*pow(epsilon_*epsilon_ - radius*radius, 3);
        case K8:
            return A_*pow(epsilon_*epsilon_ - radius*radius, 4);
        };
    }

    return 0.;
}

void Kernel::computeA()
{
    switch(type_)
    {
    case K6:
        A_ = 315./(64.*M_PI*pow(epsilon_, 9));
        break;
    case K8:
        A_ = 3465./(128.*M_PI*pow(epsilon_, 11));
        break;
    };
}
