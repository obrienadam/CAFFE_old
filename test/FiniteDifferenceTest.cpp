#include <iostream>

#include "../src/State/AdvectionDiffusion.h"
#include "../src/Math/FiniteDifference.h"

int main()
{
    using namespace std;

    double x[3], fx[3], dfxdx[3];

    x[0] = 0.;
    x[1] = 0.000002;
    x[2] = 0.000003;

    //- Test 2nd order central

    for(int i = 0; i < 3; ++i)
    {

        fx[i] = x[i]*x[i]*x[i];
        dfxdx[i] = 6.*x[i];

    }

    cout << "Analytical x^3 2nd derivative at x = 0.2: " << dfxdx[1] << endl
         << "Numerical x^3 2nd derivative at x = 0.2: " << fd2ndDeriv2ndOrderCentral(fx[0], fx[1], fx[2], x[1] - x[0], x[2] - x[1]) << endl;

    return 0;

}
