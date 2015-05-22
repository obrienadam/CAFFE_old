#include <iostream>

#include "../src/RunControl/Output.h"
#include "../src/Math/Matrix.h"
#include "Point3D.h"
#include "Interpolation.h"

int main()
{
    using namespace std;

    double a[16] = {1, 2, 4, 3, 2, 3, 4, 1, 9, 5, 3, 2, 1, 5, 3, 5};
    double b[4] = {1, 2, 3, 4};
    double c[4] = {-0.1538, 1.0000, -0.1538, -0.0769};
    Matrix A(a, 4, 4), B(b, 4, 1), C(c, 4, 1);
    Point3D points[12], pt;
    double values[12];

    try{

        cout << "Allocating 5x5 matrices..." << endl;

        cout << "A contains:" << endl;

        A.print();

        cout << "B contains:" << endl;

        B.print();

        cout << "C is the solution of AC = B and contains:" << endl;

        C.print();

        cout << "Testing the linear solver (should be the same vector as above):" << endl;

        C = solve(A, B);
        C.print();

        cout << "A still contains:" << endl;

        A.print();

        cout << "Inverting A:" << endl;

        C = inverse(A);
        C.print();

        cout << "Multiplying A*A^-1:" << endl;

        B = A*C;
        B.print();

        cout << "Multiplying A*A-1*A:" << endl;

        C = B*A;
        C.print();

        cout << "Allocating a large matrix:" << endl;

        A = random(1000, 1000, 0, 9);
        B = random(1000, 1000, -5, 9);

        cout << "Solving large matrix with a large number of rhs vectors:" << endl;

        A.solve(B);

        cout << "Setting up a least squares problem:" << endl;

        A = random(10, 3, -9, 9);
        B = random(10, 1, -9, 9);

        cout << "Matrix A:" << endl;

        A.print();

        cout << "Matrix B:" << endl;

        B.print();

        C = solveLeastSquares(A, B);

        cout << "Solution:" << endl;

        C.print();

        cout << "Solution when computed second way:" << endl;

        B = transpose(A)*B;
        A = transpose(A)*A;
        C = solve(A, B);
        C.print();

        cout << "Test the interpolation methods:" << endl;

        cout << "Testing 4 point stencil with linear:\n";

        points[0] = Point3D(1., 1., 2); values[0] = 1.;
        points[1] = Point3D(1., 3., 2); values[1] = 2.;
        points[2] = Point3D(3., 1., 2); values[2] = 3.;
        points[3] = Point3D(3., 3., 2); values[3] = 4.;

        for(int i = 0; i < 10; ++i)
        {
            pt.x = rand()/double(RAND_MAX)*2. + 1.;
            pt.y = rand()/double(RAND_MAX)*2. + 1.;
            pt.z = rand()/double(RAND_MAX)*2. + 1.;

            cout << "At " << pt << ": " << Interpolation::linear(points, values, 4, pt) << endl;
        }

        cout << "Testing 12 point stencil with linear and quadratic:\n";

        points[0] = Point3D(0, 0, 0); values[0] = 1;
        points[1] = Point3D(1, 0, 0); values[1] = 2;
        points[2] = Point3D(1, 1, 0); values[2] = 3;
        points[3] = Point3D(0, 1, 0); values[3] = 4;

        points[4] = Point3D(0, 0, 2); values[4] = 9;
        points[5] = Point3D(1, 0, 2); values[5] = 10;
        points[6] = Point3D(1, 1, 2); values[6] = 11;
        points[7] = Point3D(0, 1, 2); values[7] = 12;

        points[8] = Point3D(0, 0, 1); values[8] = 5;
        points[9] = Point3D(1, 0, 1); values[9] = 6;
        points[10] = Point3D(1, 1, 1); values[10] = 7;
        points[11] = Point3D(0, 1, 1); values[11] = 8;

        double a = 0.0001;

        for(int i = 0; i < 10; ++i)
        {
            points[i] *= a;
            values[i] *= 0.0001*a;
        }
        for(int i = 0; i < 10; ++i)
        {
            pt.x = rand()/double(RAND_MAX);
            pt.y = rand()/double(RAND_MAX);
            pt.z = rand()/double(RAND_MAX) + 1.;

            pt *= a;

            cout << "Linear at " << pt << ": " << Interpolation::linear(points, values, 12, pt) << endl;
            cout << "Quadratic at " << pt << ": " << Interpolation::quadratic(points, values, 12, pt) << endl;
        }
    }
    catch (const char* errorMessage)
    {
        cerr << errorMessage << endl;
    }

    return 0;
}
