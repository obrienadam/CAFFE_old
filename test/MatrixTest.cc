#include <iostream>

#include "../src/RunControl/Output.h"
#include "../src/Math/Matrix.h"

int main()
{
    using namespace std;

    double a[16] = {1, 2, 4, 3, 2, 3, 4, 1, 9, 5, 3, 2, 1, 5, 3, 5};
    double b[4] = {1, 2, 3, 4};
    double c[4] = {-0.1538, 1.0000, -0.1538, -0.0769};
    Matrix A(a, 4, 4), B(b, 4, 1), C(c, 4, 1);

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

    }
    catch (const char* errorMessage)
    {
        cerr << errorMessage << endl;
    }

    return 0;
}
