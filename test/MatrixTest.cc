#include <iostream>

#include "../src/RunControl/Output.h"
#include "../src/Math/Matrix.h"

using namespace std;

int main()
{
  using namespace std;

  Matrix A(10, 10), B, C, b, x;
  
  cout << "An identity matrix:" << endl;

  A = eye(5);
  A.print();

  cout << endl << "A random matrix:" << endl;

  A = random(8, 8, -5, 5);
  A.print();

  cout << endl << "Another random matrix:" << endl;

  B = random(8, 8, -6, 9);
  B.print();

  cout << endl <<"Test a known solution:" << endl;

  A.allocate(3, 3);
  A(0, 0) = 3; A(0, 1) = 2; A(0, 2) = -1;
  A(1, 0) = 2; A(1, 1) = -2; A(1, 2) = 4;
  A(2, 0) = -1; A(2, 1) = 0.5; A(2, 2) = -1;
  
  A.print();
  
  B.allocate(3, 1);
  B(0, 0) = 1;
  B(1, 0) = -2;
  B(2, 0) = 0;

  B.print();

  C = solve(A, B);
  
  C.print();

  cout << endl << "Testing the solve function:" << endl;

  C = solve(A, B);
  C.print();

  cout << endl << "Testing the least squares function by solving over-determined system:" << endl;

  B = random(10, 3, 0, 9);
  A = random(10, 8, 0, 4);
  C = solveLeastSquares(A, B);
  C.print();

}
