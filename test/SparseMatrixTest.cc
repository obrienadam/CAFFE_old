#include <iostream>
#include <cstdlib>
#include <math.h>

#include "SparseMatrix.h"
#include "SparseVector.h"

#define SIZE 1200

int main()
{
    using namespace std;

    int i, j, nIters;
    int cols[3];
    double bValue[SIZE], vals[3];

    double h = 1./double(SIZE + 1);

    SparseMatrix A(SIZE, SIZE);
    SparseVector b(SIZE), x(SIZE);

    vals[0] = 1./(h*h);
    vals[1] = -2./(h*h);
    vals[2] = 1./(h*h);

    cols[0] = 0;
    cols[1] = 1;

    i = 0;

    A.setValues(1, &i, 2, cols, &vals[1]);

    for(i = 1; i < SIZE - 1; ++i)
    {
        cols[0] = i - 1;
        cols[1] = i;
        cols[2] = i + 1;

        A.setValues(1, &i, 3, cols, vals);
    }

    for(i = 0; i < SIZE; ++i)
    {
        bValue[i] = sin(2.*M_PI*h*(i + 1));

        if(i == 0)
            bValue[0] += -1./(h*h);

        b.setValues(1, &i, &bValue[i]);
    }

    i = SIZE - 1;
    cols[0] = i - 1;
    cols[1] = i;

    A.setValues(1, &i, 2, cols, vals);
    A.assemble();
    b.assemble();

    nIters = A.solve(b, x);

    x.print();

    cout << "The problem was solved in " << nIters << " iterations." << endl;

	return 0;
}
