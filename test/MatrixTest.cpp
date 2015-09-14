#define BOOST_TEST_MODULE MatrixTest
#include <boost/test/included/unit_test.hpp>
#include "Matrix.h"

BOOST_AUTO_TEST_SUITE (MatrixTest)

BOOST_AUTO_TEST_CASE (test1)
{
    Matrix A = random(10, 8, 0, 15);
    Matrix B = A;

    BOOST_REQUIRE_EQUAL(A.nRows(), 10);
    BOOST_REQUIRE_EQUAL(A.nCols(), 8);
    BOOST_REQUIRE_EQUAL(A.nRows(), B.nRows());
    BOOST_REQUIRE_EQUAL(A.nCols(), B.nCols());
    BOOST_REQUIRE(A.data() != B.data());

    for(int i = 0; i < A.nRows(); ++i)
        for(int j = 0; j < A.nCols(); ++j)
            BOOST_REQUIRE_EQUAL(A(i, j), B(i, j));
}

BOOST_AUTO_TEST_CASE (test2)
{
    Matrix A = random(10, 10, 0, 15);
    Matrix I = A*inverse(A);

    for(int i = 0; i < A.nRows(); ++i)
        for(int j = 0; j < A.nCols(); ++j)
        {
            if(i == j)
                BOOST_REQUIRE_CLOSE(I(i, j), 1., 1e-13);
            else
                BOOST_REQUIRE_CLOSE(I(i, j) + 1, 1, 1e-13);
        }
}

BOOST_AUTO_TEST_CASE (test3)
{
    Matrix A = random(10, 10, 1, 15);
    Matrix B = random(10, 1, 1, 10);
    Matrix X1 = solve(A, B);
    Matrix X2 = inverse(A)*B;

    A.solve(B);
    BOOST_REQUIRE_EQUAL(B.nRows(), 10);
    BOOST_REQUIRE_EQUAL(B.nCols(), 1);

    for(int i = 0; i < X1.nRows(); ++i)
    {
        for(int j = 0; j < X1.nCols(); ++j)
        {
            BOOST_REQUIRE_CLOSE(X1(i, j), X2(i, j), 1e-12);
            BOOST_REQUIRE_CLOSE(X1(i, j), B(i, j), 1e-12);
        }
    }
}

BOOST_AUTO_TEST_CASE (test4)
{
    Matrix A = random(20, 10, 1, 10);
    Matrix B = random(20, 1, 1, 10);
    Matrix X1 = solveLeastSquares(A, B);
    Matrix X2;

    BOOST_REQUIRE_EQUAL(X1.nRows(), 10);
    BOOST_REQUIRE_EQUAL(X1.nCols(), 1);

    X2 = solve(transpose(A)*A, transpose(A)*B);

    for(int i = 0; i < X2.nRows(); ++i)
    {
        for(int j = 0; j < X2.nCols(); ++j)
        {
            BOOST_REQUIRE_CLOSE(X1(i, j), X2(i, j), 1e-11);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
