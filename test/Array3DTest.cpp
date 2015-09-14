#define BOOST_TEST_MODULE Array3DTest
#include <boost/test/included/unit_test.hpp>
#include "Array3D.h"

BOOST_AUTO_TEST_SUITE (Array3DTest)

BOOST_AUTO_TEST_CASE (test1)
{
    Array3D<double> array3D(100, 100, 100), array3D2;

    for (auto it = array3D.begin(); it != array3D.end(); ++it)
    {
        *it = rand();
    }

    array3D.resize(20, 20, 20);

    for(auto it = array3D.begin(); it != array3D.end(); ++it)
    {
        *it = rand();
    }

    array3D2 = array3D;
    Array3D<double> array3D3(array3D);

    BOOST_REQUIRE_EQUAL(array3D2.sizeI(), array3D.sizeI());
    BOOST_REQUIRE_EQUAL(array3D2.sizeJ(), array3D.sizeJ());
    BOOST_REQUIRE_EQUAL(array3D2.sizeK(), array3D.sizeK());
    BOOST_REQUIRE_EQUAL(array3D2.size(), array3D.size());
    BOOST_REQUIRE_EQUAL(array3D.sizeI(), array3D3.sizeI());
    BOOST_REQUIRE_EQUAL(array3D.sizeJ(), array3D3.sizeJ());
    BOOST_REQUIRE_EQUAL(array3D.sizeK(), array3D3.sizeK());
    BOOST_REQUIRE_EQUAL(array3D.size(), array3D3.size());

    for(auto it = array3D.begin(), it2 = array3D2.begin(); it != array3D.end(); ++it, ++it2)
        BOOST_REQUIRE_EQUAL(*it, *it2);

    array3D(3, 3, 3) = 5;
    BOOST_REQUIRE_EQUAL(array3D(3, 3, 3), 5);
}

BOOST_AUTO_TEST_CASE (test2)
{
    Array3D<double> array3D(20, 20, 20);

    array3D.assign(1);

    BOOST_REQUIRE_EQUAL(array3D.size(), 20*20*20);

    for(auto it = array3D.begin(); it != array3D.end(); ++it)
        BOOST_REQUIRE_EQUAL(*it, 1);
}

BOOST_AUTO_TEST_CASE (test3)
{
    Array3D<double> array3D;

    array3D.resize(10, 5, 3);
    array3D.assign(3);

    auto it = array3D.begin();

    for(int k = 0; k < array3D.sizeK(); ++k)
    {
        for(int j = 0; j < array3D.sizeJ(); ++j)
        {
            for(int i = 0; i < array3D.sizeI(); ++i)
            {
                BOOST_REQUIRE_EQUAL(array3D(i, j, k), *it);
                ++it;
            }
        }
    }

    array3D.resize(20, 30, 10);
    array3D.assign(2);

    it = array3D.begin();
    *it = 3;

    for(int k = 0; k < array3D.sizeK(); ++k)
    {
        for(int j = 0; j < array3D.sizeJ(); ++j)
        {
            for(int i = 0; i < array3D.sizeI(); ++i)
            {
                BOOST_REQUIRE_EQUAL(array3D(i, j, k), *it);
                ++it;
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
