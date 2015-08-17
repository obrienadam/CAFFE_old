#define BOOST_TEST_MODULE Vector3DTest
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include "Vector3D.h"

BOOST_AUTO_TEST_SUITE (Vector3DTest)

BOOST_AUTO_TEST_CASE (test1)
{
  Vector3D u;

  srand(time(NULL));
  u = Vector3D(rand(), rand(), rand());

  BOOST_CHECK_CLOSE(u.unitVector().mag(), 1., 1e-13);
}

BOOST_AUTO_TEST_CASE (test2)
{
    Vector3D u(rand(), rand(), rand()), v;

    u = v;
    BOOST_CHECK_EQUAL(u.x, v.x);
    BOOST_CHECK_EQUAL(u.y, v.y);
    BOOST_CHECK_EQUAL(u.z, v.z);
    BOOST_CHECK(u == v);

    v = Vector3D(u.x*1.2, rand(), rand());
    BOOST_CHECK(u != v);
}

BOOST_AUTO_TEST_SUITE_END( )
