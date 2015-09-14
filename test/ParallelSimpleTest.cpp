#define BOOST_TEST_MODULE ParallelSimpleTest
#include <boost/test/included/unit_test.hpp>

#include "ParallelHexaFvmMesh.h"
#include "Simple.h"

BOOST_AUTO_TEST_SUITE (ParallelSimpleTest)

BOOST_AUTO_TEST_CASE (test1)
{
    Parallel::initialize();

    ParallelHexaFvmMesh mesh;
    mesh.initializeCartesianMesh(1, 1, 1, 20, 20, 20);
}

BOOST_AUTO_TEST_CASE (end)
{
    Parallel::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
