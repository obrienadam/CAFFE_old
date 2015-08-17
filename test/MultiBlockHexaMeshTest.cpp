#define BOOST_TEST_MODULE MultiBlockHexaMeshTest
#include <boost/test/included/unit_test.hpp>

#include "MultiBlockHexaFvmMesh.h"
#include "Parallel.h"

BOOST_AUTO_TEST_SUITE (MultiBlockHexaMeshTest)

BOOST_AUTO_TEST_CASE (test1)
{
    Parallel::initialize();



    Parallel::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
