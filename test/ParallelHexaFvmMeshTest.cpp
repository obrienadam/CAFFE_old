#define BOOST_TEST_MODULE ParallelHexaFvmMeshTest
#include <boost/test/included/unit_test.hpp>

#include "ParallelHexaFvmMesh.h"

BOOST_AUTO_TEST_SUITE (ParallelHexaFvmMeshTest)

BOOST_AUTO_TEST_CASE (test1)
{
    Parallel::initialize();

    ParallelHexaFvmMesh mesh;
    mesh.initializeCartesianMesh(1, 1, 1, 64, 64, 64);
    mesh.name = "test1Mesh";

    BOOST_CHECK_EQUAL(mesh.nCells(), 64*64*64/8);

    //- Output a debug file
    mesh.writeTec360(0, "");
}

BOOST_AUTO_TEST_CASE (test2)
{
    ParallelHexaFvmMesh mesh;

    mesh.name = "test2Mesh";
    mesh.initializeCartesianMesh(1, 1, 1, 112, 112, 112);

    mesh.writeTec360(1, "");

    Parallel::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
