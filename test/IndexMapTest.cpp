#define BOOST_TEST_MODULE IndexMapTest
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include "ParallelHexaFvmMesh.h"

BOOST_AUTO_TEST_SUITE (IndexMapTest)

BOOST_AUTO_TEST_CASE (test1)
{
    Parallel::initialize();

    ParallelHexaFvmMesh mesh;
    mesh.changeName("indexMapTest1Mesh");
    mesh.initializeCartesianMesh(2, 1.5, 1.2, 36, 28, 52);

    BOOST_CHECK_EQUAL(mesh.iMap.nActiveGlobal(), 36*28*52);
    BOOST_CHECK_EQUAL(mesh.iMap.nActiveLocal(), 36*28*52/8);
}

BOOST_AUTO_TEST_CASE (test2)
{
    ParallelHexaFvmMesh mesh;
    int testIndex;
    mesh.changeName("indexMapTest2Mesh");
    mesh.initializeCartesianMesh(1, 1, 1, 40, 40, 40);

    testIndex = Parallel::broadcast(mesh.iMap(0, 0, 0, 0), 1);

    if(Parallel::isMainProcessor())
    {
        BOOST_CHECK_EQUAL(mesh.iMap(mesh.uCellI() + 1, 0, 0, 0), testIndex);
    }

    testIndex = Parallel::broadcast(mesh.iMap(0, 0, 0, 0), 2);

    if(Parallel::isMainProcessor())
    {
        BOOST_CHECK_EQUAL(mesh.iMap(0, mesh.uCellJ() + 1, 0, 0), testIndex);
    }

    testIndex = Parallel::broadcast(mesh.iMap(0, 0, 0, 0), 4);

    if(Parallel::isMainProcessor())
    {
        BOOST_CHECK_EQUAL(mesh.iMap(0, 0, mesh.uCellK() + 1, 0), testIndex);
    }
}

BOOST_AUTO_TEST_CASE (end)
{
    Parallel::finalize();
}

BOOST_AUTO_TEST_SUITE_END( )
