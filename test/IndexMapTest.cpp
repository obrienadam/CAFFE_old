#define BOOST_TEST_MODULE IndexMapTest
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include "ParallelHexaFvmMesh.h"
#include "SparseMatrix.h"

BOOST_AUTO_TEST_SUITE (IndexMapTest)

BOOST_AUTO_TEST_CASE (test1)
{
    Parallel::initialize();

    ParallelHexaFvmMesh mesh;
    mesh.changeName("indexMapTest1Mesh");
    mesh.initializeCartesianMesh(2, 1.5, 1.2, 36, 28, 52);

    BOOST_REQUIRE_EQUAL(mesh.iMap.nActiveGlobal(), 36*28*52);
    BOOST_REQUIRE_EQUAL(mesh.iMap.nActiveLocal(), 36*28*52/8);
}

BOOST_AUTO_TEST_CASE (test2)
{
    ParallelHexaFvmMesh mesh;
    mesh.changeName("indexMapTest3Mesh");
    mesh.initializeCartesianMesh(1, 1, 1, 20, 20, 20);
    int cellsPerProc = mesh.nCellsI()*mesh.nCellsJ()*mesh.nCellsK();

    // Check that all the global indices are correct

    for(int k = 0; k < mesh.nCellsK(); ++k)
    {
        for(int j = 0; j < mesh.nCellsJ(); ++j)
        {
            for(int i = 0; i < mesh.nCellsI(); ++i)
            {
                BOOST_REQUIRE_EQUAL(mesh.iMap(i, j, k, 0), k*mesh.nCellsJ()*mesh.nCellsI() + j*mesh.nCellsI() + i + Parallel::processNo()*cellsPerProc);
            }
        }
    }

    for(int faceNo = 0; faceNo < 6; ++faceNo)
    {
        int adjProcNo = (*mesh.getAdjProcNoPtr())[faceNo];

        if(adjProcNo != Parallel::PROC_NULL)
        {
            const Array3D<int> &adjGlobalIndices = mesh.iMap.getAdjGlobalIndices(faceNo);

            for(int k = 0; k < mesh.nCellsK(); ++k)
            {
                for(int j = 0; j < mesh.nCellsJ(); ++j)
                {
                    for(int i = 0; i < mesh.nCellsI(); ++i)
                    {
                        BOOST_REQUIRE_EQUAL(adjGlobalIndices(i, j, k), k*mesh.nCellsJ()*mesh.nCellsI() + j*mesh.nCellsI() + i + adjProcNo*cellsPerProc);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE (test3)
{
    ParallelHexaFvmMesh mesh;
    mesh.changeName("indexMapTest3Mesh");
    mesh.initializeCartesianMesh(1, 1, 1, 20, 20, 20);
    SparseMatrix A;

    A.allocate(mesh.iMap.nActiveGlobal(), mesh.iMap.nActiveGlobal(), 7);

    BOOST_REQUIRE_EQUAL(A.iLower(), mesh.iMap(0, 0, 0, 0));
    BOOST_REQUIRE_EQUAL(A.iUpper(), mesh.iMap(mesh.uCellI(), mesh.uCellJ(), mesh.uCellK(), 0) + 1);
}

BOOST_AUTO_TEST_CASE (end)
{
    Parallel::finalize();
}

BOOST_AUTO_TEST_SUITE_END( )
