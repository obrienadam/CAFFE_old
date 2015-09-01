#define BOOST_TEST_MODULE ParallelHexaFvmMeshTest
#include <boost/test/included/unit_test.hpp>

#include "ParallelHexaFvmMesh.h"
#include "Field.h"
#include "PrimitiveBoundaryCondition.h"

BOOST_AUTO_TEST_SUITE (ParallelHexaFvmMeshTest)

BOOST_AUTO_TEST_CASE (test1)
{
    Parallel::initialize();

    ParallelHexaFvmMesh mesh;
    mesh.initializeCartesianMesh(1, 1, 1, 64, 64, 64);
    mesh.changeName("test1Mesh");

    BOOST_CHECK_EQUAL(mesh.nCells(), 64*64*64/8);

    //- Output a debug file
    mesh.writeTec360(0, "");
    mesh.writeBoundaryMeshes(0, "");
}

BOOST_AUTO_TEST_CASE (test2)
{
    ParallelHexaFvmMesh mesh;

    mesh.changeName("test2Mesh");
    mesh.initializeCartesianMesh(1, 1, 1, 112, 112, 112);
}

BOOST_AUTO_TEST_CASE (test3)
{
    ParallelHexaFvmMesh mesh;
    mesh.changeName("test3Mesh");
    mesh.initializeCartesianMesh(1, 1.5, 2, 30, 40, 60);

    Field<double> phi(mesh, Field<double>::CONSERVED, "procNo");
    phi.setAll(Parallel::processNo());

    mesh.addArray3DToTecplotOutput("phi", phi.cellData());
    mesh.writeTec360(0, "");

    PrimitiveBoundaryCondition<double> bcs(phi);
    bcs.setParallelBoundaries(mesh.getAdjProcNoPtr());

    for(int i = 0; i < 6; ++i)
    {
        if((*mesh.getAdjProcNoPtr())[i] != Parallel::PROC_NULL)
            BOOST_CHECK_EQUAL(bcs.getType(i), PrimitiveBoundaryCondition<double>::PARALLEL);
        else
            BOOST_CHECK_EQUAL(bcs.getType(i), PrimitiveBoundaryCondition<double>::FIXED);
    }

    bcs.setBoundaries();

    for(int k = 0; k < phi.sizeK(); ++k)
    {
        for(int j = 0; j < phi.sizeJ(); ++j)
        {
            if((*mesh.getAdjProcNoPtr())[0] != Parallel::PROC_NULL)
                BOOST_CHECK_EQUAL(phi(phi.sizeI(), j, k), (*mesh.getAdjProcNoPtr())[0]);
            if((*mesh.getAdjProcNoPtr())[1] != Parallel::PROC_NULL)
                BOOST_CHECK_EQUAL(phi(-1, j, k), (*mesh.getAdjProcNoPtr())[1]);
        }
    }

    for(int k = 0; k < phi.sizeK(); ++k)
    {
        for(int i = 0; i < phi.sizeI(); ++i)
        {
            if((*mesh.getAdjProcNoPtr())[2] != Parallel::PROC_NULL)
                BOOST_CHECK_EQUAL(phi(i, phi.sizeJ(), k), (*mesh.getAdjProcNoPtr())[2]);
            if((*mesh.getAdjProcNoPtr())[3] != Parallel::PROC_NULL)
                BOOST_CHECK_EQUAL(phi(i, -1, k), (*mesh.getAdjProcNoPtr())[3]);
        }
    }

    for(int j = 0; j < phi.sizeJ(); ++j)
    {
        for(int i = 0; i < phi.sizeI(); ++i)
        {
            if((*mesh.getAdjProcNoPtr())[4] != Parallel::PROC_NULL)
                BOOST_CHECK_EQUAL(phi(i, j, phi.sizeK()), (*mesh.getAdjProcNoPtr())[4]);
            if((*mesh.getAdjProcNoPtr())[5] != Parallel::PROC_NULL)
                BOOST_CHECK_EQUAL(phi(i, j, -1), (*mesh.getAdjProcNoPtr())[5]);
        }
    }
}

BOOST_AUTO_TEST_CASE (end)
{
    Parallel::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
