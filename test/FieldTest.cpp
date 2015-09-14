#define BOOST_TEST_MODULE FieldTest
#include <boost/test/included/unit_test.hpp>

#include "Field.h"
#include "PrimitiveBoundaryCondition.h"

BOOST_AUTO_TEST_SUITE (FieldTest)

BOOST_AUTO_TEST_CASE (test1)
{
    HexaFvmMesh mesh;
    mesh.initializeCartesianMesh(1, 1, 1, 10, 10, 10);

    Field<double> testField(mesh, Field<double>::CONSERVED, "testField");
    Field<double> testField2(testField);

    testField.setAll(0.8);
    testField2 = testField;

    for(int i = 0; i < testField.size(); ++i)
    {
        BOOST_REQUIRE_EQUAL(testField(i), 0.8);
        BOOST_REQUIRE_EQUAL(testField2(i), 0.8);
    }
}

BOOST_AUTO_TEST_CASE (test2)
{
    HexaFvmMesh mesh;
    mesh.initializeCartesianMesh(1, 1, 1, 50, 50, 50);

    Field<double> testField(mesh, Field<double>::CONSERVED, "testField");
    Array3D<double> subArray(20, 40, 10);

    testField.setValues(13, 32, 10, 49, 20, 29, 1);
    testField.getSubfield(13, 32, 10, 49, 20, 29, subArray);

    for(int k = 0; k < testField.sizeK(); ++k)
    {
        for(int j = 0; j < testField.sizeJ(); ++j)
        {
            for(int i = 0; i < testField.sizeI(); ++i)
            {
                if(i >= 13 && i < 33
                        && j >= 10 && j < 50
                        && k >= 20 && k < 30)
                {
                    BOOST_REQUIRE_EQUAL(subArray(i - 13, j - 10, k - 20), 1);
                    BOOST_REQUIRE_EQUAL(testField(i, j, k), subArray(i - 13, j - 10, k - 20));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE (test3)
{
    HexaFvmMesh mesh;
    double patchNo[6] = {1, 2, 3, 4, 5, 6};

    mesh.initializeCartesianMesh(1, 1, 1, 20, 20, 20);
    Field<double> testField(mesh, Field<double>::CONSERVED, "testField");
    testField.setFixedBoundaryPatches(patchNo);

    testField.setEastFacesFromPatch();
    testField.setWestFacesFromPatch();
    for(int k = 0; k < testField.sizeK(); ++k)
    {
        for(int j = 0; j < testField.sizeJ(); ++j)
        {
            BOOST_REQUIRE_EQUAL(testField(20, j, k), 1);
            BOOST_REQUIRE_EQUAL(testField(-1, j, k), 2);

            BOOST_REQUIRE_EQUAL(testField.faceE(19, j, k), 1);
            BOOST_REQUIRE_EQUAL(testField.faceW(0, j, k), 2);
        }
    }

    testField.setNorthFacesFromPatch();
    testField.setSouthFacesFromPatch();
    for(int k = 0; k < testField.sizeK(); ++k)
    {
        for(int i = 0; i < testField.sizeI(); ++i)
        {
            BOOST_REQUIRE_EQUAL(testField(i, 20, k), 3);
            BOOST_REQUIRE_EQUAL(testField(i, -1, k), 4);

            BOOST_REQUIRE_EQUAL(testField.faceN(i, 19, k), 3);
            BOOST_REQUIRE_EQUAL(testField.faceS(i, 0, k), 4);
        }
    }

    testField.setTopFacesFromPatch();
    testField.setBottomFacesFromPatch();
    for(int j = 0; j < testField.sizeJ(); ++j)
    {
        for(int i = 0; i < testField.sizeI(); ++i)
        {
            BOOST_REQUIRE_EQUAL(testField(i, j, 20), 5);
            BOOST_REQUIRE_EQUAL(testField(i, j, -1), 6);

            BOOST_REQUIRE_EQUAL(testField.faceT(i, j, 19), 5);
            BOOST_REQUIRE_EQUAL(testField.faceB(i, j, 0), 6);
        }
    }
}

BOOST_AUTO_TEST_CASE(test4)
{
    HexaFvmMesh mesh;

    mesh.initializeCartesianMesh(1, 1, 1, 10, 10, 10);

    Field<double> testField(mesh, Field<double>::CONSERVED, "testField");
    PrimitiveBoundaryCondition<double> bcs(testField);

    for(int i = 0; i < 6; ++i)
        bcs.changeType(i, PrimitiveBoundaryCondition<double>::FIXED, i + 1);

    for(int k = 0; k < testField.sizeK(); ++k)
    {
        for(int j = 0; j < testField.sizeJ(); ++j)
        {
            BOOST_REQUIRE_EQUAL(testField(10, j, k), 1);
            BOOST_REQUIRE_EQUAL(testField(-1, j, k), 2);

            BOOST_REQUIRE_EQUAL(testField.faceE(9, j, k), 1);
            BOOST_REQUIRE_EQUAL(testField.faceW(0, j, k), 2);
        }
    }

    for(int k = 0; k < testField.sizeK(); ++k)
    {
        for(int i = 0; i < testField.sizeI(); ++i)
        {
            BOOST_REQUIRE_EQUAL(testField(i, 10, k), 3);
            BOOST_REQUIRE_EQUAL(testField(i, -1, k), 4);

            BOOST_REQUIRE_EQUAL(testField.faceN(i, 9, k), 3);
            BOOST_REQUIRE_EQUAL(testField.faceS(i, 0, k), 4);
        }
    }

    for(int j = 0; j < testField.sizeJ(); ++j)
    {
        for(int i = 0; i < testField.sizeI(); ++i)
        {
            BOOST_REQUIRE_EQUAL(testField(i, j, 10), 5);
            BOOST_REQUIRE_EQUAL(testField(i, j, -1), 6);

            BOOST_REQUIRE_EQUAL(testField.faceT(i, j, 9), 5);
            BOOST_REQUIRE_EQUAL(testField.faceB(j, j, 0), 6);
        }
    }
}

BOOST_AUTO_TEST_CASE(test5)
{
    HexaFvmMesh mesh;

    mesh.initializeCartesianMesh(1, 1, 1, 10, 10, 10);

    Field<double> testField(mesh, Field<double>::CONSERVED, "testField");
    PrimitiveBoundaryCondition<double> bcs(testField);

    for(int i = 0; i < 6; ++i)
        bcs.changeType(i, PrimitiveBoundaryCondition<double>::ZERO_GRADIENT, 0);

    testField.setValues(9, 9, 0, 9, 0, 9, 1);
    testField.setValues(0, 0, 0, 9, 0, 9, 1);
    testField.setValues(0, 9, 9, 9, 0, 9, 1);
    testField.setValues(0, 9, 0, 0, 0, 9, 1);
    testField.setValues(0, 9, 0, 9, 9, 9, 1);
    testField.setValues(0, 9, 0, 9, 0, 0, 1);

    bcs.setBoundaries();

    for(int k = 0; k < testField.sizeK(); ++k)
    {
        for(int j = 0; j < testField.sizeJ(); ++j)
        {
            BOOST_REQUIRE_EQUAL(testField(10, j, k), 1);
            BOOST_REQUIRE_EQUAL(testField(-1, j, k), 1);

            BOOST_REQUIRE_EQUAL(testField.faceE(9, j, k), 1);
            BOOST_REQUIRE_EQUAL(testField.faceW(0, j, k), 1);
        }
    }

    for(int k = 0; k < testField.sizeK(); ++k)
    {
        for(int i = 0; i < testField.sizeI(); ++i)
        {
            BOOST_REQUIRE_EQUAL(testField(i, 10, k), 1);
            BOOST_REQUIRE_EQUAL(testField(i, -1, k), 1);

            BOOST_REQUIRE_EQUAL(testField.faceN(i, 9, k), 1);
            BOOST_REQUIRE_EQUAL(testField.faceS(i, 0, k), 1);
        }
    }

    for(int j = 0; j < testField.sizeJ(); ++j)
    {
        for(int i = 0; i < testField.sizeI(); ++i)
        {
            BOOST_REQUIRE_EQUAL(testField(i, j, 10), 1);
            BOOST_REQUIRE_EQUAL(testField(i, j, -1), 1);

            BOOST_REQUIRE_EQUAL(testField.faceT(i, j, 9), 1);
            BOOST_REQUIRE_EQUAL(testField.faceB(j, j, 0), 1);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )
