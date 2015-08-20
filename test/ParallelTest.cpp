#define BOOST_TEST_MODULE ParallelTest
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include "Parallel.h"

BOOST_AUTO_TEST_SUITE (ParallelTest)

BOOST_AUTO_TEST_CASE (test1)
{
    Parallel::initialize();

    std::pair<int, int> range = Parallel::ownershipRange(1000);

    if(Parallel::processNo() < 1000%Parallel::nProcesses())
        BOOST_CHECK_EQUAL(1000/Parallel::nProcesses(), range.second - range.first);
    else
        BOOST_CHECK_EQUAL(1000/Parallel::nProcesses() - 1, range.second - range.first);
}

BOOST_AUTO_TEST_CASE (test2)
{
    Vector3D vec;

    if(Parallel::isMainProcessor())
        vec = Vector3D(1, 2, 3);

    vec = Parallel::broadcast(vec, Parallel::mainProcNo());

    BOOST_CHECK_EQUAL(vec.x, 1);
    BOOST_CHECK_EQUAL(vec.y, 2);
    BOOST_CHECK_EQUAL(vec.z, 3);
}

BOOST_AUTO_TEST_CASE (test3)
{
    std::vector<double> doubles(1000);
    std::vector<Vector3D> vecs(1000);

    if(Parallel::isMainProcessor())
    {
        for(int i = 0; i < 1000; ++i)
        {
            doubles[i] = i;
            vecs[i] = Vector3D(i, i + 1, i + 2);
        }
    }

    Parallel::send(Parallel::mainProcNo(), Parallel::nProcesses() - 1, doubles);
    Parallel::send(Parallel::mainProcNo(), Parallel::nProcesses() - 1, vecs);

    if(Parallel::processNo() == Parallel::nProcesses() - 1)
    {
        for(int i = 0; i < 1000; ++i)

        {
            BOOST_CHECK_EQUAL(doubles[i], i);
            BOOST_CHECK_EQUAL(vecs[i].x, i);
            BOOST_CHECK_EQUAL(vecs[i].y, i + 1);
            BOOST_CHECK_EQUAL(vecs[i].z, i + 2);
        }
    }
}

BOOST_AUTO_TEST_CASE (test4)
{
    Array3D<double> doubles(10, 20, 12);
    Array3D<Vector3D> vecs(10, 20, 12);

    if(Parallel::isMainProcessor())
    {
        for(int i = 0; i < doubles.size(); ++i)
            doubles[i] = i;

        for(int i = 0; i < vecs.size(); ++i)
            vecs[i] = Vector3D(i, i + 1, i + 2);
    }

    Parallel::send(Parallel::mainProcNo(), Parallel::nProcesses() - 1, doubles);
    Parallel::send(Parallel::mainProcNo(), Parallel::nProcesses() - 1, vecs);

    if(Parallel::processNo() == Parallel::nProcesses() - 1)
    {
        for(int i = 0; i < doubles.size(); ++i)
            BOOST_CHECK_EQUAL(doubles[i], i);

        for(int i = 0; i < vecs.size(); ++i)
        {
            BOOST_CHECK_EQUAL(vecs[i].x, i);
            BOOST_CHECK_EQUAL(vecs[i].y, i + 1);
            BOOST_CHECK_EQUAL(vecs[i].z, i + 2);
        }
    }
}

BOOST_AUTO_TEST_CASE (test5)
{
    int testInt = Parallel::processNo();
    double testDouble = Parallel::processNo()*M_PI;
    std::vector<int> allTestInts(Parallel::nProcesses());
    std::vector<double> allTestDoubles(Parallel::nProcesses());

    Parallel::allGather(testInt, allTestInts);
    Parallel::allGather(testDouble, allTestDoubles);

    for(int i = 0; i < Parallel::nProcesses(); ++i)
    {
        BOOST_CHECK_EQUAL(allTestInts[i], i);
        BOOST_CHECK_EQUAL(allTestDoubles[i], i*M_PI);
    }
}

BOOST_AUTO_TEST_CASE (test6)
{
    Parallel::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
