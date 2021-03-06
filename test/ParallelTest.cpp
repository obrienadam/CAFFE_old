#define BOOST_TEST_MODULE ParallelTest
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include "Parallel.h"
#include "Output.h"

BOOST_AUTO_TEST_SUITE (ParallelTest)

BOOST_AUTO_TEST_CASE (test1)
{
    Parallel::initialize();

    std::pair<int, int> range = Parallel::ownershipRange(1000);

    if(Parallel::processNo() < 1000%Parallel::nProcesses())
        BOOST_REQUIRE_EQUAL(1000/Parallel::nProcesses(), range.second - range.first);
    else
        BOOST_REQUIRE_EQUAL(1000/Parallel::nProcesses() - 1, range.second - range.first);
}

BOOST_AUTO_TEST_CASE (test2)
{
    Vector3D vec;

    if(Parallel::isMainProcessor())
        vec = Vector3D(1, 2, 3);

    vec = Parallel::broadcast(vec, Parallel::mainProcNo());

    BOOST_REQUIRE_EQUAL(vec.x, 1);
    BOOST_REQUIRE_EQUAL(vec.y, 2);
    BOOST_REQUIRE_EQUAL(vec.z, 3);

    Array3D<Vector3D> vecs(20, 20, 20);

    if(Parallel::isMainProcessor())
    {
        for(Vector3D &vec : vecs)
            vec = Vector3D(1, 2, 3);
    }

    Parallel::broadcast(Parallel::mainProcNo(), vecs);

    for(Vector3D &vec : vecs)
    {
        BOOST_REQUIRE_EQUAL(vec.x, 1);
        BOOST_REQUIRE_EQUAL(vec.y, 2);
        BOOST_REQUIRE_EQUAL(vec.z, 3);
    }

    std::string str;

    if(Parallel::isMainProcessor())
        str = "Broadcasted from proc 0!";

    str = Parallel::broadcast(str, 0);

    BOOST_REQUIRE(str == "Broadcasted from proc 0!");
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
            BOOST_REQUIRE_EQUAL(doubles[i], i);
            BOOST_REQUIRE_EQUAL(vecs[i].x, i);
            BOOST_REQUIRE_EQUAL(vecs[i].y, i + 1);
            BOOST_REQUIRE_EQUAL(vecs[i].z, i + 2);
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
            BOOST_REQUIRE_EQUAL(doubles[i], i);

        for(int i = 0; i < vecs.size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(vecs[i].x, i);
            BOOST_REQUIRE_EQUAL(vecs[i].y, i + 1);
            BOOST_REQUIRE_EQUAL(vecs[i].z, i + 2);
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
        BOOST_REQUIRE_EQUAL(allTestInts[i], i);
        BOOST_REQUIRE_EQUAL(allTestDoubles[i], i*M_PI);
    }
}

BOOST_AUTO_TEST_CASE (test6)
{
    Array3D<int> sourceArray(10, 10, 10), destArray;
    int i = 0;

    if(Parallel::isMainProcessor())
    {
        for(auto it = sourceArray.begin(); it != sourceArray.end(); ++it, ++i)
        {
            *it = i;
        }
    }
    Parallel::send(Parallel::mainProcNo(), Parallel::nProcesses() - 1, sourceArray, destArray);

    if(Parallel::processNo() == Parallel::processNo() - 1)
    {
        BOOST_REQUIRE_EQUAL(sourceArray.size(), destArray.size());

        i = 0;
        for(auto it = destArray.begin(); it != destArray.end(); ++it, ++i)
        {
            BOOST_REQUIRE_EQUAL(*it, i);
        }
    }
}

BOOST_AUTO_TEST_CASE (test7)
{
    Array3D<int> sourceArray(20, 20, 20), clockWise(20, 20, 20), cclockWise(20, 20, 20);
    std::vector<MPI::Request> requests;

    for(auto it = sourceArray.begin(); it != sourceArray.end(); ++it)
        *it = rand()%20;

    // Send information in two directions around a loop... meant to just test non-blocking communications
    if(Parallel::processNo() == 0)
    {
        requests.push_back(MPI::COMM_WORLD.Isend(sourceArray.data(), sourceArray.size(), MPI::INT, Parallel::processNo() + 1, 0));
        requests.push_back(MPI::COMM_WORLD.Isend(sourceArray.data(), sourceArray.size(), MPI::INT, Parallel::nProcesses() - 1, 1));
        requests.push_back(MPI::COMM_WORLD.Irecv(cclockWise.data(), cclockWise.size(), MPI::INT, Parallel::nProcesses() - 1, 0));
        requests.push_back(MPI::COMM_WORLD.Irecv(clockWise.data(), clockWise.size(), MPI::INT, Parallel::processNo() + 1, 1));
    }
    else if(Parallel::processNo() == Parallel::nProcesses() - 1)
    {
        requests.push_back(MPI::COMM_WORLD.Isend(sourceArray.data(), sourceArray.size(), MPI::INT, 0, 0));
        requests.push_back(MPI::COMM_WORLD.Isend(sourceArray.data(), sourceArray.size(), MPI::INT, Parallel::processNo() - 1, 1));
        requests.push_back(MPI::COMM_WORLD.Irecv(cclockWise.data(), cclockWise.size(), MPI::INT, Parallel::processNo() - 1, 0));
        requests.push_back(MPI::COMM_WORLD.Irecv(clockWise.data(), clockWise.size(), MPI::INT, 0, 1));
    }
    else
    {
        requests.push_back(MPI::COMM_WORLD.Isend(sourceArray.data(), sourceArray.size(), MPI::INT, Parallel::processNo() + 1, 0));
        requests.push_back(MPI::COMM_WORLD.Isend(sourceArray.data(), sourceArray.size(), MPI::INT, Parallel::processNo() - 1, 1));
        requests.push_back(MPI::COMM_WORLD.Irecv(cclockWise.data(), cclockWise.size(), MPI::INT, Parallel::processNo() - 1, 0));
        requests.push_back(MPI::COMM_WORLD.Irecv(clockWise.data(), clockWise.size(), MPI::INT, Parallel::processNo() + 1, 1));
    }

    MPI::Request::Waitall(requests.size(), requests.data());
}

BOOST_AUTO_TEST_CASE (test8)
{
    Array3D<Vector3D> sourceArray(20, 20, 20), destArray(20, 20, 20);

    for(Vector3D &vec3D : sourceArray)
        vec3D = Vector3D(Parallel::processNo(), 2*Parallel::processNo(), 3*Parallel::processNo());

    for(int i = 0; i < 3; ++i)
    {
        Output::print("Sending messages...");
        Parallel::iSend(0, 1, 0, sourceArray);
        Parallel::iSend(1, 2, 0, sourceArray);
        Parallel::iSend(1, 3, 0, sourceArray);
        Parallel::iSend(1, 0, 1, sourceArray);

        Output::print("Receiving messages...");
        Parallel::iRecv(0, 1, 0, destArray);
        Parallel::iRecv(1, 0, 1, destArray);
        Parallel::iRecv(1, 2, 0, destArray);
        Parallel::iRecv(1, 3, 0, destArray);

        Output::print("Finalizing messages...");
        Parallel::waitAll();
    }
    Output::print("Finished communication.");

    if(Parallel::processNo() == 0 || Parallel::processNo() == 2 || Parallel::processNo() == 3)
    {
        for(const Vector3D &vec : destArray)
        {
            BOOST_REQUIRE_EQUAL(vec.x, 1.);
            BOOST_REQUIRE_EQUAL(vec.y, 2.);
            BOOST_REQUIRE_EQUAL(vec.z, 3.);
        }
    }
}

BOOST_AUTO_TEST_CASE (end)
{
    Parallel::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
