#include <boost/algorithm/string.hpp>

#include "ParallelHexaFvmMesh.h"
#include "Output.h"

ParallelHexaFvmMesh::ParallelHexaFvmMesh()
    :
      adjProcNoPtr_(new std::array<int, 6>)
{
    adjProcNoPtr_->fill(Parallel::PROC_NULL);
}

void ParallelHexaFvmMesh::initialize(const std::string &filename)
{
    StructuredMesh tmpMesh;

    if(Parallel::isMainProcessor())
        tmpMesh.initialize(filename);

    initializeSubDomains(tmpMesh);
}

void ParallelHexaFvmMesh::initialize(const Array3D<Point3D> &nodes)
{
    StructuredMesh tmpMesh;

    if(Parallel::isMainProcessor())
        tmpMesh.initialize(nodes);

    initializeSubDomains(tmpMesh);
}

void ParallelHexaFvmMesh::initializeCartesianMesh(double xLength, double yLength, double zLength, int nCellsI, int nCellsJ, int nCellsK)
{
    StructuredMesh tmpMesh;

    if(Parallel::isMainProcessor())
        tmpMesh.initializeCartesianMesh(xLength, yLength, zLength, nCellsI, nCellsJ, nCellsK);

    initializeSubDomains(tmpMesh);
}

std::string ParallelHexaFvmMesh::meshStats()
{
    using namespace std;

    ostringstream sout;

    sout << "Nodes per proc in I direction -> " << nNodesI() << endl
         << "Nodes per proc in J direction -> " << nNodesJ() << endl
         << "Nodes per proc in K direction -> " << nNodesK() << endl
         << "Nodes per proc total -> " << nNodes() << endl
         << "Nodes global total -> " << Parallel::sum(nNodes()) << endl
         << "Cells per proc in I direction -> " << nCellsI() << endl
         << "Cells per proc in J direction -> " << nCellsJ() << endl
         << "Cells per proc in K direction -> " << nCellsK() << endl
         << "Cells per proc total -> " << nCells() << endl
         << "Cells global total -> " << Parallel::sum(nCells()) << endl;
    return sout.str();
}

void ParallelHexaFvmMesh::writeTec360(double time, const std::string &directory)
{
    HexaFvmMesh::writeTec360(time, directory);
}

void ParallelHexaFvmMesh::writeBoundaryMeshes(double time, const std::string &directory)
{
    using namespace std;

    Output::print("ParallelHexaFvmMesh", "Writing boundary meshes to Tec360 ASCII...");

    for(int i = 0; i < 6; ++i)
    {
        if(boundaryMeshPointer(static_cast<Direction>(i)))
            boundaryMeshPointer(static_cast<Direction>(i))->writeTec360(time, directory);
    }

    Output::print("ParallelHexaFvmMesh", "Finished wrting boundary meshes to Tec360 ASCII.");
}

void ParallelHexaFvmMesh::changeName(const std::string &newName)
{
    name_ = newName + "_subDomainNo" + std::to_string(Parallel::processNo());

    for(int i = 0; i < 6; ++i)
    {
        if(boundaryMeshPointer(static_cast<Direction>(i)))
            boundaryMeshPointer(static_cast<Direction>(i))->changeName(name_ + "_boundaryMesh" + std::to_string(i));
    }

    resetFileStream();
}

//************************* Private methods **********************************

void ParallelHexaFvmMesh::initializeSubDomains(const StructuredMesh &tmpMesh)
{
    using namespace std;

    int nCellsIGlobal, nCellsJGlobal, nCellsKGlobal;
    std::vector<int> nCellsILocal, nCellsJLocal, nCellsKLocal;
    int nSubDomainsI, nSubDomainsJ, nSubDomainsK, subDomainNo;
    vector<Point3D> vertices(8);
    Array3D<Point3D> boundaryNodesSend[6], boundaryNodesRecv[6];
    vector< vector<double> > sendBuffers(1000, vector<double>(1000)),
            recvBuffers(1000, vector<double>(1000));

    //- Determine domain decomposition
    Output::issueWarning("ParallelHexaFvmMesh", "initializeSubDomains", "number of sub-domains currently fixed. May not be suitable for all computations.");
    nSubDomainsI = 2;
    nSubDomainsJ = 2;
    nSubDomainsK = 2;

    //- Broadcast the global number of cells from the root process
    nCellsIGlobal = Parallel::broadcast(tmpMesh.nNodesI() + 1, Parallel::mainProcNo());
    nCellsJGlobal = Parallel::broadcast(tmpMesh.nNodesJ() + 1, Parallel::mainProcNo());
    nCellsKGlobal = Parallel::broadcast(tmpMesh.nNodesK() + 1, Parallel::mainProcNo());

    //- Create a vector to store the local number of cells on each process, so every process has access to this info
    nCellsILocal.resize(Parallel::nProcesses());
    nCellsJLocal.resize(Parallel::nProcesses());
    nCellsKLocal.resize(Parallel::nProcesses());
    Parallel::allGather(nCellsIGlobal/nSubDomainsI, nCellsILocal);
    Parallel::allGather(nCellsJGlobal/nSubDomainsJ, nCellsJLocal);
    Parallel::allGather(nCellsKGlobal/nSubDomainsK, nCellsKLocal);

    //- Construct a structured mesh where the vertices outline each sub domain. Broadcast this to all processes
    if(Parallel::isMainProcessor())
    {
        vertices[0] = tmpMesh.node(0, 0, 0);
        vertices[1] = tmpMesh.node(tmpMesh.uNodeI(), 0, 0);
        vertices[2] = tmpMesh.node(tmpMesh.uNodeI(), tmpMesh.uNodeJ(), 0);
        vertices[3] = tmpMesh.node(0, tmpMesh.uNodeJ(), 0);
        vertices[4] = tmpMesh.node(0, 0, tmpMesh.uNodeK());
        vertices[5] = tmpMesh.node(tmpMesh.uNodeI(), 0, tmpMesh.uNodeK());
        vertices[6] = tmpMesh.node(tmpMesh.uNodeI(), tmpMesh.uNodeJ(), tmpMesh.uNodeK());
        vertices[7] = tmpMesh.node(0, tmpMesh.uNodeJ(), tmpMesh.uNodeK());
    }
    Parallel::broadcast(Parallel::mainProcNo(), vertices);
    HexaFvmMesh::initialize(vertices, nSubDomainsI, nSubDomainsJ, nSubDomainsK);

    std::array<int, 6> &adjProcNo = *adjProcNoPtr_;

    //- Using the subdomain outline, allocate cells one each local domain, and create adjacency list
    for(int k = 0; k < nSubDomainsK; ++k)
    {
        for(int j = 0; j < nSubDomainsJ; ++j)
        {
            for(int i = 0; i < nSubDomainsI; ++i)
            {
                subDomainNo = k*nSubDomainsI*nSubDomainsJ + j*nSubDomainsI + i;

                if(subDomainNo == Parallel::processNo())
                {
                    vertices[0] = StructuredMesh::node(i, j, k);
                    vertices[1] = StructuredMesh::node(i + 1, j, k);
                    vertices[2] = StructuredMesh::node(i + 1, j + 1, k);
                    vertices[3] = StructuredMesh::node(i, j + 1, k);
                    vertices[4] = StructuredMesh::node(i, j, k + 1);
                    vertices[5] = StructuredMesh::node(i + 1, j, k + 1);
                    vertices[6] = StructuredMesh::node(i + 1, j + 1, k + 1);
                    vertices[7] = StructuredMesh::node(i, j + 1, k + 1);

                    HexaFvmMesh::initialize(vertices, nCellsILocal[subDomainNo], nCellsJLocal[subDomainNo], nCellsKLocal[subDomainNo]);

                    if(i < nSubDomainsI - 1)
                        adjProcNo[EAST] = subDomainNo + 1;
                    if(i > 0)
                        adjProcNo[WEST] = subDomainNo - 1;
                    if(j < nSubDomainsJ - 1)
                        adjProcNo[NORTH] = subDomainNo + nSubDomainsI;
                    if(j > 0)
                        adjProcNo[SOUTH] = subDomainNo - nSubDomainsI;
                    if(k < nSubDomainsK - 1)
                        adjProcNo[TOP] = subDomainNo + nSubDomainsI*nSubDomainsJ;
                    if(k > 0)
                        adjProcNo[BOTTOM] = subDomainNo - nSubDomainsI*nSubDomainsJ;
                }
            }
        }
    }

    //- With the mesh initialization complete, extract nodes around boundaries and send them to the apropriate process
    for(int i = 0; i < 6; ++i)
    {
        if(adjProcNo[i] != Parallel::PROC_NULL)
        {
            switch(i)
            {
            case EAST:
                boundaryNodesSend[0].resize(2, nNodesJ(), nNodesK());
                boundaryNodesRecv[0].resize(2, nNodesJ(), nNodesK());

                for(int k = 0; k < nNodesK(); ++k)
                    for(int j = 0; j < nNodesJ(); ++j)
                        for(int i = 0; i < 2; ++i)
                            boundaryNodesSend[0](i, j, k) = StructuredMesh::node(uNodeI() - 1 + i, j, k);

                Parallel::iSend(Parallel::processNo(), adjProcNo[EAST], 0, boundaryNodesSend[0]);
                Parallel::iRecv(adjProcNo[EAST], Parallel::processNo(), 1, boundaryNodesRecv[0]);
                break;

            case WEST:
                boundaryNodesSend[1].resize(2, nNodesJ(), nNodesK());
                boundaryNodesRecv[1].resize(2, nNodesJ(), nNodesK());

                for(int k = 0; k < nNodesK(); ++k)
                    for(int j = 0; j < nNodesJ(); ++j)
                        for(int i = 0; i < 2; ++i)
                            boundaryNodesSend[1](i, j, k) = StructuredMesh::node(i, j, k);

                Parallel::iSend(Parallel::processNo(), adjProcNo[WEST], 1, boundaryNodesSend[1]);
                Parallel::iRecv(adjProcNo[WEST], Parallel::processNo(), 0, boundaryNodesRecv[1]);
                break;

            case NORTH:
                boundaryNodesSend[2].resize(nNodesI(), 2, nNodesK());
                boundaryNodesRecv[2].resize(nNodesI(), 2, nNodesK());

                for(int k = 0; k < nNodesK(); ++k)
                    for(int j = 0; j < 2; ++j)
                        for(int i = 0; i < nNodesI(); ++i)
                            boundaryNodesSend[2](i, j, k) = StructuredMesh::node(i, uNodeJ() - 1 + j, k);

                Parallel::iSend(Parallel::processNo(), adjProcNo[NORTH], 2, boundaryNodesSend[2]);
                Parallel::iRecv(adjProcNo[NORTH], Parallel::processNo(), 3, boundaryNodesRecv[2]);
                break;

            case SOUTH:
                boundaryNodesSend[3].resize(nNodesI(), 2, nNodesK());
                boundaryNodesRecv[3].resize(nNodesI(), 2, nNodesK());

                for(int k = 0; k < nNodesK(); ++k)
                    for(int j = 0; j < 2; ++j)
                        for(int i = 0; i < nNodesI(); ++i)
                            boundaryNodesSend[3](i, j, k) = StructuredMesh::node(i, j, k);

                Parallel::iSend(Parallel::processNo(), adjProcNo[SOUTH], 3, boundaryNodesSend[3]);
                Parallel::iRecv(adjProcNo[SOUTH], Parallel::processNo(), 2, boundaryNodesRecv[3]);
                break;

            case TOP:
                boundaryNodesSend[4].resize(nNodesI(), nNodesJ(), 2);
                boundaryNodesRecv[4].resize(nNodesI(), nNodesJ(), 2);

                for(int k = 0; k < 2; ++k)
                    for(int j = 0; j < nNodesJ(); ++j)
                        for(int i = 0; i < nNodesI(); ++i)
                            boundaryNodesSend[4](i, j, k) = StructuredMesh::node(i, j, uNodeK() - 1 + k);

                Parallel::iSend(Parallel::processNo(), adjProcNo[TOP], 4, boundaryNodesSend[4]);
                Parallel::iRecv(adjProcNo[TOP], Parallel::processNo(), 5, boundaryNodesRecv[4]);
                break;

            case BOTTOM:
                boundaryNodesSend[5].resize(nNodesI(), nNodesJ(), 2);
                boundaryNodesRecv[5].resize(nNodesI(), nNodesJ(), 2);

                for(int k = 0; k < 2; ++k)
                    for(int j = 0; j < nNodesJ(); ++j)
                        for(int i = 0; i < nNodesI(); ++i)
                            boundaryNodesSend[5](i, j, k) = StructuredMesh::node(i, j, k);

                Parallel::iSend(Parallel::processNo(), adjProcNo[BOTTOM], 5, boundaryNodesSend[5]);
                Parallel::iRecv(adjProcNo[BOTTOM], Parallel::processNo(), 4, boundaryNodesRecv[5]);
                break;
            };
        }
    }
    Parallel::waitAll();

    for(int i = 0; i < 6; ++i)
    {
        if(adjProcNo[i] != Parallel::PROC_NULL)
        {
            std::shared_ptr<HexaFvmMesh> meshPtr(new HexaFvmMesh);

            meshPtr->initialize(boundaryNodesRecv[i]);
            addBoundaryMesh(meshPtr, static_cast<Direction>(i));
        }
    }

    //- Call this method to update the boundary mesh names
    changeName(name_);

    Output::print("ParallelHexaFvmMesh", "finished allocating sub-domains.");

    iMap.initialize(nCellsI(), nCellsJ(), nCellsK());
    iMap.generateLocalIndices();
    iMap.generateGlobalIndices(adjProcNoPtr_);
}
