#include <boost/algorithm/string.hpp>

#include "ParallelHexaFvmMesh.h"
#include "Output.h"

ParallelHexaFvmMesh::ParallelHexaFvmMesh()
    :
      boundaryMeshes_(6)
{

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

void ParallelHexaFvmMesh::writeTec360(double time, const std::string &directory)
{
    // make sure the name contains the apropriate sub-string in the name, to ensure the file-streams don't end up being the same
    if(!(foutTec360_.is_open() || boost::contains(name, "subDomain" + std::to_string(Parallel::processNo()))))
        name += "subDomain" + std::to_string(Parallel::processNo());

    HexaFvmMesh::writeTec360(time, directory);
}

//************************* Private methods **********************************

void ParallelHexaFvmMesh::initializeSubDomains(const StructuredMesh &tmpMesh)
{
    int nCellsIGlobal, nCellsJGlobal, nCellsKGlobal;
    std::vector<int> nCellsILocal, nCellsJLocal, nCellsKLocal;
    int nSubDomainsI, nSubDomainsJ, nSubDomainsK, subDomainNo;
    std::vector<Point3D> vertices(8);

    name += "_subDomain" + std::to_string(Parallel::processNo());

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

                    for(int i = 0; i < 6; ++i)
                        adjacentSubDomainProcNo_[i] = -1;

                    if(i < nSubDomainsI - 1)
                        adjacentSubDomainProcNo_[EAST] = subDomainNo + 1;
                    if(i > 0)
                        adjacentSubDomainProcNo_[WEST] = subDomainNo - 1;
                    if(j < nSubDomainsJ - 1)
                        adjacentSubDomainProcNo_[NORTH] = subDomainNo + nSubDomainsI;
                    if(j > 0)
                        adjacentSubDomainProcNo_[SOUTH] = subDomainNo - nSubDomainsI;
                    if(k < nSubDomainsK - 1)
                        adjacentSubDomainProcNo_[TOP] = subDomainNo + nSubDomainsI*nSubDomainsJ;
                    if(k > 0)
                        adjacentSubDomainProcNo_[BOTTOM] = subDomainNo - nSubDomainsI*nSubDomainsJ;
                }
            }
        }
    }

    Output::print("ParallelHexaFvmMesh", "finished allocating sub-domains.");

    iMap.initialize(nCellsI(), nCellsJ(), nCellsK());
    iMap.generateLocalIndices();
    iMap.generateGlobalIndices(adjacentSubDomainProcNo_);
}
