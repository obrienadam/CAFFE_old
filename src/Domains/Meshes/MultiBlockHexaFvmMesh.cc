#include <fstream>

#include "MultiBlockHexaFvmMesh.h"
#include "InputStringProcessing.h"

MultiBlockHexaFvmMesh::MultiBlockHexaFvmMesh()
{

}

void MultiBlockHexaFvmMesh::initialize(Input &input)
{
    using namespace std;

    int i, j, k, l, nI, nJ, nK, processNo, nBlocks, lNodeI, uNodeI, lNodeJ, nIThisProc;
    string buffer, name;
    ifstream fin;
    vector<string> bufferVec;
    Array3D<Point3D> nodes;
    double tmp;

    nBlocks = Parallel::nProcesses();
    allocate(nBlocks);

    fin.open("mesh/structuredMesh.dat");

    //- All processes read the header
    StructuredMesh::readTecplotMeshHeader(fin, name, nI, nJ, nK);

    Parallel::ownershipRange(nI, lNodeI, uNodeI, nIThisProc);

    nodes.allocate(nIThisProc, nJ, nK);

    for(l = 0; l < 3; ++l)
    {
        for(k = 0; k < nK; ++k)
        {
            for(j = 0; j < nJ; ++j)
            {
                for(i = 0; i < nI; ++i)
                {
                    fin >> tmp;

                    if(i >= lNodeI && i <= uNodeI + 1)
                    {
                        nodes(i - lNodeI, j, k)(l) = tmp;
                    }
                }
            }
        }
    }

    blocks_[Parallel::processNo()].addVectorField("u", CONSERVED);
    blocks_[Parallel::processNo()].addScalarField("phi", CONSERVED);
    blocks_[Parallel::processNo()].initialize(nodes);
    blocks_[Parallel::processNo()].name += "_block" + to_string(Parallel::processNo());
}

void MultiBlockHexaFvmMesh::allocate(int nBlocks)
{
    blocks_.resize(nBlocks);
    Parallel::ownershipRange(nBlocks, lBlockNo_, uBlockNo_, nBlocksThisProc_);
}

void MultiBlockHexaFvmMesh::writeTec360(double time, std::string directoryName)
{
    blocks_[Parallel::processNo()].writeTec360(time, directoryName);
}
