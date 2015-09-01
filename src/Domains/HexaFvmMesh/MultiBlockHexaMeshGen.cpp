#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

#include "MultiBlockHexaMeshGen.h"
#include "Output.h"

MultiBlockHexaMeshGen::MultiBlockHexaMeshGen()
{

}

void MultiBlockHexaMeshGen::readFile(const std::string &directory)
{
    using namespace boost::property_tree;

    ptree multiBlockMeshParameters;

    mesh.initialize(directory + "structuredMesh.dat");
    read_info(directory + "multiBlockStructuredMesh.info", multiBlockMeshParameters);

    nBlocksI_ = multiBlockMeshParameters.get<int>("Decomposition.nBlocksI");
    nBlocksJ_ = multiBlockMeshParameters.get<int>("Decomposition.nBlocksJ");
    nBlocksK_ = multiBlockMeshParameters.get<int>("Decomposition.nBlocksK");

    Output::print("MultiBlockHexaMeshGen", "Successfully read file \"mesh/multiBlockStructuredMesh.info\".");
}

void MultiBlockHexaMeshGen::writeMeshFiles(const std::string &directory)
{
    using namespace boost::filesystem;

    int i, j, k, ii, jj, kk, subMeshNo, lNodeI, uNodeI, lNodeJ, uNodeJ, lNodeK, uNodeK, nLocalNodesI, nLocalNodesJ, nLocalNodesK;
    int eastBlockNo, westBlockNo, northBlockNo, southBlockNo, topBlockNo, bottomBlockNo;
    Array3D<Point3D> localNodes;
    StructuredMesh subMesh;
    std::ofstream fout;

    Output::print("MultiBlockHexaMeshGen", "Creating blocks...");

    create_directory(directory + "multiBlockStructuredMesh");
    fout.open(directory + "multiBlockStructuredMesh/multiBlockStructuredMeshConnectivity.info");
    fout << "; block connectivity file for the multi-block structured mesh" << std::endl
         << std::endl
         << "numberOfBlocks " << nBlocksI_*nBlocksJ_*nBlocksK_ << std::endl;

    subMeshNo = 0;
    for(k = 0; k < nBlocksK_; ++k)
    {
        for(j = 0; j < nBlocksJ_; ++j)
        {
            for(i = 0; i < nBlocksI_; ++i)
            {
                getOwnershipRange(mesh.nNodesI(), nBlocksI_, i, lNodeI, uNodeI, nLocalNodesI);
                getOwnershipRange(mesh.nNodesJ(), nBlocksJ_, j, lNodeJ, uNodeJ, nLocalNodesJ);
                getOwnershipRange(mesh.nNodesK(), nBlocksK_, k, lNodeK, uNodeK, nLocalNodesK);

                //- Add additional nodes such that adjacent blocks have common nodes
                if(i != 0)
                {
                    --lNodeI;
                    ++nLocalNodesI;
                }
                if(j != 0)
                {
                    --lNodeJ;
                    ++nLocalNodesJ;
                }
                if(k != 0)
                {
                    --lNodeK;
                    ++nLocalNodesK;
                }

                localNodes.resize(nLocalNodesI, nLocalNodesJ, nLocalNodesK);

                for(kk = lNodeK; kk <= uNodeK; ++kk)
                {
                    for(jj = lNodeJ; jj <= uNodeJ; ++jj)
                    {
                        for(ii = lNodeI; ii <= uNodeI; ++ii)
                        {
                            localNodes(ii - lNodeI, jj - lNodeJ, kk - lNodeK) = mesh.node(ii, jj, kk);
                        }
                    }
                }

                subMesh.changeName("structured_block" + std::to_string(subMeshNo));
                subMesh.initialize(localNodes);
                subMesh.writeTec360(0, directory + "multiBlockStructuredMesh");
                subMesh.resetFileStream();

                if(i < nBlocksI_ - 1)
                    eastBlockNo = k*nBlocksI_*nBlocksJ_ + j*nBlocksI_ + i + 1;
                else
                    eastBlockNo = -1;

                if(i > 0)
                    westBlockNo = k*nBlocksI_*nBlocksJ_ + j*nBlocksI_ + i - 1;
                else
                    westBlockNo = -1;

                if(j < nBlocksJ_ - 1)
                    northBlockNo = k*nBlocksI_*nBlocksJ_ + (j + 1)*nBlocksI_ + i;
                else
                    northBlockNo = -1;

                if(j > 0)
                    southBlockNo = k*nBlocksI_*nBlocksJ_ + (j - 1)*nBlocksI_ + i;
                else
                    southBlockNo = -1;

                if(k < nBlocksK_ - 1)
                    topBlockNo = (k + 1)*nBlocksI_*nBlocksJ_ + j*nBlocksI_ + i;
                else
                    topBlockNo = -1;

                if(k > 0)
                    bottomBlockNo = (k - 1)*nBlocksI_*nBlocksJ_ + j*nBlocksI_ + i;
                else
                    bottomBlockNo = -1;

                fout << "Block" << subMeshNo << std::endl
                     << "{" << std::endl
                     << "   eastBlockNo " << eastBlockNo << std::endl
                     << "   westBlockNo " << westBlockNo << std::endl
                     << "   northBlockNo " << northBlockNo << std::endl
                     << "   southBlockNo " << southBlockNo << std::endl
                     << "   topBlockNo " << topBlockNo << std::endl
                     << "   bottomBlockNo " << bottomBlockNo << std::endl
                     << "}" << std::endl
                     << std::endl;

                ++subMeshNo;
            }
        }
    }

    fout.close();
    Output::print("MultiBlockHexaMeshGen", "Finished creating blocks.");
}

//******************************* Private methods ***************************************

void MultiBlockHexaMeshGen::getOwnershipRange(int nEntities, int nProcs, int procNo, int &iLower, int &iUpper, int &nEntitiesThisProc)
{
    int nRemainingEntities;

    nEntitiesThisProc = nEntities/nProcs;
    nRemainingEntities = nEntities%nProcs;

    if(procNo < nRemainingEntities)
    {
        ++nEntitiesThisProc;
        iLower = procNo*nEntitiesThisProc;
        iUpper = iLower + nEntitiesThisProc - 1;
    }
    else
    {
        iLower = (nEntitiesThisProc + 1)*nRemainingEntities + (procNo - nRemainingEntities)*nEntitiesThisProc;
        iUpper = iLower + nEntitiesThisProc - 1;
    }
}
