#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

#include "MultiBlockHexaFvmMesh.h"
#include "Parallel.h"

MultiBlockHexaFvmMesh::MultiBlockHexaFvmMesh()
    :
      blocks_(Parallel::nProcesses())
{

}

void MultiBlockHexaFvmMesh::initialize()
{
    using namespace boost::property_tree;
    using namespace std;

    ptree connectivity;

    read_info("mesh/multiBlockStructuredMesh/multiBlockStructuredMeshConnectivity.info", connectivity);

    if(Parallel::nProcesses() != connectivity.get<int>("numberOfBlocks"))
        Output::raiseException("MultiBlockHexaFvmMesh", "initialize", "number of blocks is not equal to the number of processes.");

    Parallel::barrier();

    blocks_[Parallel::processNo()].initialize("mesh/multiBlockStructuredMesh/structuredMesh_block" + std::to_string(Parallel::processNo()) + ".dat");

    adjacentBlockNo_[0] = connectivity.get<int>("Block" + std::to_string(Parallel::processNo()) + ".eastBlockNo");
    adjacentBlockNo_[1] = connectivity.get<int>("Block" + std::to_string(Parallel::processNo()) + ".westBlockNo");
    adjacentBlockNo_[2] = connectivity.get<int>("Block" + std::to_string(Parallel::processNo()) + ".northBlockNo");
    adjacentBlockNo_[3] = connectivity.get<int>("Block" + std::to_string(Parallel::processNo()) + ".southBlockNo");
    adjacentBlockNo_[4] = connectivity.get<int>("Block" + std::to_string(Parallel::processNo()) + ".topBlockNo");
    adjacentBlockNo_[5] = connectivity.get<int>("Block" + std::to_string(Parallel::processNo()) + ".bottomBlockNo");

    for(int i = 0; i < 6; ++i)
    {
        if(adjacentBlockNo_[i] < 0)
            continue;

        blocks_[adjacentBlockNo_[i]].initialize("mesh/multiBlockStructuredMesh/structuredMesh_block" + std::to_string(adjacentBlockNo_[i]) + ".dat");

        blocks_[Parallel::processNo()].addBoundaryMesh(blocks_[adjacentBlockNo_[i]], (HexaFvmMesh::Direction)i);
    }

    blocks_[Parallel::processNo()].iMap.generateIndices();
}

const HexaFvmMesh& MultiBlockHexaFvmMesh::operator ()() const
{
    return blocks_[Parallel::processNo()];
}

void MultiBlockHexaFvmMesh::writeTec360(double time, const std::string &directoryName)
{
    blocks_[Parallel::processNo()].writeTec360(time, directoryName);
}
