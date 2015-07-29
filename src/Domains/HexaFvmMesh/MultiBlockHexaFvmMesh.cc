#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

#include "MultiBlockHexaFvmMesh.h"
#include "Parallel.h"

MultiBlockHexaFvmMesh::MultiBlockHexaFvmMesh()
{

}

void MultiBlockHexaFvmMesh::initialize()
{
    using namespace std;

    int blockNo = Parallel::processNo();

    nProcesses_ = Parallel::nProcesses();
    blocks_.resize(nProcesses_);

    Parallel::barrier();

    blocks_[blockNo].initialize("mesh/multiBlockStructuredMesh/structuredMesh_block" + std::to_string(blockNo) + ".dat");
    blocks_[blockNo].writeTec360(0, "solution");
}
