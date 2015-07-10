#ifndef MULTI_BLOCK_HEXA_FVM_MESH_H
#define MULTI_BLOCK_HEXA_FVM_MESH_H

#include <vector>
#include <string>

#include "Input.h"
#include "HexaFvmMesh.h"
#include "Parallel.h"

class MultiBlockHexaFvmMesh
{

private:

    int lBlockNo_, uBlockNo_, nBlocksThisProc_;
    std::vector<HexaFvmMesh> blocks_;

public:

    MultiBlockHexaFvmMesh();

    std::string name;

    void initialize(Input& input);
    void allocate(int nBlocks);

    void writeTec360(double time = 0, std::string directoryName = "");
};

#endif
