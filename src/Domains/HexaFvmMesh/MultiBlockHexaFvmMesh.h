#ifndef MULTI_BLOCK_HEXA_FVM_MESH_H
#define MULTI_BLOCK_HEXA_FVM_MESH_H

#include <vector>
#include <string>

#include "Input.h"
#include "HexaFvmMesh.h"

class MultiBlockHexaFvmMesh
{
public:

    MultiBlockHexaFvmMesh();

    void initialize();

    int nBlocks() const { return blocks_.size(); }

private:

    std::vector<HexaFvmMesh> blocks_;

    int nProcesses_, lProcess_, uProcess_;
};

#endif
