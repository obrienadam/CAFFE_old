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

    const HexaFvmMesh& operator()() const;

    void writeTec360(double time, const std::string &directoryName);

private:

    std::vector<HexaFvmMesh> blocks_;

    int nProcesses_, adjacentBlockNo_[6];
};

#endif
