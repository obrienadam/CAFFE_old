#include <cstdlib>
#include <algorithm>

#include "StructuredMesh.h"
#include "Output.h"

StructuredMesh::StructuredMesh()
{
    std::fill_n(facePatches_, 6, INTERIOR);
}

StructuredMesh::~StructuredMesh()
{

}

void StructuredMesh::allocate(Input &input)
{

    int nI = input.inputInts["nI"],
            nJ = input.inputInts["nJ"],
            nK = input.inputInts["nK"];

    nodes_.allocate(nI, nJ, nK);

    Output::printToScreen("Initialized nodes of StructuredMesh.");

}

void StructuredMesh::size()
{

    return nodes_.size();

}
