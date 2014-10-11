#include <cstdlib>
#include <algorithm>

#include "PrimitiveMesh.h"
#include "Output.h"

PrimitiveMesh::PrimitiveMesh()
{
    std::fill_n(facePatches_, 6, INTERIOR);
}

PrimitiveMesh::~PrimitiveMesh()
{

}

void PrimitiveMesh::allocate(Input &input)
{

    int nI = input.inputInts["nI"],
            nJ = input.inputInts["nJ"],
            nK = input.inputInts["nK"];

    nodes_.allocate(nI, nJ, nK);

    Output::printToScreen("Initialized nodes of PrimitiveMesh.");

    allocateFields(input);

}
