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

    ScalarFieldMap::const_iterator scalarItr;
    VectorFieldMap::const_iterator vectorItr;
    TensorFieldMap::const_iterator tensorItr;

    nodes_.allocate(nI, nJ, nK);

    Output::printToScreen("Initialized nodes of PrimitiveMesh.");

    for(scalarItr = scalarFields_.begin(); scalarItr != scalarFields_.end(); ++scalarItr)
    {

        scalarItr->second->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized scalar field \"" + scalarItr->second->fieldName + "\".");

    }

    for(vectorItr = vectorFields_.begin(); vectorItr != vectorFields_.end(); ++vectorItr)
    {

        vectorItr->second->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized vector field \"" + vectorItr->second->fieldName + "\".");

    }

    for(tensorItr = tensorFields_.begin(); tensorItr != tensorFields_.end(); ++tensorItr)
    {

        tensorItr->second->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized tensor field \"" + tensorItr->second->fieldName + "\".");

    }

    Output::printLine();

}
