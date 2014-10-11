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

    for(scalarFieldMapItr_ = scalarFields_.begin(); scalarFieldMapItr_ != scalarFields_.end(); ++scalarFieldMapItr_)
    {

        scalarFieldMapItr_->second->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized scalar field \"" + scalarFieldMapItr_->second->fieldName + "\".");

    }

    for(vectorFieldMapItr_ = vectorFields_.begin(); vectorFieldMapItr_ != vectorFields_.end(); ++vectorFieldMapItr_)
    {

        vectorFieldMapItr_->second->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized vector field \"" + vectorFieldMapItr_->second->fieldName + "\".");

    }

    for(tensorFieldMapItr_ = tensorFields_.begin(); tensorFieldMapItr_ != tensorFields_.end(); ++tensorFieldMapItr_)
    {

        tensorFieldMapItr_->second->allocate(nI, nJ, nK);

        Output::printToScreen("Initialized tensor field \"" + tensorFieldMapItr_->second->fieldName + "\".");

    }

    Output::printLine();

}
