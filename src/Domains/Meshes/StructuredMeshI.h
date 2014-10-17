#include <cstdlib>
#include <algorithm>
#include <sstream>

#include "StructuredMesh.h"
#include "Output.h"

template <class STATE_TYPE>
StructuredMesh<STATE_TYPE>::StructuredMesh()
{
    std::fill_n(facePatches_, 6, INTERIOR);
}

template <class STATE_TYPE>
StructuredMesh<STATE_TYPE>::~StructuredMesh()
{

}

template <class STATE_TYPE>
void StructuredMesh<STATE_TYPE>::initialize(Input &input)
{

    int nI = input.inputInts["nI"],
            nJ = input.inputInts["nJ"],
            nK = input.inputInts["nK"];

    nodes_.allocate(nI, nJ, nK);

    Output::printToScreen("Initialized nodes of StructuredMesh.");

    states_.allocate(nI, nJ, nK);

    Output::printToScreen("Initialized solution states of StructuredMesh.");

    Output::printToScreen(meshStats());

}

template <class STATE_TYPE>
int StructuredMesh<STATE_TYPE>::size()
{

    return nodes_.size();

}

template <class STATE_TYPE>
std::string StructuredMesh<STATE_TYPE>::meshStats()
{

    std::ostringstream stats;

    stats << "Mesh statistics:" << "\n"
          << "----------------" << "\n"
          << "Nodes in I direction-> " << nodes_.sizeI() << "\n"
          << "Nodes in J direction-> " << nodes_.sizeJ() << "\n"
          << "Nodes in K direction-> " << nodes_.sizeK() << "\n"
          << "Nodes total-> " << nodes_.size();

    return stats.str();

}

template <class STATE_TYPE>
typename
Field<STATE_TYPE>::iterator StructuredMesh<STATE_TYPE>::begin()
{

    return states_.begin();

}

template <class STATE_TYPE>
typename
Field<STATE_TYPE>::iterator StructuredMesh<STATE_TYPE>::end()
{

    return states_.end();

}
