#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <vector>

#include "StructuredMesh.h"
#include "Output.h"
#include "InputStringProcessing.h"

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
void StructuredMesh<STATE_TYPE>::initialize(std::string filename)
{
    using namespace std;

    string buffer;
    vector<string> bufferVec;
    ifstream fin;

    int i, j, k;
    int nI(0), nJ(0), nK(0);

    fin.open(filename.c_str());

    while(!fin.eof())
    {

        getline(fin, buffer);
        buffer = InputStringProcessing::processBuffer(buffer);

        if(buffer.empty())
            continue;

        bufferVec = InputStringProcessing::partition(buffer, "=");

        if(bufferVec[0] == "nNodesI")
        {

            nI = stoi(bufferVec[1]);

        }
        else if(bufferVec[0] == "nNodesJ")
        {

            nJ = stoi(bufferVec[1]);

        }
        else if(bufferVec[0] == "nNodesK")
        {

            nK = stoi(bufferVec[1]);

        }
        else
        {

            throw ("Unrecognized structured mesh dimension \"" + bufferVec[0] + "\".").c_str();

        }

    }

    // Check to see if a dimension was not allocated

    if(nI == 0 || nJ == 0 || nK == 0)
    {

        throw ("One or more structured mesh dimensions were not found in file \"" + filename + "\".").c_str();

    }

    i = j = k = 0;

    while(!fin.eof())
    {

        getline(fin, buffer);
        buffer = InputStringProcessing::processBuffer(buffer, false);

        if(buffer.empty())
            continue;

        bufferVec = InputStringProcessing::partition(buffer, " ");

        for(i = 0; i < nI; ++i)
        {

            nodes_(i, j, k) = Point3D(bufferVec[i]);

        }

        ++j;

        if(j == nJ)
        {

            j = 0;
            ++nK;

        }

    }

    fin.close();

}

template <class STATE_TYPE>
void StructuredMesh<STATE_TYPE>::computeTimeDerivatives(STATE_TYPE* timeDerivatives)
{



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

template <class STATE_TYPE>
void StructuredMesh<STATE_TYPE>::outputData(double time)
{



}
