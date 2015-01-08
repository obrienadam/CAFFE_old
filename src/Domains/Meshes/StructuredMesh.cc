#include <cstdlib>
#include <algorithm>
#include <sstream>

#include "StructuredMesh.h"
#include "Output.h"
#include "InputStringProcessing.h"

// ************* Private Methods *************

// ************* Public Methods *************

StructuredMesh::StructuredMesh()
    :
      name("UnnamedMesh")
{

}

StructuredMesh::~StructuredMesh()
{

    if(foutRestart_.is_open())
        foutRestart_.close();

    if(foutTec360_.is_open())
        foutTec360_.close();

}

void StructuredMesh::initialize(Input &input)
{

    Output::printToScreen("Initializing structured mesh...");

    initialize(input.inputStrings["domainFile"]);

    Output::printToScreen("Initialization of structured mesh complete.");

    Output::printToScreen(meshStats());

}

int StructuredMesh::size()
{

    return nodes_.size();

}

std::string StructuredMesh::meshStats()
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

void StructuredMesh::initialize(std::string filename)
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
            break;

        }
        else
        {

            Output::raiseException("StructuredMesh", "initialize", "unrecognized structured mesh dimension \"" + bufferVec[0] + "\".");

        }

    }

    nodes_.allocate(nI, nJ, nK);

    // Check to see if a dimension was not allocated

    if(nI == 0 || nJ == 0 || nK == 0)
    {

        Output::raiseException("StructuredMesh", "initialize", "one or more structured mesh dimensions was not found in file \"" + filename + "\".");

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
            ++k;

        }

    }

    fin.close();

    write();
    writeTec360();

}

void StructuredMesh::write(double time)
{

    int nI, nJ, nK;
    int i, j, k;

    nI = nodes_.sizeI();
    nJ = nodes_.sizeJ();
    nK = nodes_.sizeK();

    if(DomainInterface::nOutputs_ == 0)
    {

        foutRestart_.open((name + ".msh").c_str());

        foutRestart_ << "nNodesI=" << nI << "\n"
                     << "nNodesJ=" << nJ << "\n"
                     << "nNodesK=" << nK << "\n";

    }

    for(k = 0; k < nK; ++k)
    {

        for(j = 0; j < nJ; ++j)

        {

            for(i = 0; i < nI; ++i)

            {

                 foutRestart_ << "(" << nodes_(i, j, k) << ") ";

            } // end for i

            foutRestart_ << "\n";

        } // end for j
    } // end for k

    foutRestart_.close();

}

void StructuredMesh::writeTec360(double time)
{

    int nI, nJ, nK;
    int i, j, k;

    nI = nodes_.sizeI();
    nJ = nodes_.sizeJ();
    nK = nodes_.sizeK();

    if(nTec360Outputs_ == 0)
    {

        foutTec360_.open((name + ".dat").c_str());

        foutTec360_ << "TITLE=" << name << "\n"
                    << "STRANDID=1, SOLUTIONTIME=" << time << "\n"
                    << "VARIABLES = " << "\"x\", \"y\", \"z\", " << "\n"
                    << "ZONE I=" << nI << ", J=" << nJ << ", K=" << nK << "\n"
                    << "ZONETYPE=ORDERED, DATAPACKING=POINT" << "\n\n";

    }

    for(k = 0; k < nK; ++k)
    {

        for(j = 0; j < nJ; ++j)

        {

            for(i = 0; i < nI; ++i)

            {

                 foutTec360_ << nodes_(i, j, k).x << " " << nodes_(i, j, k).y << " " << nodes_(i, j, k).z << "\n";

            } // end for i
        } // end for j
    } // end for k

    ++nTec360Outputs_;

    foutTec360_.close();

}
