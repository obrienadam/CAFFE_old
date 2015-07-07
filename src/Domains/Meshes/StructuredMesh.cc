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
    Output::print("StructuredMesh", "Initializing structured mesh...");

    initialize("mesh/structuredMesh.dat");

    Output::print("StructuredMesh", "Initialization of structured mesh complete.");
    Output::print(meshStats());
}

void StructuredMesh::initialize(Array3D<Point3D> &nodes)
{
    nodes_.allocate(nodes.sizeI(), nodes.sizeJ(), nodes.sizeK());

    for(int i = 0; i < nodes_.size(); ++i)
        nodes_(i) = nodes(i);
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
          << "Nodes total-> " << nodes_.size() << "\n";

    return stats.str();
}

void StructuredMesh::initialize(std::string filename)
{
    using namespace std;

    string buffer;
    vector<string> bufferVec;
    ifstream fin;

    int i, j, k, l;
    int nI, nJ, nK;

    fin.open(filename.c_str());

    nI = nJ = nK = 0;

    while(!fin.eof())
    {
        getline(fin, buffer);
        buffer = InputStringProcessing::processBuffer(buffer, false);

        if(buffer.empty())
            continue;

        bufferVec = InputStringProcessing::partition(buffer, " ,=\"");

        if(bufferVec[0] == "TITLE")
            name = bufferVec[1];
        else if(bufferVec[0] == "VARIABLES")
            continue;
        else if(bufferVec[0] == "ZONE")
        {
            for(i = 1; i < 7; i += 2)
            {
                if(bufferVec[i] == "I")
                    nI = stoi(bufferVec[i + 1]);
                else if(bufferVec[i] == "J")
                    nJ = stoi(bufferVec[i + 1]);
                else if(bufferVec[i] == "K")
                    nK = stoi(bufferVec[i + 1]);
                else
                    Output::raiseException("StructuredMesh", "initialize", "invalid zone specifier \"" + bufferVec[i] + "\".");
            }
        }
        else if(bufferVec[0] == "FILETYPE")
            continue;
        else if(bufferVec[0] == "DATAPACKING")
            break;
        else
            Output::raiseException("StructuredMesh", "initialize", "invalid mesh header \"" + bufferVec[0] + "\".");
    }

    // Check to see if a dimension was not allocated
    if(nI == 0 || nJ == 0 || nK == 0)
        Output::raiseException("StructuredMesh", "initialize", "one or more structured mesh dimensions was not found in file \"" + filename + "\".");

    nodes_.allocate(nI, nJ, nK);

    for(l = 0; l < 3; ++l)
    {
        for(k = 0; k < nK; ++k)
        {
            for(j = 0; j < nJ; ++j)
            {
                for(i = 0; i < nI; ++i)
                {
                    fin >> nodes_(i, j, k)(l);
                }
            }
        }
    }

    fin.close();
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
