#include <cstdlib>
#include <algorithm>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include "StructuredMesh.h"
#include "Output.h"


// ************* Public Methods *************

StructuredMesh::StructuredMesh()
    :
      name("UnnamedMesh"),
      nTec360Outputs_(0)
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
}

void StructuredMesh::initialize(Array3D<Point3D> &nodes)
{
    nodes_.allocate(nodes.sizeI(), nodes.sizeJ(), nodes.sizeK());

    for(int i = 0; i < nodes_.size(); ++i)
        nodes_(i) = nodes(i);
}

std::string StructuredMesh::meshStats()
{
    std::ostringstream stats;

    stats << "Mesh statistics:" << "\n"
          << "----------------" << "\n"
          << "Nodes in I direction -> " << nodes_.sizeI() << "\n"
          << "Nodes in J direction -> " << nodes_.sizeJ() << "\n"
          << "Nodes in K direction -> " << nodes_.sizeK() << "\n"
          << "Nodes total -> " << nodes_.size() << "\n";

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

    readTecplotMeshHeader(fin, name, nI, nJ, nK);

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

void StructuredMesh::resetFileStream()
{
    if(foutTec360_.is_open())
        foutTec360_.close();
}

void StructuredMesh::writeTec360(double time)
{
    int nI, nJ, nK;
    int i, j, k;

    nI = nodes_.sizeI();
    nJ = nodes_.sizeJ();
    nK = nodes_.sizeK();

    if(!foutTec360_.is_open())
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

// ************* Protected Methods *************

void StructuredMesh::readTecplotMeshHeader(std::ifstream &fin, std::string& name, int &nI, int &nJ, int &nK)
{
    using namespace std;
    using namespace boost::algorithm;

    int i;
    string buffer;
    vector<string> splitVec;

    while(!fin.eof())
    {
        getline(fin, buffer);

        splitVec.clear();
        splitVec = split(splitVec, buffer, is_any_of(" ,=\""), token_compress_on);

        if(splitVec.empty())
            continue;

        if(splitVec[0] == "TITLE")
            name = splitVec[1];
        else if(splitVec[0] == "VARIABLES")
            continue;
        else if(splitVec[0] == "ZONE")
        {
            for(i = 1; i < 7; i += 2)
            {
                if(splitVec[i] == "I")
                    nI = stoi(splitVec[i + 1]);
                else if(splitVec[i] == "J")
                    nJ = stoi(splitVec[i + 1]);
                else if(splitVec[i] == "K")
                    nK = stoi(splitVec[i + 1]);
                else
                    Output::raiseException("StructuredMesh", "initialize", "invalid zone specifier \"" + splitVec[i] + "\".");
            }
        }
        else if(splitVec[0] == "FILETYPE")
            continue;
        else if(splitVec[0] == "DATAPACKING")
            break;
        else
            Output::raiseException("StructuredMesh", "readTecplotMeshHeader", "invalid mesh header \"" + splitVec[0] + "\".");
    }
}
