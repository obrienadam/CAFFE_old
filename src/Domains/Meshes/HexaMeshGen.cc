#include <fstream>
#include <boost/algorithm/string/erase.hpp>

#include "HexaMeshGen.h"

HexaMeshGen::HexaMeshGen()
    :
      metricConversion_(1.),
      vertices_(2, 2, 2)
{

}

HexaMeshGen::HexaMeshGen(int argc, const char *argv[])
    :
      HexaMeshGen()
{

    argsList_.readArgs(argc, argv);

}

void HexaMeshGen::processBuffer(std::string &buffer)
{
    using namespace boost::algorithm;

    //- Remove whitespace from the buffer

    erase_all(buffer, " ");

    // Check if it is a comment line, if so discard the input

    if(buffer[0] == '#')
    {

        buffer.clear();

    }

    // Remove any comments on the line

    buffer = buffer.substr(0, buffer.find("#"));

}

void HexaMeshGen::readVertices(std::ifstream& inFile)
{

}

void HexaMeshGen::readResolution(std::ifstream& inFile)
{

}

void HexaMeshGen::readMeshInputFile()
{
    using namespace std;

    string buffer;
    ifstream inputFile(argsList_.varsMap_["file"].as<string>().c_str());

    if(!inputFile.is_open())
    {

        throw "Mesh input file not found.";

    }

    while(!inputFile.eof())
    {

        //- Get a line from the buffer and process it

        getline(inputFile, buffer);
        processBuffer(buffer);

        //- Check to see if the buffer is empty

        if(buffer == "")
            continue;

        //- Check the contents of the buffer, which must be a header. Pass the inputfile to the apropriate method

        if(buffer.substr(0, buffer.find("=")) == "MetricConversion")
        {



        }
        else if(buffer == "Vertices")
        {

            readVertices(inputFile);

        }
        else if(buffer == "Resolution")
        {

            readResolution(inputFile);

        }
        else
        {

            throw "A problem occurred while trying to read a field header input.";

        }

    }

}

void HexaMeshGen::writeMeshFile()
{

}

void HexaMeshGen::generateMesh()
{

}

void HexaMeshGen::generateBoxMesh(double dx, double dy, double dz)
{



}
