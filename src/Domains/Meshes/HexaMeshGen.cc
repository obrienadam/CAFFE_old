#include <fstream>
#include <boost/algorithm/string.hpp>

#include "HexaMeshGen.h"
#include "Output.h"

HexaMeshGen::HexaMeshGen()
    :
      metricConversion_(1.),
      vertices_(8)
{
    
}

HexaMeshGen::HexaMeshGen(int argc, const char *argv[])
    :
      HexaMeshGen()
{
    
    argsList_.readArgs(argc, argv);
    
}

void HexaMeshGen::processBuffer(std::string &buffer, bool removeAllWhitespace)
{
    using namespace boost::algorithm;
    
    //- Trim the leading/lagging whitespace

    trim(buffer);

    //- Remove all whitespace from buffer if desired (default)
    
    if(removeAllWhitespace)
        erase_all(buffer, " ");
    
    // Check if it is a comment line, if so discard the input
    
    if(buffer[0] == '#')
    {
        
        buffer.clear();
        
    }
    
    // Remove any comments on the line
    
    buffer = buffer.substr(0, buffer.find("#"));
    
}

void HexaMeshGen::readVertices(std::ifstream& inputFile)
{
    using namespace std;

    Point3D tempVertex;
    string buffer;

    while(!inputFile.eof())
    {
        
        //- Get a line from the file and process it

        getline(inputFile, buffer);
        processBuffer(buffer);
        
        if(buffer.empty())
        {

            continue;

        }
        else if(buffer != "{")
        {

            throw ("Expected a \"{\", but received a \"" + buffer + "\".").c_str();

        }

        //- Input is good, break

        break;

    }

    //- Begin reading the vertices

    while(true)
    {

        getline(inputFile, buffer);
        processBuffer(buffer, false);

        if(buffer.empty())
            continue;

        if(buffer == "}")
            break;


        //- Extract the bracketed coordinate

        buffer = buffer.substr(buffer.find_first_of("("), buffer.find_first_of(")") + 1);

        //- Begin extracting the vertex coordinates from the buffer

        tempVertex.x = getNextElement(buffer);
        tempVertex.y = getNextElement(buffer);
        tempVertex.z = getNextElement(buffer);
        vertices_.push_back(tempVertex);

        if(buffer != ")")
        {

            throw "Expected a \")\"";

        }

    }

    Output::printToScreen("HexaMeshGen: Successfully initialized domain vertices.");

}


void HexaMeshGen::readResolution(std::ifstream& inputFile)
{
    using namespace std;

    int nI, nJ, nK;
    string buffer;

    while(!inputFile.eof())
    {

        //- Get a line from the file and process it

        getline(inputFile, buffer);
        processBuffer(buffer);

        if(buffer.empty())
        {

            continue;

        }
        else if(buffer != "{")
        {

            throw ("Expected a \"{\", but received a \"" + buffer + "\".").c_str();

        }

        //- Input is good, break

        break;

    }

    //- Begin reading the vertices

    while(true)
    {

        getline(inputFile, buffer);
        processBuffer(buffer, false);

        if(buffer.empty())
            continue;

        if(buffer == "}")
            break;


        //- Extract the bracketed coordinate

        buffer = buffer.substr(buffer.find_first_of("("), buffer.find_first_of(")") + 1);

        //- Begin extracting the vertex coordinates from the buffer

        nI = int(getNextElement(buffer));
        nJ = int(getNextElement(buffer));
        nK = int(getNextElement(buffer));
        nodes_.allocate(nI, nJ, nK);

        if(buffer != ")")
        {

            throw "Expected a \")\"";

        }

    }

    Output::printToScreen("HexaMeshGen: Successfully allocated mesh nodes.");

}

double HexaMeshGen::getNextElement(std::string &buffer)
{
    using namespace boost::algorithm;

    double element;

    //- Ensure the buffer is trimmed

    trim_left_if(buffer, is_any_of("( "));

    //- Extract a double element from the string

    element = stod(buffer.substr(0, buffer.find_first_of(" )")));

    //- Remove the extracted element from the string

    buffer = buffer.substr(buffer.find_first_of(" )"), buffer.back());

    return element;

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
        
        if(buffer.empty())
            continue;
        
        //- Check the contents of the buffer, which must be a header. Pass the inputfile to the apropriate method
        
        if(buffer.substr(0, buffer.find("=")) == "MetricConversion")
        {
            
            //- Extract the metric conversion floating point number and store
            
            buffer = buffer.substr(buffer.find("=") + 1, buffer.back());
            
            metricConversion_ = stod(buffer);
            
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
            
            throw ("Unrecognized input field header " + buffer + ".").c_str();
            
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
