#include <fstream>

#include "HexaMeshGen.h"
#include "Output.h"
#include "InputStringProcessing.h"
#include "Geometry.h"

HexaMeshGen::HexaMeshGen()
    :
      metricConversion_(1.)
{
    
}

void HexaMeshGen::readVertices(std::ifstream& inputFile)
{
    using namespace std;

    string buffer;

    // Make sure the vertex list is empty since push_back is used

    vertices_.clear();

    while(!inputFile.eof())
    {
        // Get a line from the file and process it

        getline(inputFile, buffer);
        buffer = InputStringProcessing::processBuffer(buffer);
        
        if(buffer.empty())
        {
            continue;
        }
        else if(buffer != "{")
        {
            Output::raiseException("HexaMeshGen", "readVertices", "Expected a \"{\", but received a \"" + buffer + "\".");
        }

        // Input is good, break

        break;
    }

    // Begin reading the vertices

    while(true)
    {
        getline(inputFile, buffer);
        buffer = InputStringProcessing::processBuffer(buffer, false);

        if(buffer.empty())
            continue;

        if(buffer == "}")
            break;

        // buffer should contain a bracketed vector

        vertices_.push_back(Point3D(buffer)*metricConversion_);
    }

    Output::print("HexaMeshGen: Successfully initialized domain vertices.");

    checkMesh();
}

void HexaMeshGen::readResolution(std::ifstream& inputFile)
{
    using namespace std;

    int nI, nJ, nK;
    string buffer;

    while(!inputFile.eof())
    {
        // Get a line from the file and process it

        getline(inputFile, buffer);
        buffer = InputStringProcessing::processBuffer(buffer);

        if(buffer.empty())
        {
            continue;
        }
        else if(buffer != "{")
        {
            Output::raiseException("HexaMeshGen", "readResolution", "Expected a \"{\", but received a \"" + buffer + "\".");
        }

        // Input is good, break

        break;
    }

    // Begin reading the vertices

    while(true)
    {
        getline(inputFile, buffer);
        buffer = InputStringProcessing::processBuffer(buffer, false);

        if(buffer.empty())
            continue;

        if(buffer == "}")
            break;


        // Extract the bracketed coordinate

        buffer = buffer.substr(buffer.find_first_of("("), buffer.find_first_of(")") + 1);

        // Begin extracting the vertex coordinates from the buffer

        nI = int(InputStringProcessing::getNextElement(buffer));
        nJ = int(InputStringProcessing::getNextElement(buffer));
        nK = int(InputStringProcessing::getNextElement(buffer));
        nodes_.allocate(nI, nJ, nK);

        if(buffer != ")")
        {
            Output::raiseException("HexaMeshGen", "readResolution", "expected a \")\"");
        }
    }

    Output::print("HexaMeshGen: Successfully allocated mesh nodes.");
}

void HexaMeshGen::readMeshInputFile(std::string filename)
{
    using namespace std;
    
    string buffer;
    ifstream inputFile(filename.c_str());
    
    if(!inputFile.is_open())
    {      
        Output::raiseException("HexaMeshGen", "readMeshInputFile", "mesh input file not found.");
    }
    
    while(!inputFile.eof())
    {
        // Get a line from the buffer and process it
        
        getline(inputFile, buffer);
        buffer = InputStringProcessing::processBuffer(buffer);
        
        // Check to see if the buffer is empty
        
        if(buffer.empty())
            continue;
        
        // Check the contents of the buffer, which must be a header. Pass the inputfile to the apropriate method
        
        if(buffer.substr(0, buffer.find("=")) == "MetricConversion")
        {
            // Extract the metric conversion floating point number and store
            
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
            Output::raiseException("HexaMeshGen", "readMeshInputFile", "Unrecognized input field header " + buffer + ".");
        }
    }
}

void HexaMeshGen::writeMeshFile()
{
    using namespace std;

    int i, j, k;

    ofstream fout("mesh/structuredMesh.msh");

    fout << "# CAFFE structured mesh file\n\n"
         << "nNodesI=" << nodes_.sizeI() << endl
         << "nNodesJ=" << nodes_.sizeJ() << endl
         << "nNodesK=" << nodes_.sizeK() << endl;

    for(k = 0; k < nodes_.sizeK(); ++k)
    {
        for(j = 0; j < nodes_.sizeJ(); ++j)
        {
            for(i = 0; i < nodes_.sizeI(); ++i)
            {
                fout << "(" << nodes_(i, j, k) << ") ";
            } // end for i

            fout << endl;
        } // end for j
    } // end for k

    fout.close();
}

void HexaMeshGen::generateMesh()
{
    using namespace std;

    int i, j, k, upperI(nodes_.sizeI() - 1), upperJ(nodes_.sizeJ() - 1), upperK(nodes_.sizeK() - 1);
    double s;
    Vector3D tmp1, tmp2;

    // Generate surface meshes on the west and east sides

    for(k = 0; k <= upperK; ++k)
    {
        s = double(k)/double(upperK);

        nodes_(0, 0, k) = (vertices_[4] - vertices_[0])*s + vertices_[0];
        nodes_(0, upperJ, k) = (vertices_[7] - vertices_[3])*s + vertices_[3];
        nodes_(upperI, 0, k) = (vertices_[5] - vertices_[1])*s + vertices_[1];
        nodes_(upperI, upperJ, k) = (vertices_[6] - vertices_[2])*s + vertices_[2];

        tmp1 = nodes_(0, upperJ, k) - nodes_(0, 0, k);
        tmp2 = nodes_(upperI, upperJ, k) - nodes_(upperI, 0, k);

        for(j = 0; j <= upperJ; ++j)
        {
            s = double(j)/double(upperJ);

            nodes_(0, j, k) = tmp1*s + nodes_(0, 0, k);
            nodes_(upperI, j, k) = tmp2*s + nodes_(upperI, 0, k);
        } // end for j
    } // end for k

    // Generate the 3D mesh using the opposing surface meshes

    for(k = 0; k <= upperK; ++k)
    {
        for(j = 0; j <= upperJ; ++j)
        {
            tmp1 = nodes_(upperI, j, k) - nodes_(0, j, k);

            for(i = 0; i <= upperI; ++i)
            {
                s = double(i)/double(upperI);

                nodes_(i, j, k) = tmp1*s + nodes_(0, j, k);
            } // end for i
        } // end for j
    } // end for k

    Output::print("HexaMeshGen: Mesh generation complete.");
}

void HexaMeshGen::generateBoxMesh(double dx, double dy, double dz)
{
    int i, j, k;

    for(k = 0; k < nodes_.sizeK(); ++k)
    {
        for(j = 0; j < nodes_.sizeJ(); ++j)
        {
            for(i = 0; i < nodes_.sizeI(); ++i)
            {

                nodes_(i, j, k) = Point3D( dx*double(i)/double(nodes_.sizeI() - 1),
                                           dy*double(j)/double(nodes_.sizeJ() - 1),
                                           dz*double(k)/double(nodes_.sizeK() - 1) );
            } // end for i
        } // end for j
    } // end for k
}

void HexaMeshGen::checkMesh()
{
    // Checks for non-planar surfaces

    if(!Geometry::checkHexahedronSurfacesIsPlanar(vertices_.data()))
        Output::raiseException("HexaMeshGen", "checkMesh", "one or more surfaces of the mesh geometry is not planar, which is not currently supported.");
}
