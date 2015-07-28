#include <fstream>

#include "HexaMeshGen.h"
#include "Output.h"
#include "Geometry.h"

HexaMeshGen::HexaMeshGen()
{
    using namespace boost::property_tree;

    read_info("mesh/structuredMesh.info", meshParameters_);
    meshName_ = meshParameters_.get<std::string>("MeshName");
    metricConversion_ = meshParameters_.get<double>("MetricConversion");

    for(int i = 1; i <= 8; ++i)
        vertices_[i - 1] = metricConversion_*std::stov(meshParameters_.get<std::string>("Vertices.v" + std::to_string(i)));

    resolution_ = std::stov(meshParameters_.get<std::string>("Resolution"));

    nodes_.allocate((int)resolution_.x, (int)resolution_.y, (int)resolution_.z);

    Output::print("HexaMeshGen", "Successfully read file \"mesh/structuredMesh.info\".");
}

//******************** Public methods **************************

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

    Output::print("HexaMeshGen", "Mesh generation complete.");
}

void HexaMeshGen::writeMeshFile()
{
    using namespace std;

    int i, j, k, l;

    ofstream fout("mesh/structuredMesh.dat");

    Output::print("HexaMeshGen", "Writing mesh file...");

    fout << "TITLE = " << "\"" + meshName_ + "\"" << endl
         << "VARIABLES = \"x\", \"y\", \"z\"" << endl
         << "FILETYPE = GRID" << endl
         << "ZONE I = " << nodes_.sizeI() << ", J = " << nodes_.sizeJ() << ", K = " << nodes_.sizeK() << endl
         << "DATAPACKING = BLOCK" << endl;

    for(l = 0; l < 3; ++l)
    {
        for(k = 0; k < nodes_.sizeK(); ++k)
        {
            for(j = 0; j < nodes_.sizeJ(); ++j)
            {
                for(i = 0; i < nodes_.sizeI(); ++i)
                {
                    fout << nodes_(i, j, k)(l) << " ";
                }

                fout << endl;
            }
        }
    }

    fout.close();

    Output::print("HexaMeshGen", "Finished writing mesh file.");
}

void HexaMeshGen::checkMesh()
{
    // Checks for non-planar surfaces

    if(!Geometry::checkHexahedronSurfacesIsPlanar(vertices_))
        Output::raiseException("HexaMeshGen", "checkMesh", "one or more surfaces of the mesh geometry is not planar, which is not currently supported.");
}
