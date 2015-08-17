#ifndef STRUCTURED_MESH_H
#define STRUCTURED_MESH_H

#include <string>
#include <fstream>
#include <vector>

#include "Input.h"
#include "Array3D.h"
#include "Point3D.h"

class StructuredMesh
{
public:

    //- Constructors and destructor
    StructuredMesh();
    ~StructuredMesh();

    //- Initialization
    virtual void initialize(const std::string &filename);
    virtual void initialize(const Array3D<Point3D> &nodes);

    /**
     * @brief Method used to initialize a simple cartesian mesh, primarily for testing
     */
    virtual void initializeCartesianMesh(double xLength, double yLength, double zLength, int nNodesI, int nNodesJ, int nNodesK);

    int nNodes() const { return nodes_.size(); }
    int nNodesI() const { return nodes_.sizeI(); }
    int nNodesJ() const { return nodes_.sizeJ(); }
    int nNodesK() const { return nodes_.sizeK(); }

    Point3D node(int i, int j, int k) const { return nodes_(i, j, k); }

    virtual std::string meshStats();

    //- Output methods
    virtual void writeTec360(double time, const std::string &directoryName);
    void resetFileStream();
    static void readTecplotMeshHeader(std::ifstream &fin, std::string &name, int& nI, int& nJ, int& nK);

    std::string name;

protected:

    //- Structured data. In this class, only the geometric mesh is stored
    Array3D<Point3D> nodes_;

    //- File output
    std::ofstream foutTec360_;
};

#endif
