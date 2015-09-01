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
    virtual void initialize(const std::vector<Point3D> &vertices, int nNodesI, int nNodesJ, int nNodesK);

    /**
     * @brief Method used to initialize a simple cartesian mesh, primarily for testing
     */
    virtual void initializeCartesianMesh(double xLength, double yLength, double zLength, int nNodesI, int nNodesJ, int nNodesK);

    int nNodes() const { return nodes_.size(); }
    int nNodesI() const { return nodes_.sizeI(); }
    int nNodesJ() const { return nodes_.sizeJ(); }
    int nNodesK() const { return nodes_.sizeK(); }

    int uNodeI() const { return nNodesI() - 1; }
    int uNodeJ() const { return nNodesJ() - 1; }
    int uNodeK() const { return nNodesK() - 1; }

    Point3D node(int i, int j, int k) const { return nodes_(i, j, k); }

    virtual std::string meshStats();

    //- Output methods
    virtual void writeTec360(double time, const std::string &directory);
    void resetFileStream();
    static void readTecplotMeshHeader(std::ifstream &fin, std::string &name, int& nI, int& nJ, int& nK);

    virtual void changeName(const std::string &newName);
    const std::string& name() const { return name_; }

protected:

    //- Structured data. In this class, only the geometric mesh is stored
    Array3D<Point3D> nodes_;

    //- File output
    const int OUTPUT_PRECISION;
    std::ofstream foutTec360_;
    std::string name_;
};

#endif
