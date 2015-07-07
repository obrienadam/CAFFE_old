#ifndef STRUCTURED_MESH_H
#define STRUCTURED_MESH_H

#include <string>
#include <fstream>
#include <vector>

#include "DomainInterface.h"
#include "Point3D.h"
#include "Array3D.h"

class StructuredMesh : public DomainInterface
{
protected:

    //- Structured data. In this class, only the geometric mesh is stored

    Array3D<Point3D> nodes_;

    //- File output object

    std::ofstream foutRestart_, foutTec360_;

public:

    //- Constructors and destructor

    StructuredMesh();
    ~StructuredMesh();

    //- Mesh name

    std::string name;

    //- Initialization

    virtual void initialize(Input& input);
    virtual void initialize(Array3D<Point3D>& nodes);
    virtual int size();
    virtual std::string meshStats();

    //- Initialize from input file

    virtual void initialize(std::string filename);

    //- Output methods
    virtual void writeTec360(double time = 0.);
};

#endif
