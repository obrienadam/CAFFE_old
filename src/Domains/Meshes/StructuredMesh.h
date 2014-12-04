#ifndef STRUCTURED_MESH_H
#define STRUCTURED_MESH_H

#include <string>
#include <fstream>

#include "DomainInterface.h"
#include "Point3D.h"
#include "Array3D.h"

enum Patch{BOUNDARY, INTERIOR};
enum Face{EAST = 0, WEST = 1, NORTH = 2, SOUTH = 3, TOP = 4, BOTTOM = 5};

template <class STATE_TYPE>
class StructuredMesh : public DomainInterface<STATE_TYPE>
{

protected:

    //- Structured data

    Array3D<Point3D> nodes_;
    Array3D<STATE_TYPE> states_;

    //- Patches used for boundary conditions and communication

    Patch facePatches_[6];

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
    virtual int size();
    virtual std::string meshStats();

    //- Initialize from input file

    virtual void initialize(std::string filename);

    //- Compute the time derivatives

    virtual void computeTimeDerivatives(STATE_TYPE* timeDerivatives);

    //- Iterators methods

    typedef typename Array3D<STATE_TYPE>::iterator iterator;

    iterator begin();
    iterator end();

    //- Output methods

    virtual void write(double time = 0.);
    virtual void writeTec360(double time = 0.);

};

#include "StructuredMeshI.h"

#endif
