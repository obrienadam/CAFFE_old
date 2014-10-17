#ifndef STRUCTURED_MESH_H
#define STRUCTURED_MESH_H

#include <string>

#include "DomainInterface.h"
#include "Point3D.h"
#include "Array3D.h"
#include "Field.h"

enum Patch{BOUNDARY, INTERIOR};
enum Face{EAST = 0, WEST = 1, NORTH = 2, SOUTH = 3, TOP = 4, BOTTOM = 5};

template <class STATE_TYPE>
class StructuredMesh : public DomainInterface<STATE_TYPE>
{

protected:

    //- Structured data

    Array3D<Point3D> nodes_;
    Field<STATE_TYPE> states_;

    //- Patches used for boundary conditions and communication

    Patch facePatches_[6];

public:

    //- Constructors and destructor

    StructuredMesh();
    ~StructuredMesh();

    //- Initialization

    virtual void initialize(Input& input);
    virtual int size();
    virtual std::string meshStats();

    typedef typename Field<STATE_TYPE>::iterator iterator;

    iterator begin();
    iterator end();

};

#include "StructuredMeshI.h"

#endif
