#ifndef LINEAR_ADVECTION_H
#define LINEAR_ADVECTION_H

#include "FvScheme.h"
#include "Field.h"
#include "Vector3D.h"

class LinearAdvection : public FvScheme
{

private:

    Field<double>* phi_;
    Field<Vector3D>* a_;

public:

    void setMeshPointer(HexaFvmMesh* mesh);

};

#endif
