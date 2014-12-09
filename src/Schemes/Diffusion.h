#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "FvScheme.h"
#include "Field.h"
#include "Vector3D.h"

class Diffusion : public FvScheme
{

private:

    Field<double>* phi_;
    Field<Vector3D>* mu_;

public:

    Diffusion();

    void setMeshPointer(HexaFvmMesh* mesh);

};

#endif
