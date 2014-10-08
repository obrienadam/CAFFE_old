#ifndef VECTOR_FIELD_H
#define VECTOR_FIELD_H

#include <string>

#include "FieldInterface.h"
#include "SmartPointer3D.h"
#include "Vector3D.h"

class VectorField : public FieldInterface, public SmartPointer3D<Vector3D>
{

private:

public:

    VectorField(std::string fieldName = "V(xi, yj, zk)",
                int nI = 0,
                int nJ = 0,
                int nK = 0);

};

#endif
