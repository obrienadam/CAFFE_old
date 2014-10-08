#include "VectorField.h"

VectorField::VectorField(std::string fieldName,
                         int nI,
                         int nJ,
                         int nK)
    :
      FieldInterface(fieldName),
      SmartPointer3D<Vector3D>(nI, nJ, nK)
{

}
