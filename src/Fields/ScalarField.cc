#include "ScalarField.h"

ScalarField::ScalarField(std::string fieldName,
                         int nI,
                         int nJ,
                         int nK)
    :
      FieldInterface(fieldName),
      SmartPointer3D<double>(nI, nJ, nK)
{

}
