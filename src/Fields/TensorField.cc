#include "TensorField.h"

TensorField::TensorField(std::string fieldName,
                         int nI,
                         int nJ,
                         int nK)
    :
      FieldInterface(fieldName),
      SmartPointer3D<Tensor3D>(nI, nJ, nK)
{

}
