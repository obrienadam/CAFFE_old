#ifndef TENSOR_FIELD_H
#define TENSOR_FIELD_H

#include <string>

#include "FieldInterface.h"
#include "SmartPointer3D.h"
#include "Tensor3D.h"

class TensorField : public FieldInterface, public SmartPointer3D<Tensor3D>
{

private:

public:

    TensorField(std::string fieldName = "T(3x3)",
                int nI = 0,
                int nJ = 0,
                int nK = 0);

};

#endif
