#include "BoundaryCondition.h"

template<class T>
BoundaryCondition<T>::BoundaryCondition(const Input &input, const std::string &patchLocation, Field<T> &internalField)
    :
      patchLocation_(patchLocation),
      internalField_(internalField)
{

}
