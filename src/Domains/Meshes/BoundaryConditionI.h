#include "BoundaryCondition.h"

template<class T>
BoundaryCondition<T>::BoundaryCondition(const Input &input, Location patchLocation, Field<T> &internalField)
    :
      patchLocation_(patchLocation),
      internalField_(internalField)
{

}
