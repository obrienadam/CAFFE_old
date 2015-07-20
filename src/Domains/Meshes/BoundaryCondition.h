#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <string>

#include "Input.h"
#include "Field.h"

template<class T>
class BoundaryCondition
{
public:

    BoundaryCondition();
    BoundaryCondition(const Input &input, const std::string &patchLocation, Field<T> &internalField);

    virtual void setBoundaryPatch() = 0;

protected:

    std::string patchLocation_;
    Field<T> &internalField_;
};

#include "BoundaryConditionI.h"

#endif
