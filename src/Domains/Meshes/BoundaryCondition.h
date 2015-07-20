#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include "Input.h"
#include "Field.h"

template<class T>
class BoundaryCondition
{
public:

    enum Location{EAST, WEST, NORTH, SOUTH, TOP, BOTTOM};

    BoundaryCondition();
    BoundaryCondition(const Input &input, Location patchLocation, Field<T> &internalField);

    virtual void setBoundaryPatch() = 0;

protected:

    Location patchLocation_;
    Field<T> &internalField_;
};

#include "BoundaryConditionI.h"

#endif
