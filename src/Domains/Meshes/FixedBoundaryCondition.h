#ifndef FIXED_BOUNDARY_CONDITION_H
#define FIXED_BOUNDARY_CONDITION_H

#include "BoundaryCondition.h"

template<class T>
class FixedBoundaryCondition : public BoundaryCondition<T>
{
public:

    FixedBoundaryCondition(const Input &input, const std::string &patchLocation, Field<T>& field);

private:

    T refValue_;

};

#endif
