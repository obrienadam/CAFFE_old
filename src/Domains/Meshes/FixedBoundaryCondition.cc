#include "FixedBoundaryCondition.h"

template<>
FixedBoundaryCondition<double>::FixedBoundaryCondition(const Input &input, const std::string &patchLocation, Field<double> &field)
    :
      BoundaryCondition<double>::BoundaryCondition(input, patchLocation, field)
{
    refValue_ = input.caseParameters.get<double>("Boundaries." + patchLocation_ + ".refValue");
}

template<>
FixedBoundaryCondition<Vector3D>::FixedBoundaryCondition(const Input &input, const std::string &patchLocation, Field<Vector3D> &field)
    :
      BoundaryCondition<Vector3D>::BoundaryCondition(input, patchLocation, field)
{
    refValue_ = std::stov(input.caseParameters.get<std::string>("Boundaries." + patchLocation_ + ".refVector"));
}
