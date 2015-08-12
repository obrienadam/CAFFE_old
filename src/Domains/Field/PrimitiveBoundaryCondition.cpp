#include "PrimitiveBoundaryCondition.h"
#include "Vector3D.h"

template<>
PrimitiveBoundaryCondition<Vector3D>::PrimitiveBoundaryCondition(const Input &input, Field<Vector3D> &internalField_)
    :
      internalField_(internalField_),
      mesh_(internalField_.getMesh())
{
    using namespace std;

    int i;
    string locationStr[6] = {"east", "west", "north", "south", "top", "bottom"}, typeStr;

    for(i = 0; i < 6; ++i)
    {
        typeStr = input.caseParameters.get<string>("Boundaries." + locationStr[i] + ".type");

        if(typeStr == "fixed")
            types_[i] = FIXED;
        else if(typeStr == "zeroGradient")
            types_[i] = ZERO_GRADIENT;
        else if(typeStr == "empty")
            types_[i] = EMPTY;
        else
            Output::raiseException("PrimitiveBoundaryCondition", "PrimitiveBoundaryCondition", "invalid boundary type \"" + typeStr + "\".");

        refValues_[i] = std::stov(input.caseParameters.get<string>("Boundaries." + locationStr[i] + ".refVector"));
    }

    setFixedBoundaries();
}
