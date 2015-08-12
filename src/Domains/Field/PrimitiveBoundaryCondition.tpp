#include "PrimitiveBoundaryCondition.h"

template<class T>
PrimitiveBoundaryCondition<T>::PrimitiveBoundaryCondition(const Input &input, Field<T> &field)
    :
      internalField_(field),
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
            Output::raiseException("FlowBoundaryCondition", "FlowBoundaryCondition", "invalid boundary type \"" + typeStr + "\".");

        refValues_[i] = input.caseParameters.get<T>("Boundaries." + locationStr[i] + ".refValue");
    }

    setFixedBoundaries();
}

template<class T>
PrimitiveBoundaryCondition<T>::PrimitiveBoundaryCondition(Field<T> &field)
    :
      internalField_(field),
      mesh_(internalField_.getMesh())
{

}

template<class T>
void PrimitiveBoundaryCondition<T>::setBoundaries()
{
    //- East boundary
    switch(types_[0])
    {
    case FIXED:
        break;

    case ZERO_GRADIENT: case EMPTY:
        internalField_.setZeroGradientBoundaryEast();
        break;
    }

    //- West boundary
    switch(types_[1])
    {
    case FIXED:
        break;

    case ZERO_GRADIENT: case EMPTY:
        internalField_.setZeroGradientBoundaryWest();
        break;
    }

    //- North boundary
    switch(types_[2])
    {
    case FIXED:
        break;

    case ZERO_GRADIENT: case EMPTY:
        internalField_.setZeroGradientBoundaryNorth();
        break;
    }

    //- South boundary
    switch(types_[3])
    {
    case FIXED:
        break;

    case ZERO_GRADIENT: case EMPTY:
        internalField_.setZeroGradientBoundarySouth();
        break;
    }

    //- Top boundary
    switch(types_[4])
    {
    case FIXED:
        break;

    case ZERO_GRADIENT: case EMPTY:
        internalField_.setZeroGradientBoundaryTop();
        break;
    }

    //- Bottom boundary
    switch(types_[5])
    {
    case FIXED:
        break;

    case ZERO_GRADIENT: case EMPTY:
        internalField_.setZeroGradientBoundaryBottom();
        break;
    }
}

template<class T>
void PrimitiveBoundaryCondition<T>::setImplicitBoundaryCoefficients(int i, int j, int k, double a[], T &b)
{
    //- East/west boundaries
    if(i == mesh_.uCellI())
    {
        switch(types_[0])
        {
        case FIXED:
            b -= a[1]*internalField_.eastBoundaryPatch(0, j, k);
            break;

        case ZERO_GRADIENT:
            a[0] += a[1];
            break;

        case EMPTY:
            break;
        }
        a[1] = 0.;
    }
    if(i == 0)
    {
        switch(types_[1])
        {
        case FIXED:
            b -= a[2]*internalField_.westBoundaryPatch(0, j, k);
            break;

        case ZERO_GRADIENT:
            a[0] += a[2];
            break;

        case EMPTY:
            break;
        }
        a[2] = 0.;
    }

    //- North/south boundaries
    if(j == mesh_.uCellJ())
    {
        switch(types_[2])
        {
        case FIXED:
            b -= a[3]*internalField_.northBoundaryPatch(i, 0, k);
            break;

        case ZERO_GRADIENT:
            a[0] += a[3];
            break;

        case EMPTY:
            break;
        }
        a[3] = 0.;
    }
    if(j == 0)
    {
        switch(types_[3])
        {
        case FIXED:
            b -= a[4]*internalField_.southBoundaryPatch(i, 0, k);
            break;

        case ZERO_GRADIENT:
            a[0] += a[4];
            break;

        case EMPTY:
            break;
        }
        a[4] = 0.;
    }

    //- Top/bottom boundaries
    if(k == mesh_.uCellK())
    {
        switch(types_[4])
        {
        case FIXED:
            b -= a[5]*internalField_.topBoundaryPatch(i, j, 0);
            break;

        case ZERO_GRADIENT:
            a[0] += a[5];
            break;

        case EMPTY:
            break;
        }
        a[5] = 0.;
    }
    if(k == 0)
    {
        switch(types_[5])
        {
        case FIXED:
            b -= a[6]*internalField_.bottomBoundaryPatch(i, j, 0);
            break;

        case ZERO_GRADIENT:
            a[0] += a[6];
            break;

        case EMPTY:
            break;
        }
        a[6] = 0.;
    }
}

template<class T>
void PrimitiveBoundaryCondition<T>::changeType(int i, Type newType, T refValue)
{
    types_[i] = newType;
    refValues_[i] = refValue;

    if(newType == FIXED)
        setFixedBoundaries(); // This is a bit wasteful
}

//************************ Protected methods ***********************************

template<class T>
void PrimitiveBoundaryCondition<T>::setFixedBoundaries()
{
    internalField_.setFixedBoundaryPatches(refValues_);
    internalField_.setEastFacesFromPatch();
    internalField_.setWestFacesFromPatch();
    internalField_.setNorthFacesFromPatch();
    internalField_.setSouthFacesFromPatch();
    internalField_.setTopFacesFromPatch();
    internalField_.setBottomFacesFromPatch();
}
