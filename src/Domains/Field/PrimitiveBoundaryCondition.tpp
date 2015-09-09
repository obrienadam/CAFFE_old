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
    for(int i = 0; i < 6; ++i)
        types_[i] = FIXED;
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
    case PARALLEL:
        internalField_.setZeroGradientBoundaryEast();
        Parallel::iSend(Parallel::processNo(), (*adjProcNoPtr_)[0], 0, internalField_.eastBoundaryPatch);
        Parallel::iRecv((*adjProcNoPtr_)[0], Parallel::processNo(), 1, internalField_.eastBoundaryPatch);
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
    case PARALLEL:
        internalField_.setZeroGradientBoundaryWest();
        Parallel::iSend(Parallel::processNo(), (*adjProcNoPtr_)[1], 1, internalField_.westBoundaryPatch);
        Parallel::iRecv((*adjProcNoPtr_)[1], Parallel::processNo(), 0, internalField_.westBoundaryPatch);
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
    case PARALLEL:
        internalField_.setZeroGradientBoundaryNorth();
        Parallel::iSend(Parallel::processNo(), (*adjProcNoPtr_)[2], 2, internalField_.northBoundaryPatch);
        Parallel::iRecv((*adjProcNoPtr_)[2], Parallel::processNo(), 3, internalField_.northBoundaryPatch);
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
    case PARALLEL:
        internalField_.setZeroGradientBoundarySouth();
        Parallel::iSend(Parallel::processNo(), (*adjProcNoPtr_)[3], 3, internalField_.southBoundaryPatch);
        Parallel::iRecv((*adjProcNoPtr_)[3], Parallel::processNo(), 2, internalField_.southBoundaryPatch);
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
    case PARALLEL:
        internalField_.setZeroGradientBoundaryTop();
        Parallel::iSend(Parallel::processNo(), (*adjProcNoPtr_)[4], 4, internalField_.topBoundaryPatch);
        Parallel::iRecv((*adjProcNoPtr_)[4], Parallel::processNo(), 5, internalField_.topBoundaryPatch);
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
    case PARALLEL:
        internalField_.setZeroGradientBoundaryBottom();
        Parallel::iSend(Parallel::processNo(), (*adjProcNoPtr_)[5], 5, internalField_.bottomBoundaryPatch);
        Parallel::iRecv((*adjProcNoPtr_)[5], Parallel::processNo(), 4, internalField_.bottomBoundaryPatch);
        break;
    }

    Parallel::waitAll();
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
            a[1] = 0.;
            break;

        case ZERO_GRADIENT:
            a[0] += a[1];
            a[1] = 0.;
            break;

        case EMPTY:
            a[1] = 0.;
            break;

        case PARALLEL:
            break;
        }
    }
    if(i == 0)
    {
        switch(types_[1])
        {
        case FIXED:
            b -= a[2]*internalField_.westBoundaryPatch(0, j, k);
            a[2] = 0.;
            break;

        case ZERO_GRADIENT:
            a[0] += a[2];
            a[2] = 0.;
            break;

        case EMPTY:
            a[2] = 0.;
            break;

        case PARALLEL:
            break;
        }
    }

    //- North/south boundaries
    if(j == mesh_.uCellJ())
    {
        switch(types_[2])
        {
        case FIXED:
            b -= a[3]*internalField_.northBoundaryPatch(i, 0, k);
            a[3] = 0.;
            break;

        case ZERO_GRADIENT:
            a[0] += a[3];
            a[3] = 0.;
            break;

        case EMPTY:
            a[3] = 0.;
            break;

        case PARALLEL:
            break;
        }
    }
    if(j == 0)
    {
        switch(types_[3])
        {
        case FIXED:
            b -= a[4]*internalField_.southBoundaryPatch(i, 0, k);
            a[4] = 0.;
            break;

        case ZERO_GRADIENT:
            a[0] += a[4];
            a[4] = 0.;
            break;

        case EMPTY:
            a[4] = 0.;
            break;

        case PARALLEL:
            break;
        }
    }

    //- Top/bottom boundaries
    if(k == mesh_.uCellK())
    {
        switch(types_[4])
        {
        case FIXED:
            b -= a[5]*internalField_.topBoundaryPatch(i, j, 0);
            a[5] = 0.;
            break;

        case ZERO_GRADIENT:
            a[0] += a[5];
            a[5] = 0.;
            break;

        case EMPTY:
            a[5] = 0.;
            break;

        case PARALLEL:
            break;
        }
    }
    if(k == 0)
    {
        switch(types_[5])
        {
        case FIXED:
            b -= a[6]*internalField_.bottomBoundaryPatch(i, j, 0);
            a[6] = 0.;
            break;

        case ZERO_GRADIENT:
            a[0] += a[6];
            a[6] = 0.;
            break;

        case EMPTY:
            a[6] = 0.;
            break;

        case PARALLEL:
            break;
        }
    }
}

template<class T>
void PrimitiveBoundaryCondition<T>::changeType(int i, Type newType, const T &refValue)
{
    types_[i] = newType;
    refValues_[i] = refValue;

    if(newType == FIXED)
        setFixedBoundaries(); // This is a bit wasteful
}

template <class T>
void PrimitiveBoundaryCondition<T>::setParallelBoundaries(std::shared_ptr<std::array<int, 6> > adjProcNoPtr)
{
    if(!adjProcNoPtr)
        return;

    adjProcNoPtr_ = adjProcNoPtr;

    for(int i = 0; i < 6; ++i)
    {
        if((*adjProcNoPtr_)[i] != Parallel::PROC_NULL)
            changeType(i, PARALLEL, T());
    }
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
