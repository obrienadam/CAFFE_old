#ifndef PRIMITIVE_BOUNDARY_CONDITION_H
#define PRIMITIVE_BOUNDARY_CONDITION_H

#include <memory>
#include <array>

#include "Parallel.h"
#include "Field.h"

template<class T>
class PrimitiveBoundaryCondition
{

public:

    enum Type{FIXED, ZERO_GRADIENT, EMPTY, PARALLEL};

    PrimitiveBoundaryCondition(const Input &input, Field<T> &field);
    PrimitiveBoundaryCondition(Field<T> &field);

    virtual void setBoundaries();
    virtual void setImplicitBoundaryCoefficients(int i, int j, int k, double *a, T &b);

    Type getType(int i) const { return types_[i]; }
    Type getTypeEast() const { return types_[0]; }
    Type getTypeWest() const { return types_[1]; }
    Type getTypeNorth() const { return types_[2]; }
    Type getTypeSouth() const { return types_[3]; }
    Type getTypeTop() const { return types_[4]; }
    Type getTypeBottom() const { return types_[5]; }

    void changeType(int i, Type newType, const T &refValue);
    void setParallelBoundaries(std::shared_ptr< std::array<int, 6> > adjProcNoPtr);

protected:

    virtual void setFixedBoundaries();

    Field<T> &internalField_;

    Type types_[6];
    T refValues_[6];

    const HexaFvmMesh &mesh_;
    std::shared_ptr< std::array<int, 6> > adjProcNoPtr_;

};

#include "PrimitiveBoundaryCondition.tpp"

#endif
