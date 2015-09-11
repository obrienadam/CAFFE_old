#include "SimpleBoundaryCondition.h"

SimpleBoundaryCondition::SimpleBoundaryCondition(const Input &input,
                                                 Field<Vector3D> &uField,
                                                 Field<double> &pField,
                                                 Field<double> &rhoField,
                                                 Field<double> &muField,
                                                 Field<Vector3D> &hField,
                                                 Field<double> &dField,
                                                 Field<double> &pCorrField)
    :
      FlowBoundaryConditions(input, uField, pField, rhoField, muField, hField, dField),
      pCorrFieldBcs(pCorrField)
{
    using namespace std;

    pCorrFieldBcs.changeType(0, pFieldBcs.getTypeEast(), 0.);
    pCorrFieldBcs.changeType(1, pFieldBcs.getTypeWest(), 0.);
    pCorrFieldBcs.changeType(2, pFieldBcs.getTypeNorth(), 0.);
    pCorrFieldBcs.changeType(3, pFieldBcs.getTypeSouth(), 0.);
    pCorrFieldBcs.changeType(4, pFieldBcs.getTypeTop(), 0.);
    pCorrFieldBcs.changeType(5, pFieldBcs.getTypeBottom(), 0.);
}

void SimpleBoundaryCondition::setParallelBoundaries(std::shared_ptr<std::array<int, 6> > adjProcNoPtr)
{
    if(!adjProcNoPtr)
        return;

    FlowBoundaryConditions::setParallelBoundaries(adjProcNoPtr);
    pCorrFieldBcs.setParallelBoundaries(adjProcNoPtr);
}

bool SimpleBoundaryCondition::massCorrectionRequiredEast() const
{
    return types_[0] == OUTLET || types_[0] == PARALLEL;
}

bool SimpleBoundaryCondition::massCorrectionRequiredWest() const
{
    return types_[1] == OUTLET || types_[1] == PARALLEL;
}

bool SimpleBoundaryCondition::massCorrectionRequiredNorth() const
{
    return types_[2] == OUTLET || types_[2] == PARALLEL;
}

bool SimpleBoundaryCondition::massCorrectionRequiredSouth() const
{
    return types_[3] == OUTLET || types_[3] == PARALLEL;
}

bool SimpleBoundaryCondition::massCorrectionRequiredTop() const
{
    return types_[4] == OUTLET || types_[4] == PARALLEL;
}

bool SimpleBoundaryCondition::massCorrectionRequiredBottom() const
{
    return types_[5] == OUTLET || types_[5] == PARALLEL;
}

