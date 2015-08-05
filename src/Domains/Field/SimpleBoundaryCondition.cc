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
