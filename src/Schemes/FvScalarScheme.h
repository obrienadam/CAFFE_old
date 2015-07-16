#ifndef FV_SCALAR_SCHEME_H
#define FV_SCALAR_SCHEME_H

#include "FvScheme.h"

class FvScalarScheme : public FvScheme
{
public:

    static void computeCellCenteredGradients(GradientEvaluationMethod method, const Field<double> &field, Field<Vector3D> &gradField);
    static void extrapolateInteriorFaces(GradientEvaluationMethod method, Field<double>& field, Field<Vector3D>& gradField);
};

#endif
