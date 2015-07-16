#ifndef FV_VECTOR_SCHEME_H
#define FV_VECTOR_SCHEME_H

#include "FvScheme.h"
#include "Vector3D.h"
#include "Tensor3D.h"

class FvVectorScheme : public FvScheme
{
public:

    static void computeCellCenteredGradients(GradientEvaluationMethod method, const Field<Vector3D> &field, Field<Tensor3D> &gradField);
    static void extrapolateInteriorFaces(GradientEvaluationMethod method, Field<Vector3D>& field, Field<Tensor3D>& gradField);
};

#endif
