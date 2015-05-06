/**
 * @file    FvScheme.h
 * @author  Adam O'Brien <obrienadam89@gmail.com>
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * https://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * This file contains the interface for class FvScheme, which is an
 * abstract interface class used for deriving specific finite volume
 * spatial schemes.
 */

#ifndef FV_SCHEME_H
#define FV_SCHEME_H

#include <vector>
#include <string>

#include "Input.h"
#include "Vector3D.h"
#include "Tensor3D.h"
#include "Field.h"
#include "HexaFvmMesh.h"

enum {ADD, REPLACE};
enum {LEAST_SQUARES, DIVERGENCE_THEOREM};
enum {VOLUME_WEIGHTED, DISTANCE_WEIGHTED, NON_WEIGHTED};

class FvScheme
{
protected:

    std::string conservedFieldName_;
    HexaFvmMesh* meshPtr_;
    int nCellsI_, nCellsJ_, nCellsK_, nFacesI_, nFacesJ_, nFacesK_;

public:

    FvScheme();

    virtual void initialize(Input& input, HexaFvmMesh& mesh, std::string conservedFieldName = "phi");
    virtual int nConservedVariables() = 0;

    virtual void discretize(std::vector<double>& timeDerivatives_) = 0;
    virtual void copySolution(std::vector<double>& original) = 0;
    virtual void updateSolution(std::vector<double>& timeDerivatives_, int method) = 0;

    /**
     * @brief This method is used for computing a weighted averaging coefficient based on a specified criteria.
     * @return An interpolation factor.
     */
    double getAlpha(int i, int j, int k, int direction);

    template <class FIELD>
    void interpolateInteriorFaces(FIELD& field, int method = NON_WEIGHTED);

    /**
     * @brief Compute the gradient of a scalar field at the cell center using a least-squares reconstruction method.
     * @param phiField A reference to the scalar field.
     * @param gradPhiField A reference to the vector field that will contain the cell-centered gradients.
     */
    virtual void computeCellCenteredGradients(Field<double>& phiField, Field<Vector3D>& gradPhiField, int method = LEAST_SQUARES);

    /**
     * @brief Compute the Jacobian of a vector field at the cell center using a least-squares reconstruction method.
     * @param vecField A reference to a vector field.
     * @param tensorField A reference to a tensor field that will contain the computed jacobians.
     */
    virtual void computeCellCenteredJacobians(Field<Vector3D>& vecField, Field<Tensor3D>& tensorField, int method = LEAST_SQUARES);

    /**
     * @brief Compute the gradient of a scalar field at the face center. An apropriate cell-centered gradient computation method should be called first.
     * @param phiField A reference to the scalar field.
     * @param gradPhiField A reference to the vector field that contains the computed cell-centered gradients.
     */
    virtual void computeFaceCenteredGradients(Field<double>& phiField, Field<Vector3D>& gradPhiField);

    void getMeshStencil(int i, int j, int k, int direction, Vector3D& faceNorm, Vector3D& cellRelVec, double &alpha);

    virtual void displayUpdateMessage();
};

#include "FvSchemeI.h"

#endif
