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
    int uCellI_, uCellJ_, uCellK_, uFaceI_, uFaceJ_, uFaceK_;

public:

    FvScheme();

    virtual void initialize(Input& input, HexaFvmMesh& mesh, std::string conservedFieldName = "phi");
    virtual int nConservedVariables() = 0;

    virtual void discretize(std::vector<double>& timeDerivatives_) = 0;
    virtual void copySolution(std::vector<double>& original) = 0;
    virtual void updateSolution(std::vector<double>& timeDerivatives_, int method) = 0;

    template <class T>
    void interpolateInteriorFaces(Field<T>& field, int method);

    template <class T, class GRAD_T>
    void extrapolateInteriorFaces(Field<T>& field, Field<GRAD_T>& gradField);

    /**
     * @brief Compute the gradient of a scalar field at the cell center.
     * @param phiField A reference to the scalar field.
     * @param gradPhiField A reference to the vector field that will contain the cell-centered gradients.
     */
    void computeCellCenteredGradients(Field<double>& phiField, Field<Vector3D>& gradPhiField, int method);

    /**
     * @brief Compute the gradient (Jacobian) of a vector field at the cell center.
     * @param vecField A reference to the scalar field.
     * @param gradVecField A reference to the tensor field that will contain the cell-centered gradients (Jacobians).
     */
    void computeCellCenteredGradients(Field<Vector3D>& vecField, Field<Tensor3D>& gradVecField, int method);

    /**
     * @brief Compute the gradient of a scalar field at the face center. An apropriate cell-centered gradient computation method should be called first.
     * @param phiField A reference to the scalar field.
     * @param gradPhiField A reference to the vector field that contains the computed cell-centered gradients.
     */
    virtual void computeFaceCenteredGradients(Field<double>& phiField, Field<Vector3D>& gradPhiField);

    void getMeshStencil(int i, int j, int k, int direction, Vector3D& faceNorm, Vector3D& cellRelVec, double &alpha);

    Vector3D sfE(int i, int j, int k){ return meshPtr_->fAreaNormE(i, j, k); }
    Vector3D sfW(int i, int j, int k){ return meshPtr_->fAreaNormW(i, j, k); }
    Vector3D sfN(int i, int j, int k){ return meshPtr_->fAreaNormN(i, j, k); }
    Vector3D sfS(int i, int j, int k){ return meshPtr_->fAreaNormS(i, j, k); }
    Vector3D sfT(int i, int j, int k){ return meshPtr_->fAreaNormT(i, j, k); }
    Vector3D sfB(int i, int j, int k){ return meshPtr_->fAreaNormB(i, j, k); }

    Vector3D dsE(int i, int j, int k){ return meshPtr_->rCellE(i, j, k); }
    Vector3D dsW(int i, int j, int k){ return meshPtr_->rCellW(i, j, k); }
    Vector3D dsN(int i, int j, int k){ return meshPtr_->rCellN(i, j, k); }
    Vector3D dsS(int i, int j, int k){ return meshPtr_->rCellS(i, j, k); }
    Vector3D dsT(int i, int j, int k){ return meshPtr_->rCellT(i, j, k); }
    Vector3D dsB(int i, int j, int k){ return meshPtr_->rCellB(i, j, k); }

    virtual void displayUpdateMessage();
};

#include "FvSchemeI.h"

#endif
