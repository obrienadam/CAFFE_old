/**
 * @file    Diffusion.h
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
 * This file contains the interface for class Diffusion, which is a
 * class containing methods for discretizing diffusion problems.
 */

#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "FvScheme.h"
#include "Array3D.h"
#include "Point3D.h"
#include "Vector3D.h"
#include "Matrix.h"

class Diffusion : public FvScheme
{
private:

    Field<double>* phiFieldPtr_;
    Field<Vector3D> gradPhiField_;
    Array3D<double> stencil_;

public:
    /** @brief Helper function that computes the gradient at cell centers using a least squares method.
     */
    void computeCellCenteredGradients();

    /** @brief Helper function that computes the gradient at cell faces using the method from Mathur and Murthy.
     */
    void computeFaceCenteredGradients();

    /**
     * @brief Helper function that computes all of the face fluxes after the reconstructions are complete.
     */
    void computeFaceFluxes();

public:

    Diffusion();
    ~Diffusion();

    void initialize(HexaFvmMesh &mesh, std::string conservedFieldName);
    int nConservedVariables();

    void discretize(std::vector<double>& timeDerivatives);
    void copySolution(std::vector<double>& original);
    void updateSolution(std::vector<double>& update, int method);
};

#endif
