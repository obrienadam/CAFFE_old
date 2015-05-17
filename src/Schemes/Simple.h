/**
 * @file    Simple.h
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
 * This file contains the interface for class Simple, which contains
 * schemes for using the SIMPLE method for solving incompressible flow
 * problems.
 */

#ifndef SIMPLE_H
#define SIMPLE_H

#include "FvScheme.h"
#include "Tensor3D.h"
#include "SparseVector.h"

enum FlowBoundary{INLET, OUTLET, WALL};

class Simple : public FvScheme
{
private:

    Field<Vector3D>* uFieldPtr_;
    Field<double>* pFieldPtr_;
    Field<Vector3D> uStar_;
    Field<double> a0P_, aP_, aE_, aW_, aN_, aS_, aT_, aB_;
    Field<Vector3D> bP_;
    Field<Vector3D> hField_;
    Field<double> dField_;

    Field<double> massFlow_;
    Field<double> pCorr_;
    Field<Vector3D> gradPCorr_;

    Field<Tensor3D> gradUField_;
    Field<Vector3D> gradPField_;

    Field<Tensor3D> gradVecField_;
    Field<Vector3D> gradScalarField_;

    bool timeAccurate_;

    double relaxationFactorMomentum_, relaxationFactorPCorr_, rho_, mu_, nu_;
    int gradReconstructionMethod_;

    int momentumSorItrs_, pCorrSorItrs_, maxMomentumSorIters_, maxPCorrSorIters_;
    double momentumSorToler_, pCorrSorToler_, momentumSorConvergence_, pCorrSorConvergence_, sorOmega_;
    Vector3D momentumResidual_;

    int uxStartI_, uyStartI_, uzStartI_, pStartI_;
    int uxEndI_, uyEndI_, uzEndI_, pEndI_;

    /**
     * @brief Store the solution from the previous time-step/iteration.
     * @param uField A reference to the velocity field.
     */
    void storeUStar(Field<Vector3D>& uField);

    /**
     * @brief Compute a predicted momentum using the latest available pressure field.
     */
    void computeMomentum(double timeStep, Field<Vector3D>& uField, Field<double>& pField);

    /**
     * @brief Interpolate the face centered velocities using the Rhie-Chow interpolation method
     * @param uField A reference to the velocity field.
     * @param pField A reference to the pressure field.
     * @param dField A reference to the D field.
     */
    void rhieChowInterpolateInteriorFaces(Field<Vector3D>& uField, Field<double>& pField);

    /**
     * @brief Compute the mass fluxes at the interior faces from the known velocity field.
     * @param uField A reference to the velocity field.
     */
    void computeMassFlowFaces(Field<Vector3D>& uField);

    /**
     * @brief Compute the pressure corrections.
     */
    void computePCorr(Field<Vector3D>& uField, Field<double> &pField);

    /**
     * @brief Correct the mass flow, pressure and velocity fields using the computed pressure corrections.
     * @param uField A reference to the velocity field.
     * @param pField A reference to the pressure field.
     */
    void correctContinuity(Field<Vector3D>& uField, Field<double>& pField);

    Vector3D computeResidual(Field<Vector3D>& uField);

public:

    Simple();

    void initialize(Input &input, HexaFvmMesh &mesh);
    void setBoundaryConditions(Input& input);

    int nConservedVariables();

    void storeUField(Field<Vector3D>& uField);

    void discretize(double timeStep, std::vector<double>& timeDerivatives);
    void copySolution(std::vector<double>& original);
    void updateSolution(std::vector<double>& update, int method);

    void displayUpdateMessage();
};

#endif
