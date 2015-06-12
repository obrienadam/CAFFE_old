/**
 * @file    MultiphaseSimple.h
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
 * This file contains the interface for class MultiphaseSimple, which contains
 * schemes for using the SIMPLE method for solving multi-phase incompressible
 * flow problems.
 */

#ifndef MULTIPHASE_SIMPLE_H
#define MULTIPHASE_SIMPLE_H

#include "Simple.h"
#include "Vector3D.h"

class MultiphaseSimple : public Simple
{
private:

    Field<double>* alphaFieldPtr_;
    Field<double> alphaField0_;
    Field<Vector3D> gradAlphaField_;
    Field<Vector3D> interfaceNormals_;
    Field<double> kField_;
    Field<Vector3D> bF_;

    double rho1_, rho2_, mu1_, mu2_;

    double sigma_;

    int alphaGmresIters_;

    void computePhysicalConstants(Field<double>& alphaField, Field<double>& rhoField, Field<double>& muField);

    void computeCurvature(Field<double>& alphaField);
    void computeSurfaceTensionSource();
    void advectAlphaField(Field<double> &rhoField, Field<double> &massFlowField, Field<Vector3D>& uField, double timeStep, Field<double>& alphaField);

public:

    MultiphaseSimple();

    void initialize(Input& input, HexaFvmMesh& mesh);

    void storeAlphaField(Field<double>& alphaField);

    void discretize(double timeStep, std::vector<double>& timeDerivatives);
};

#endif
