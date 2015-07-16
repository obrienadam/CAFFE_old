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

#include "Solver.h"
#include "Field.h"

class Simple : public Solver
{
public:

    Simple(const Input &input, const HexaFvmMesh &mesh);

    virtual double solve(double timeStep);
    virtual void displayUpdateMessage();

public:

    void setBoundaries(const Input &input);
    void setConstantFields(const Input &input);

    void computeMomentum(double timeStep);
    void computePCorr();
    void correct();
    double computeContinuityError();
    void rhieChowInterpolateFaces();

    //- Primary fields
    Field<Vector3D> uField_;
    Field<double> pField_;
    Field<double> rhoField_;
    Field<double> muField_;
    Field<double> massFlowField_;
    Field<double> pCorrField_;

    //- Gradient fields
    Field<Tensor3D> gradUField_;
    Field<Vector3D> gradPField_;
    Field<Vector3D> gradPCorrField_;

    //- Misc auxillary fields
    Field<Vector3D> uField0_;
    Field<Vector3D> uFieldStar_;
    Field<double> dField_;
    Field<Vector3D> hField_;

    //- Misc gradient fields to be used for reconstructions of misc fields
    Field<Vector3D> gradScalarField_;
    Field<Tensor3D> gradVectorField_;

    int nInnerIters_;
    double omegaMomentum_, omegaPCorr_;
    double continuityError_;
};

#endif
