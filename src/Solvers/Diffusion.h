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

#include "PrimitiveBoundaryCondition.h"
#include "Solver.h"
#include "Field.h"

class Diffusion : public Solver
{
public:

    Diffusion(const Input &input, const HexaFvmMesh &mesh);
    ~Diffusion();

    virtual double solve(double timeStep);

protected:

    Field<double> phiField_;
    Field<Vector3D> gradPhiField_;

    PrimitiveBoundaryCondition<double> bcs_;
};

#endif
