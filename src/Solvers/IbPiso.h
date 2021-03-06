/**
 * @file    IbPiso.h
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
 * This file contains the interface for class IbSimple, which contains
 * schemes for using the SIMPLE method for solving incompressible flow
 * problems, in conjunction with the immersed boundary method for
 * handling more complex geometries.
 */

#ifndef IB_PISO_H
#define IB_PISO_H

#include "Piso.h"
#include "Sphere.h"

class IbPiso : public Piso
{

public:

    IbPiso(const Input &input, const HexaFvmMesh &mesh);

    virtual double solve(double timeStep);

private:

    enum CellType{FLUID, IB, SOLID};

    void computeIbField();
    void computeIbCoeffs();

    Field<CellType> ibField_;
    Sphere ibSphere_;
};

#endif
