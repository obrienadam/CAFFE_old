/**
 * @file    IbSimple.h
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

#ifndef IB_SIMPLE_H
#define IB_SIMPLE_H

#include "Simple.h"
#include "Sphere.h"

enum CellType{FLUID, IB, SOLID};

class IbSimple : public Simple
{
private:

    Field<CellType> ibField_;
    Field<Vector3D> ibSourceField_;

    Sphere ibSphere_;

    void computeIbField(Field<Vector3D>& uField, Field<double>& pField);
    void setIbCells(Field<Vector3D> &uField, Field<double>& pField);

public:

    IbSimple();

    void initialize(Input& input, HexaFvmMesh& mesh);

    void discretize(double timeStep, std::vector<double>& timeDerivatives);
};

#endif
