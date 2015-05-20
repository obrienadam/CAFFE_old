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
    Field<Vector3D> interfaceNormals_;
    Field<Vector3D> kField_;

    double tau_;

public:

    void initialize(Input& input, HexaFvmMesh& mesh);
};

#endif
