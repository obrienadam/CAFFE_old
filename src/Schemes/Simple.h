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

class Simple : public FvScheme
{
private:

    Field<double> *uFieldPtr_, *vFieldPtr_, *wFieldPtr_, *pFieldPtr_;
    Field<Vector3D> gradUField_, gradVField_, gradWField_, gradPField_;

public:

    void initialize(Input &input, HexaFvmMesh &mesh);

    int nConservedVariables();

    void discretize(std::vector<double>& timeDerivatives_);
    void copySolution(std::vector<double>& original);
    void updateSolution(std::vector<double>& update, int method);
};

#endif
