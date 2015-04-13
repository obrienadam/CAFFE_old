/**
 * @file    FvScheme.cc
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
 * This file contains the implementations of methods for class FvScheme.
 */

#include <cstdlib>
#include <math.h>

#include "FvScheme.h"
#include "Output.h"

FvScheme::FvScheme()
    :
      conservedFieldName_("phi"),
      meshPtr_(NULL)
{

}

void FvScheme::initialize(Input &input, HexaFvmMesh &mesh, std::string conservedFieldName)
{
    meshPtr_ = &mesh;
    conservedFieldName_ = conservedFieldName;
}

double FvScheme::getAlpha(int i, int j, int k, int direction)
{
    int deltaI = 0, deltaJ = 0, deltaK = 0;

    switch (direction)
    {
    case EAST:
        ++deltaI;
        break;
    case WEST:
        --deltaI;
        break;
    case NORTH:
        ++deltaJ;
        break;
    case SOUTH:
        --deltaJ;
        break;
    case TOP:
        ++deltaK;
        break;
    case BOTTOM:
        --deltaK;
        break;
    default:
        Output::raiseException("FvScheme", "alpha", "Invalid direction specified.");
    };

    return meshPtr_->cellVol(i, j, k)/(meshPtr_->cellVol(i, j, k) + meshPtr_->cellVol(i + deltaI, j + deltaJ, k + deltaK));
}
