/**
 * @file    InitialConditionsI.h
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
 * This file contains the templated implementations for the methods of class
 * InitialConditions.
 */

#include <vector>

#include "InitialConditions.h"
#include "Output.h"

template <class T>
void InitialConditions::createUniform(T value, Field<T>& field)
{
    int i, j, k;
    const HexaFvmMesh &mesh = field.getMesh();

    Output::print("InitialConditions", "Setting uniform initial conditions for field \"" + field.name + "\".");

    for(k = 0; k < mesh.nCellsK(); ++k)
    {
        for(j = 0; j < mesh.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh.nCellsI(); ++i)
            {
                field(i, j, k) = value;
            }
        }
    }
}

template <class T>
void InitialConditions::createSphere(double radius, Point3D center, T sphereInnerValue, Field<T> &field)
{
    int i, j, k;
    const HexaFvmMesh &mesh = field.getMesh();

    for(k = 0; k < mesh.nCellsK(); ++k)
    {
        for(j = 0; j < mesh.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh.nCellsI(); ++i)
            {
                if((mesh.cellXc(i, j, k) - center).mag() <= radius)
                {
                    field(i, j, k) = sphereInnerValue;
                }
            }
        }
    }
}

template <class T>
void InitialConditions::createBox(double xLength, double yLength, double zLength, Point3D center, T boxInnerValue, Field<T>& field)
{
    int i, j, k;
    const HexaFvmMesh &mesh = field.getMesh();

    for(k = 0; k < mesh.nCellsK(); ++k)
    {
        for(j = 0; j < mesh.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh.nCellsI(); ++i)
            {
                if(mesh.cellXc(i, j, k).x >= center.x - 0.5*xLength && mesh.cellXc(i, j, k).x <= center.x + 0.5*xLength
                        && mesh.cellXc(i, j, k).y >= center.y - 0.5*yLength && mesh.cellXc(i, j, k).y <= center.y + 0.5*yLength
                        && mesh.cellXc(i, j, k).z >= center.z - 0.5*zLength && mesh.cellXc(i, j, k).z <= center.z + 0.5*zLength)
                {
                    field(i, j, k) = boxInnerValue;
                }
            }
        }
    }
}
