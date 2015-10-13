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
#include "Kernel.h"

template <class T>
void InitialConditions::createUniform(T value, Field<T>& field)
{
    Output::print("InitialConditions", "Setting uniform initial conditions for field \"" + field.name + "\".");
    field.setInitialCondition([value](Point3D){ return value; });
}

template <class T>
void InitialConditions::createSphere(double radius, Point3D center, T sphereInnerValue, Field<T> &field)
{
    Output::print("InitialConditions", "Generating sphere at " + std::to_string(center) + " with radius " + std::to_string(radius) + " for field \"" + field.name + "\".");

    auto icFunc = [radius, center, sphereInnerValue](Point3D point)
    {
        return (point - center).mag() <= radius ? sphereInnerValue : T();
    };

    field.setInitialCondition(icFunc);

    smootheField(field);
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

template <class T>
void InitialConditions::smootheField(Field<T> &field)
{
    using namespace std;

    const HexaFvmMesh &mesh = field.getMesh();
    Field<T> smoothedField(mesh, Field<T>::AUXILLARY, "smoothedField");
    Kernel kernel(Kernel::K6, 2.*(mesh.faceXcI(0, 0, 0) - mesh.faceXcI(1, 0, 0)).mag());

    Output::print("InitialConditions", "Smoothing field \"" + field.name + "\".");

    for(int k = 0; k < mesh.nCellsK(); ++k)
    {
        for(int j = 0; j < mesh.nCellsJ(); ++j)
        {
            for(int i = 0; i < mesh.nCellsI(); ++i)
            {
                Point3D center = mesh.cellXc(i, j, k);

                for(int kk = max(k - 4, 0), endK = min(k + 4, mesh.uCellK()); kk < endK; ++kk)
                {
                    for(int jj = max(j - 4, 0), endJ = min(j + 4, mesh.uCellJ()); jj < endJ; ++jj)
                    {
                        for(int ii = max(i - 4, 0), endI = min(i + 4, mesh.uCellI()); i < endI; ++ii)
                        {
                            smoothedField(i, j, k) += field(ii, jj, kk)*kernel.value(center, mesh.cellXc(ii, jj, kk))*mesh.cellVol(ii, jj, kk);
                        }
                    }
                }
            }
        }
    }

    field = smoothedField;
}
