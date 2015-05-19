/**
 * @file    InitialConditions.cc
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
 * This file contains the implementations for the methods of class
 * InitialConditions.
 */

#include "InitialConditions.h"
#include "Output.h"

InitialConditions::InitialConditions()
    :
      meshPtr_(NULL)
{

}

void InitialConditions::initialize(HexaFvmMesh &mesh)
{
    meshPtr_ = &mesh;

    nCellsI_ = mesh.nCellsI();
    nCellsJ_ = mesh.nCellsJ();
    nCellsK_ = mesh.nCellsK();
}

void InitialConditions::openInputFile(std::string filename, std::string directory)
{
    inputFile_.open((directory + "/" + filename).c_str());

    if(!inputFile_.good())
        Output::raiseException("InitialConditions", "openInputFile", "file \"" + directory + "/" + filename + "\" does not exist.");
}

void InitialConditions::createUniform(double value, Field<double> &phiField)
{
    int i, j, k;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                phiField(i, j, k) = value;
            }
        }
    }
}

void InitialConditions::createSphere(double radius, Point3D center, double sphereInnerValue, Field<double> &phiField)
{
    HexaFvmMesh& mesh = *meshPtr_;
    int i, j, k;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if((mesh.cellXc(i, j, k) - center).mag() <= radius)
                {
                    phiField(i, j, k) = sphereInnerValue;
                }
            }
        }
    }
}

void InitialConditions::createBox(double xLength, double yLength, double zLength, Point3D center, double boxInnerValue, Field<double> &phiField)
{
    HexaFvmMesh& mesh = *meshPtr_;
    int i, j, k;

    for(k = 0; k < nCellsK_; ++k)
    {
        for(j = 0; j < nCellsJ_; ++j)
        {
            for(i = 0; i < nCellsI_; ++i)
            {
                if(mesh.cellXc(i, j, k).x >= center.x - 0.5*xLength && mesh.cellXc(i, j, k).x <= center.x + 0.5*xLength &&
                        mesh.cellXc(i, j, k).y >= center.y - 0.5*yLength && mesh.cellXc(i, j, k).y <= center.y + 0.5*yLength &&
                        mesh.cellXc(i, j, k).z >= center.z - 0.5*zLength && mesh.cellXc(i, j, k).z <= center.z + 0.5*zLength)
                {
                    phiField(i, j, k) = boxInnerValue;
                }
            }
        }
    }
}
