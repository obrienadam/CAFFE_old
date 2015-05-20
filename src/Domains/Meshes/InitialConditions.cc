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

#include <vector>

#include "InitialConditions.h"
#include "Output.h"
#include "InputStringProcessing.h"

InitialConditions::InitialConditions()
    :
      meshPtr_(NULL),
      metricConversion_(1.)
{

}

void InitialConditions::findOpeningBrace()
{
    using namespace std;

    string buffer;

    while(!inputFile_.eof())
    {
        getline(inputFile_, buffer);

        InputStringProcessing::processBuffer(buffer, true);

        if(buffer.empty())
            continue;
        else if(buffer == "{")
            break;
        else
            Output::raiseException("InitialConditions", "setInitialConditions", "received \"" + buffer + "\" but expected \"{\".");
    }
}

void InitialConditions::initialize(HexaFvmMesh &mesh)
{
    meshPtr_ = &mesh;

    nCellsI_ = mesh.nCellsI();
    nCellsJ_ = mesh.nCellsJ();
    nCellsK_ = mesh.nCellsK();
}

void InitialConditions::readInputFile(std::string filename)
{
    using namespace std;

    HexaFvmMesh& mesh = *meshPtr_;
    string buffer;

    inputFile_.open(filename.c_str());

    if(!inputFile_.good())
        Output::raiseException("InitialConditions", "openInputFile", "file \"" + filename + "\" does not exist.");

    while(!inputFile_.eof())
    {
        getline(inputFile_, buffer);

        InputStringProcessing::processBuffer(buffer, true);

        if(buffer.empty())
            continue;
        else if(buffer.substr(0, buffer.find_first_of("=")) == "MetricConversion")
        {
            metricConversion_ = std::stod(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            continue;
        }

        //- This is a bit shady. May want to look at improving this
        try
        {
            setInitialConditions(mesh.findScalarField(buffer));
        }
        catch(const char* errorMessage)
        {
            setInitialConditions(mesh.findVectorField(buffer));
        }
    }

    inputFile_.close();
}

void InitialConditions::setInitialConditions(Field<double> &scalarField)
{
    std::string buffer;

    findOpeningBrace();

    while(!inputFile_.eof())
    {
        getline(inputFile_, buffer);

        InputStringProcessing::processBuffer(buffer, true);

        if(buffer.empty())
            continue;
        else if(buffer == "Sphere")
            setSphere(scalarField);
        else if(buffer == "Box")
            setBox(scalarField);
        else if(buffer == "Uniform")
            setUniform(scalarField);
        else if(buffer == "}")
            break;
        else
            Output::raiseException("InitialConditions", "setInitialConditions", "unrecognized initial condition \"" + buffer + "\"");
    }
}

void InitialConditions::setInitialConditions(Field<Vector3D> &vectorField)
{

}

void InitialConditions::setSphere(Field<double> &scalarField)
{
    using namespace std;

    string buffer;
    vector<string> partitionedBuffer;
    double radius = 0., sphereInnerValue = 0.;
    Point3D center;
    bool radiusSet = false, sphereInnerValueSet = false, centerSet = false;

    findOpeningBrace();

    while(!inputFile_.eof())
    {
        getline(inputFile_, buffer);

        InputStringProcessing::processBuffer(buffer, true);

        if(buffer.empty())
            continue;
        else if(buffer.substr(0, buffer.find_first_of("=")) == "radius")
        {
            radius = stod(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            radiusSet = true;
        }
        else if(buffer.substr(0, buffer.find_first_of("=")) == "value")
        {
            sphereInnerValue = stod(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            sphereInnerValueSet = true;
        }
        else if(buffer.substr(0, buffer.find_first_of("=")) == "center")
        {
            buffer = buffer.substr(buffer.find_first_of("=") + 1, buffer.length());
            partitionedBuffer = InputStringProcessing::partition(buffer, "(,)");
            center.x = stod(partitionedBuffer[1]);
            center.y = stod(partitionedBuffer[2]);
            center.z = stod(partitionedBuffer[3]);
            centerSet = true;
        }
        else if(buffer == "}")
            break;
        else
            Output::raiseException("InitialConditions", "setSphere", "unrecognized initial condition \"" + buffer + "\"");
    }

    if(!(radiusSet && sphereInnerValueSet && centerSet))
        Output::raiseException("InitialConditions", "setSphere", "one or more required parameters not specified.");

    radius *= metricConversion_;
    center *= metricConversion_;

    createSphere(radius, center, sphereInnerValue, scalarField);
    Output::print("InitialConditions", "set spherical initial conditions for field \"" + scalarField.name + "\".");
}

void InitialConditions::setBox(Field<double> &scalarField)
{

}

void InitialConditions::setUniform(Field<double> &scalarField)
{

}

void InitialConditions::setSphere(Field<Vector3D> &scalarField)
{

}

void InitialConditions::setBox(Field<Vector3D> &scalarField)
{

}

void InitialConditions::setUniform(Field<Vector3D> &scalarField)
{

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

void InitialConditions::createSphere(double radius, Point3D center, double sphereInnerValue, Field<double> &scalarField)
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
                    scalarField(i, j, k) = sphereInnerValue;
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
