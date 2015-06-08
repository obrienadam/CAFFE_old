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
 * This file contains the implementations for the methods of class
 * InitialConditions.
 */

#include "InitialConditions.h"

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

template <>
void InitialConditions::setUniform(Field<double>& field)
{
    using namespace std;

    string buffer;
    double value;
    bool valueSet = false;

    findOpeningBrace();

    while(!inputFile_.eof())
    {
        getline(inputFile_, buffer);
        InputStringProcessing::processBuffer(buffer, true);

        if(buffer.empty())
            continue;
        else if(buffer.substr(0, buffer.find_first_of("=")) == "value")
        {
            value = stod(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            valueSet = true;
        }
    }

    if(!valueSet)
        Output::raiseException("InitialConditions", "setUniform", "uniform field value not set.");

    createUniform(value, field);
}

template <>
void InitialConditions::setUniform(Field<Vector3D>& field)
{
    using namespace std;

    string buffer;
    Vector3D value;
    bool valueSet = false;

    findOpeningBrace();

    while(!inputFile_.eof())
    {
        getline(inputFile_, buffer);
        InputStringProcessing::processBuffer(buffer, true);

        if(buffer.empty())
            continue;
        else if(buffer.substr(0, buffer.find_first_of("=")) == "value")
        {
            value = stov(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            valueSet = true;
        }
    }

    if(!valueSet)
        Output::raiseException("InitialConditions", "setUniform", "uniform field value not set.");

    createUniform(value, field);
    Output::print("InitialConditions", "set uniform initial conditions for field \"" + field.name + "\".");
}

template <>
void InitialConditions::setSphere(Field<double>& field)
{
    using namespace std;

    string buffer;
    double radius = 0.;
    double sphereInnerValue = 0.;
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
            center = stov(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
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

    createSphere(radius, center, sphereInnerValue, field);
    Output::print("InitialConditions", "set spherical initial conditions for field \"" + field.name + "\".");
}

template <>
void InitialConditions::setSphere(Field<Vector3D>& field)
{
    using namespace std;

    string buffer;
    double radius = 0.;
    Vector3D sphereInnerValue;
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
            sphereInnerValue = stov(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            sphereInnerValueSet = true;
        }
        else if(buffer.substr(0, buffer.find_first_of("=")) == "center")
        {
            center = stov(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
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

    createSphere(radius, center, sphereInnerValue, field);
    Output::print("InitialConditions", "set spherical initial conditions for field \"" + field.name + "\".");
}

template <>
void InitialConditions::setBox(Field<double>& field)
{
    using namespace std;

    string buffer;
    double xLength = 0., yLength = 0., zLength = 0.;
    double boxInnerValue = 0.;
    Point3D center;
    bool xLengthSet = false, yLengthSet = false, zLengthSet = false, centerSet = false, boxInnerValueSet = false;

    findOpeningBrace();

    while(!inputFile_.eof())
    {
        getline(inputFile_, buffer);
        InputStringProcessing::processBuffer(buffer, true);

        if(buffer.empty())
            continue;
        else if(buffer.substr(0, buffer.find_first_of("=")) == "xLength")
        {
            xLength = stod(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            xLengthSet = true;
        }
        else if(buffer.substr(0, buffer.find_first_of("=")) == "yLength")
        {
            yLength = stod(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            yLengthSet = true;
        }
        else if(buffer.substr(0, buffer.find_first_of("=")) == "zLength")
        {
            zLength = stod(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            zLengthSet = true;
        }
        else if(buffer.substr(0, buffer.find_first_of("=")) == "center")
        {
            center = stov(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            centerSet = true;
        }
        else if(buffer.substr(0, buffer.find_first_of("=")) == "boxInnerValue")
        {
            boxInnerValue = stod(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            boxInnerValueSet = true;
        }
        else if(buffer == "}")
            break;
        else
            Output::raiseException("InitialConditions", "setBox", "unrecognized initial condition \"" + buffer + "\"");
    }

    if(!(xLengthSet && yLengthSet && zLengthSet && centerSet && boxInnerValueSet))
        Output::raiseException("InitialConditions", "setBox", "one or more required parameters not specified.");

    xLength *= metricConversion_;
    yLength *= metricConversion_;
    zLength *= metricConversion_;
    center *= metricConversion_;

    createBox(xLength, yLength, zLength, center, boxInnerValue, field);
    Output::print("InitialConditions", "set box initial conditions for field \"" + field.name + "\".");
}

template <>
void InitialConditions::setBox(Field<Vector3D>& field)
{
    using namespace std;

    string buffer;
    double xLength = 0., yLength = 0., zLength = 0.;
    Vector3D boxInnerValue = Vector3D(0., 0., 0.);
    Point3D center;
    bool xLengthSet = false, yLengthSet = false, zLengthSet = false, centerSet = false, boxInnerValueSet = false;

    findOpeningBrace();

    while(!inputFile_.eof())
    {
        getline(inputFile_, buffer);
        InputStringProcessing::processBuffer(buffer, true);

        if(buffer.empty())
            continue;
        else if(buffer.substr(0, buffer.find_first_of("=")) == "xLength")
        {
            xLength = stod(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            xLengthSet = true;
        }
        else if(buffer.substr(0, buffer.find_first_of("=")) == "yLength")
        {
            yLength = stod(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            yLengthSet = true;
        }
        else if(buffer.substr(0, buffer.find_first_of("=")) == "zLength")
        {
            zLength = stod(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            zLengthSet = true;
        }
        else if(buffer.substr(0, buffer.find_first_of("=")) == "center")
        {
            center = stov(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            centerSet = true;
        }
        else if(buffer.substr(0, buffer.find_first_of("=")) == "boxInnerValue")
        {
            boxInnerValue = stov(buffer.substr(buffer.find_first_of("=") + 1, buffer.length()));
            boxInnerValueSet = true;
        }
        else if(buffer == "}")
            break;
        else
            Output::raiseException("InitialConditions", "setBox", "unrecognized initial condition \"" + buffer + "\"");
    }

    if(!(xLengthSet && yLengthSet && zLengthSet && centerSet && boxInnerValueSet))
        Output::raiseException("InitialConditions", "setBox", "one or more required parameters not specified.");

    xLength *= metricConversion_;
    yLength *= metricConversion_;
    zLength *= metricConversion_;
    center *= metricConversion_;

    createBox(xLength, yLength, zLength, center, boxInnerValue, field);
    Output::print("InitialConditions", "set box initial conditions for field \"" + field.name + "\".");
}
