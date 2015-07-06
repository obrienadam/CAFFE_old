/**
 * @file    Input.cc
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
 * This file contains all of the implementations for the methods of
 * class Input.
 */

#include <iostream>

#include <boost/algorithm/string/erase.hpp>

#include "Input.h"
#include "InputStringProcessing.h"
#include "Output.h"

Input::Input()
{
    //- Set-up default input values
    inputDoubles["simStartTime"] = 0.;
    inputDoubles["maxSimTime"] = 1.;
    inputDoubles["timeStep"] = 1e-4;
    inputDoubles["maxRealTimeHours"] = 48.;
    inputDoubles["maxRealTimeMinutes"] = 0.;
    inputDoubles["maxRealTimeSeconds"] = 0.;

    inputInts["maxItrs"] = 1000;
    inputInts["fileWriteInterval"] = 50;
    inputInts["screenWriteInterval"] = 50;

    inputStrings["solver"] = "Euler";
    inputStrings["terminationCondition"] = "simTime";

    inputStrings["boundaryTypeEast"] = "fixed";
    inputStrings["boundaryTypeWest"] = "fixed";
    inputStrings["boundaryTypeNorth"] = "fixed";
    inputStrings["boundaryTypeSouth"] = "fixed";
    inputStrings["boundaryTypeTop"] = "fixed";
    inputStrings["boundaryTypeBottom"] = "fixed";

    inputStrings["boundaryRefValueEast"] = "(0., 0., 0.)";
    inputStrings["boundaryRefValueWest"] = "(0., 0., 0.)";
    inputStrings["boundaryRefValueNorth"] = "(0., 0., 0.)";
    inputStrings["boundaryRefValueSouth"] = "(0., 0., 0.)";
    inputStrings["boundaryRefValueTop"] = "(0., 0., 0.)";
    inputStrings["boundaryRefValueBottom"] = "(0., 0., 0.)";

    //- Linear solver parameters
    inputInts["kspMaxItrs"] = 5000;
    inputDoubles["kspRelativeTolerance"] = 1e-12;
    inputDoubles["kspAbsoluteTolerance"] = 1e-6;

    //- Simple parameters
    inputStrings["timeAccurate"] = "OFF";
    inputDoubles["relaxationFactorMomentum"] = 0.3;
    inputDoubles["relaxationFactorPCorr"] = 0.1;
    inputDoubles["rho"] = 998.;
    inputDoubles["mu"] = 8.94e-4;
    inputInts["maxInnerIters"] = 10;

    //- Piso parameters
    inputInts["numberOfPCorrections"] = 2;

    //- multiphaseSimple parameters
    inputDoubles["sigma"] = 0.07;
    inputDoubles["rho2"] = 782.;
    inputDoubles["mu2"] = 8.94e-3;

    //- ibSimple parameters
    inputDoubles["ibSphereRadius"] = 0.5;
    inputStrings["ibSphereCenter"] = "(0.5, 0.5, 0.5)";
}

Input::Input(std::string filename)
    :
      Input()
{
    openInputFile(filename);
}

Input::~Input()
{
    if(fin_.is_open())
        fin_.close();
}

void Input::openInputFile(std::string filename)
{
    using namespace std;

    string buffer, buffer2;

    if(fin_.is_open())
        fin_.close();

    filename_ = filename;
    fin_.open(filename.c_str());

    if(!fin_.is_open())
        Output::raiseException("Input", "openInputFile", "Input file \"" + filename_ + "\" was not found.");

    while(!fin_.eof())
    {
        getline(fin_, buffer);

        // Process the buffer

        buffer = InputStringProcessing::processBuffer(buffer);

        // Make sure the buffer is not empty, in case the last line was blank or a comment

        if(buffer.empty())
            continue;

        if(inputInts.find(buffer) != inputInts.end())
        {
            getline(fin_, buffer2);
            buffer2 = InputStringProcessing::processBuffer(buffer2);

            try
            {
                inputInts[buffer] = stoi(buffer2);
            }
            catch(...)
            {
                Output::raiseException("Input", "openInputFile", "Input field \"" + buffer + "\" expects an integer value.");
            }
        }

        else if (inputDoubles.find(buffer) != inputDoubles.end())
        {
            getline(fin_, buffer2);
            buffer2 = InputStringProcessing::processBuffer(buffer2);

            try
            {
                inputDoubles[buffer] = stod(buffer2);
            }
            catch(...)
            {
                Output::raiseException("Input", "openInputFile", "Input field \"" + buffer + "\" expects a double value.");
            }
        }

        else if (inputStrings.find(buffer) != inputStrings.end())
        {
            getline(fin_, buffer2);
            buffer2 = InputStringProcessing::processBuffer(buffer2, false);

            inputStrings[buffer] = buffer2;
        }
        else
        {
            Output::raiseException("Input", "openInputFile", "Unrecognized input field \"" + buffer + "\".");
        }
    }
}

void Input::print()
{
    using namespace std;

    map<string, int>::iterator intMapItr;
    map<string, double>::iterator doubleMapItr;
    map<string, string>::iterator stringMapItr;

    cout << "Input Integers:\n";

    for(intMapItr = inputInts.begin(); intMapItr != inputInts.end(); ++intMapItr)
    {
        cout << intMapItr->first << ": " << intMapItr->second << endl;
    }

    cout << "Input Doubles:\n";

    for(doubleMapItr = inputDoubles.begin(); doubleMapItr != inputDoubles.end(); ++doubleMapItr)
    {
        cout << doubleMapItr->first << ": " << doubleMapItr->second << endl;
    }

    cout << "Input Strings:\n";

    for(stringMapItr = inputStrings.begin(); stringMapItr != inputStrings.end(); ++stringMapItr)
    {
        cout << stringMapItr->first << ": " << stringMapItr->second << endl;
    }
}
