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
    inputDoubles_["simStartTime"] = 0.;
    inputDoubles_["maxSimTime"] = 1.;
    inputDoubles_["timeStep"] = 1e-4;
    inputDoubles_["maxRealTimeHours"] = 48.;
    inputDoubles_["maxRealTimeMinutes"] = 0.;
    inputDoubles_["maxRealTimeSeconds"] = 0.;
    inputInts_["maxNumberOfIterations"] = 1000;
    inputInts_["fileWriteInterval"] = 50;
    inputStrings_["terminationCondition"] = "simTime";
    inputStrings_["timeAccurate"] = "OFF";

    //- Boundary conditions
    inputStrings_["boundaryTypeEast"] = "fixed";
    inputStrings_["boundaryTypeWest"] = "fixed";
    inputStrings_["boundaryTypeNorth"] = "fixed";
    inputStrings_["boundaryTypeSouth"] = "fixed";
    inputStrings_["boundaryTypeTop"] = "fixed";
    inputStrings_["boundaryTypeBottom"] = "fixed";

    inputStrings_["boundaryRefVectorEast"] = "(0., 0., 0.)";
    inputStrings_["boundaryRefVectorWest"] = "(0., 0., 0.)";
    inputStrings_["boundaryRefVectorNorth"] = "(0., 0., 0.)";
    inputStrings_["boundaryRefVectorSouth"] = "(0., 0., 0.)";
    inputStrings_["boundaryRefVectorTop"] = "(0., 0., 0.)";
    inputStrings_["boundaryRefVectorBottom"] = "(0., 0., 0.)";

    inputDoubles_["boundaryRefValueEast"] = 0.;
    inputDoubles_["boundaryRefValueWest"] = 0.;
    inputDoubles_["boundaryRefValueNorth"] = 0.;
    inputDoubles_["boundaryRefValueSouth"] = 0.;
    inputDoubles_["boundaryRefValueTop"] = 0.;
    inputDoubles_["boundaryRefValueBottom"] = 0.;

    //- Simple parameters
    inputDoubles_["relaxationFactorMomentum"] = 0.3;
    inputDoubles_["relaxationFactorPCorr"] = 0.1;
    inputDoubles_["rho"] = 998.;
    inputDoubles_["mu"] = 8.94e-4;
    inputInts_["numberOfInnerIterations"] = 1;

    //- Piso parameters
    inputInts_["numberOfPressureCorrections"] = 2;

    //- ibPiso parameters
    inputDoubles_["ibSphereRadius"] = 0.5;
    inputStrings_["ibSphereCenter"] = "(0.5, 0.5, 0.5)";
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

        if(inputInts_.find(buffer) != inputInts_.end())
        {
            getline(fin_, buffer2);
            buffer2 = InputStringProcessing::processBuffer(buffer2);

            try
            {
                inputInts_[buffer] = stoi(buffer2);
            }
            catch(...)
            {
                Output::raiseException("Input", "openInputFile", "Input field \"" + buffer + "\" expects an integer value.");
            }
        }

        else if (inputDoubles_.find(buffer) != inputDoubles_.end())
        {
            getline(fin_, buffer2);
            buffer2 = InputStringProcessing::processBuffer(buffer2);

            try
            {
                inputDoubles_[buffer] = stod(buffer2);
            }
            catch(...)
            {
                Output::raiseException("Input", "openInputFile", "Input field \"" + buffer + "\" expects a double value.");
            }
        }

        else if (inputStrings_.find(buffer) != inputStrings_.end())
        {
            getline(fin_, buffer2);
            buffer2 = InputStringProcessing::processBuffer(buffer2, false);

            inputStrings_[buffer] = buffer2;
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

    for(intMapItr = inputInts_.begin(); intMapItr != inputInts_.end(); ++intMapItr)
    {
        cout << intMapItr->first << ": " << intMapItr->second << endl;
    }

    cout << "Input Doubles:\n";

    for(doubleMapItr = inputDoubles_.begin(); doubleMapItr != inputDoubles_.end(); ++doubleMapItr)
    {
        cout << doubleMapItr->first << ": " << doubleMapItr->second << endl;
    }

    cout << "Input Strings:\n";

    for(stringMapItr = inputStrings_.begin(); stringMapItr != inputStrings_.end(); ++stringMapItr)
    {
        cout << stringMapItr->first << ": " << stringMapItr->second << endl;
    }
}
