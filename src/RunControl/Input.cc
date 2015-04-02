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

Input::Input()
{
    inputDoubles["simStartTime"] = 0.;                  ///< Start time in simulation time
    inputDoubles["maxSimTime"] = 1.;                    ///< End time in simulation time
    inputDoubles["timeStep"] = 1e-4;                    ///< Fixed time step
    inputDoubles["maxCpuTime"] = 100;                   ///< Maximum allowed CPU time
    inputDoubles["maxRealTime"] = 48;                   ///< Maximum allowed real time

    inputStrings["solver"] = "Euler";                   ///< The solver type
    inputStrings["simTimeUnits"] = "seconds";           ///< Simulation time units
    inputStrings["cpuTimeUnits"] = "years";             ///< CPU time units
    inputStrings["reaTimeUnits"] = "hours";             ///< Real time units

    inputInts["maxItrs"] = 30000;                       ///< Maximum allowed iterations

    inputStrings["terminationCondition"] = "simTime";   ///< The simulation termination condition
    inputInts["fileWriteInterval"] = 50;                ///< Number of iterations per file write
    inputInts["screenWriteInterval"] = 50;              ///< Number of interations before console output

    inputStrings["domainFile"] = "domain.in";           ///< Name of the domain input file
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

    string buffer;

    if(fin_.is_open())
        fin_.close();

    filename_ = filename;
    fin_.open(filename.c_str());

    if(!fin_.is_open())
    {
        throw ("Input file \"" + filename_ + "\" was not found.").c_str();
    }

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
            if(!(fin_ >> inputInts[buffer]))
                throw "A problem occurred while attempting to read an input integer.";
        }

        else if (inputDoubles.find(buffer) != inputDoubles.end())
        {
            if(!(fin_ >> inputDoubles[buffer]))
                throw "A problem occurred while attempting to read an input double.";
        }

        else if (inputStrings.find(buffer) != inputStrings.end())
        {
            if(!(fin_ >> inputStrings[buffer]))
                throw "A problem occurred while attempting to read an input string.";
        }
        else
        {
            throw ("Unrecognized input parameter type \"" + buffer + "\".").c_str();
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
