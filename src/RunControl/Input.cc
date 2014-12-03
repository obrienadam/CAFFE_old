#include <boost/algorithm/string/erase.hpp>

#include "Input.h"
#include "InputStringProcessing.h"

Input::Input()
{
    //- Load all of the default input values, which may be overridden later

    // Time input data

    inputDoubles["simStartTime"] = 0.;
    inputDoubles["maxSimTime"] = 1.;
    inputDoubles["timeStep"] = 1e-4;
    inputDoubles["maxCpuTime"] = 100;
    inputDoubles["maxRealTime"] = 48;

    inputStrings["solver"] = "Euler";
    inputStrings["simTimeUnits"] = "seconds";
    inputStrings["cpuTimeUnits"] = "years";
    inputStrings["reaTimeUnits"] = "hours";

    // Iteration related data

    inputInts["maxItrs"] = 30000;

    // Run control related data

    inputStrings["terminationCondition"] = "simTime";
    inputInts["fileWriteInterval"] = 50;
    inputInts["screenWriteInterval"] = 50;

    // Domain related data

    inputInts["nI"] = 30;
    inputInts["nJ"] = 30;
    inputInts["nK"] = 30;

    // Scheme related data

    inputStrings["scheme"] = "finiteDifference";

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

    // Open an input file

    if(fin_.is_open())
        fin_.close();

    filename_ = filename;
    fin_.open(filename.c_str());

    //- Check if the file was found

    if(!fin_.is_open())
    {

        string errorMessage("Input file \"" + filename_ + "\" was not found.");

        throw errorMessage.c_str();

    }

    while(!fin_.eof())
    {

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
