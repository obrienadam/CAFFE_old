#include <boost/algorithm/string/erase.hpp>
#include "Input.h"

Input::Input()
{
  //- Load all of the accepted file input parameters

  // Time input data

  inputDoubles["simStartTime="] = 0.;
  inputDoubles["simEndTime="] = 1.;
  inputDoubles["timeStep="] = 1e-4;
  inputDoubles["maxCpuTime="] = 100;
  inputDoubles["maxRealTime="] = 48;

  inputStrings["simTimeUnits="] = "seconds";
  inputStrings["cpuTimeUnits="] = "years";
  inputStrings["reaTimeUnits="] = "hours";

  // Iteration related data
  
  inputInts["maxItrs="] = 30000;
  
  // Run control related data

  inputStrings["terminationCondition="] = "simTime";

  filename_ = "None";
  
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

void Input::processBuffer(std::string& buffer)
{
  while(!fin_.eof())
    {
      std::getline(fin_, buffer);

      // Remove the whitespace

      boost::algorithm::erase_all(buffer, " ");
      
      // Check if it is a comment line

      if(buffer[0] == '#')
	continue;

      // Remove any comments on the line

      buffer = buffer.substr(0, buffer.find("#"));

      break;
    }
}

void Input::openInputFile(std::string filename)
{
  using namespace std;

  string buffer;

  if(fin_.is_open())
    fin_.close();

  filename_ = filename;
  fin_.open(filename.c_str());

  while(!fin_.eof())
    {
      processBuffer(buffer);
      
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
	  string errorMessage("Unrecognized input parameter type " + buffer + ".");

	  throw errorMessage.c_str();
	}
    }
}
