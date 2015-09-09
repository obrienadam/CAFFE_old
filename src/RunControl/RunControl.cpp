/**
 * @file    RunControl.cpp
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
 * This file contains the implementations for class RunControl.
 */

#include <boost/filesystem.hpp>

#include "RunControl.h"
#include "Output.h"

RunControl::RunControl()
    :
      itrs_(0),
      simTime_(0.),
      terminationCondition_("iterations")
{

}

void RunControl::initialize(const Input &input)
{
    using namespace std;

    terminationCondition_ = input.caseParameters.get<string>("RunControl.terminationCondition");
    maxIters_ = input.caseParameters.get<int>("RunControl.maxNumberOfIterations");
    timeStep_ = input.caseParameters.get<double>("RunControl.timeStep");
    maxSimTime_ = input.caseParameters.get<double>("RunControl.maxSimTime");
    fileWriteInterval_ = input.caseParameters.get<int>("RunControl.fileWriteInterval");

    //- Create a directory for the solution output
    createDirectory("solution");
}

bool RunControl::continueRun()
{

    ++itrs_;
    simTime_ += timeStep_;

    if(terminationCondition_ == "iterations")
    {
        if(itrs_ >= maxIters_)
            return false;
    }
    else if(terminationCondition_ == "simTime")
    {
        if(simTime_ >= maxSimTime_)
            return false;
    }
    else if (terminationCondition_ == "realTime")
    {
        Output::issueWarning("RunControl", "continueRun", "termination condition for real time not yet implemented.");
    }
    else
    {
        Output::raiseException("RunControl", "continueRun", "invalid termination condition \"" + terminationCondition_ + "\" selected.");
    }

    return true;
}

bool RunControl::writeToFile()
{
    if(itrs_%fileWriteInterval_ == 0)
        return true;

    return false;
}

void RunControl::reset()
{
    Output::raiseException("RunControl", "reset", "Method not yet implemented.");
}

void RunControl::createDirectory(std::string directoryName)
{
    boost::filesystem::path dir(directoryName);

    if(!boost::filesystem::exists(dir))
    {
        if(!boost::filesystem::create_directory(dir))
            Output::raiseException("RunControl", "createDirectory", "Creation of directory \"" + directoryName + "\" failed.");
    }
}

void RunControl::displayStartMessage()
{
    using namespace std;

    Output::print("RunControl", "Beginning simulation. Terminating on condition: " + terminationCondition_);
    Output::print("RunControl", "Iterations beginning on " + time_.currentTime() + ".");
}

void RunControl::displayUpdateMessage()
{
    using namespace std;
    using namespace boost::posix_time;

    ostringstream message;
    double completionPercentage;

    if(terminationCondition_ == "iterations")
    {
        completionPercentage = double(itrs_)/double(maxIters_)*100.;
    }
    else if(terminationCondition_ == "simTime")
    {
        completionPercentage = simTime_/maxSimTime_*100.;
    }
    else if (terminationCondition_ == "realTime")
    {
        //completionPercentage = elapsedRealTime_.total_microseconds()/maxElapsedRealTime_.total_microseconds()*100.;
    }

    message << "Simulation completion (%) |      " << completionPercentage << endl
            << "Iterations completed      |      " << itrs_ << endl
            << "Simulation time (sec)     |      " << simTime_ << endl
            << "Elapsed time (hh:mm:ss)   |      " << time_.elapsedTime() << endl
            // << "CPU time (microseconds)   |      " << time_.elapsedCpuTimeMicroseconds() << endl
            << "Residual norm             |      " << residualNorm;

    Output::print("RunControl", message.str());
    Output::printLine();
}

void RunControl::displayEndMessage()
{
    using namespace std;

    ostringstream message;

    Output::printLine();

    message << "Iterations completed on " << time_.currentTime() << endl
            << "Iterations completed: " << itrs_ << endl
            << "Simulation time: " << simTime_ << endl
            << "Elapsed time: " << time_.elapsedTime();

    Output::print("RunControl", message.str());
    Output::printLine();
}
