/**
 * @file    RunControl.cc
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

#include <sstream>

#include "RunControl.h"
#include "Output.h"

RunControl::RunControl()
    :
      itrs_(0),
      simTime_(0.),
      startRealTime_(boost::posix_time::microsec_clock::local_time()),
      maxElapsedRealTime_(boost::posix_time::hours(48)),
      terminationCondition_("iterations")
{

}

RunControl::RunControl(int argc, const char* argv[])
    :
      RunControl()
{
    argsList_.readArgs(argc, argv);
    input_.openInputFile(argsList_.inputFilename_);

    terminationCondition_ = input_.inputStrings["terminationCondition"];
    maxItrs_ = input_.inputInts["maxItrs"];
    maxSimTime_ = input_.inputDoubles["maxSimTime"];
}

bool RunControl::continueRun(double timeStep)
{
    using namespace boost::posix_time;

    ++itrs_;
    simTime_ += timeStep;
    elapsedRealTime_ = microsec_clock::local_time() - startRealTime_;

    if(terminationCondition_ == "iterations")
    {
        if(itrs_ >= maxItrs_)
            return false;
    }
    else if(terminationCondition_ == "simTime")
    {
        if(simTime_ >= maxSimTime_)
            return false;
    }
    else if (terminationCondition_ == "realTime")
    {
        if(elapsedRealTime_ >= maxElapsedRealTime_)
            return false;
    }
    else
    {
        Output::raiseException("RunControl", "continueRun", "invalid termination condition \"" + terminationCondition_ + "\" selected.");
    }

    return true;
}

void RunControl::displayStartMessage()
{
    using namespace std;

    ostringstream message;

    Output::printLine();

    message << "Beginning simulation. Terminating on condition: " << terminationCondition_ << "." << endl
            << "Iterations beginning on " << startRealTime_ << ".";

    Output::print(message.str());
}

void RunControl::displayUpdateMessage()
{
    using namespace std;
    using namespace boost::posix_time;

    ostringstream message;
    double completionPercentage;

    elapsedRealTime_ = microsec_clock::local_time() - startRealTime_;

    if(terminationCondition_ == "iterations")
    {
        completionPercentage = double(itrs_)/double(maxItrs_)*100.;
    }
    else if(terminationCondition_ == "simTime")
    {
        completionPercentage = simTime_/maxSimTime_*100.;
    }
    else if (terminationCondition_ == "realTime")
    {
        completionPercentage = elapsedRealTime_.total_microseconds()/maxElapsedRealTime_.total_microseconds()*100.;
    }

    message << "Simulation completion: " << completionPercentage << "%" << endl
            << "Iterations completed: " << itrs_ << endl
            << "Simulation time: " << simTime_ << endl
            << "Elapsed time: " << elapsedRealTime_;

    Output::print(message.str());
}

void RunControl::displayEndMessage()
{
    using namespace std;
    using namespace boost::posix_time;

    elapsedRealTime_ = microsec_clock::local_time() - startRealTime_;

    ostringstream message;

    Output::printLine();

    message << "Iterations completed on " << startRealTime_ + elapsedRealTime_ << endl
            << "Iterations completed: " << itrs_ << endl
            << "Simulation time: " << simTime_ << endl
            << "Elapsed time: " << elapsedRealTime_;

    Output::print(message.str());
    Output::printLine();
}

void RunControl::initializeCase(Solver& solver, DomainInterface& domain)
{
    domain.initialize(input_);
    solver.initialize(input_);
}
