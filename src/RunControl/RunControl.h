/**
 * @file    RunControl.h
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
 * This file contains the interface for class RunControl, which is
 * responsible for setting up and controlling the flow of simulations.
 */

#ifndef RUN_CONTROL_H
#define RUN_CONTROL_H

#include <string>

#include <boost/date_time.hpp>
#include <boost/chrono.hpp>

#include "ArgsList.h"
#include "Input.h"

#include "DomainInterface.h"
#include "Solver.h"

class RunControl
{
    typedef boost::posix_time::ptime RealTime;
    typedef boost::posix_time::time_duration RealTimeDuration;
    typedef boost::chrono::process_cpu_clock::duration CpuTimeDuration;

private:

    ArgsList argsList_;
    Input input_;

    //- Simulation control

    std::string terminationCondition_;
    int itrs_, maxItrs_;
    double simTime_, maxSimTime_;

    //- Time related objects

    RealTime startTime_;
    RealTimeDuration elapsedTime_, maxElapsedTime_;
    CpuTimeDuration cpuTime_, maxCpuTime_;

    //- Private constructor only used to initialize defaults

    RunControl();

public:

    /** Constructor that accepts command line arguments.
     * @param argc number of command line arguments. Expects two.
     * @param argv command line arguments, should be a filename.
     */
    RunControl(int argc, const char* argv[]);

    /** Decide whether or not to continue running.
     * @param timeStep the fixed timeStep.
     * @retval terminate if false, continue if true.
     */
    bool continueRun(double timeStep = 0.);

    //- Output messages

    void displayStartMessage();
    void displayUpdateMessage();
    void displayEndMessage();

    /** Initializes a Solver and Domain by passing it the input.
     * @param solver the solver.
     * @param domain the computational domain.
     */
    void initializeCase(Solver& solver, DomainInterface& domain);
};

#endif
