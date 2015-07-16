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

#include "Input.h"
#include "DomainInterface.h"
#include "Solver.h"
#include "Time.h"

class RunControl
{
public:
    RunControl();

    void initialize(Input& input);

    /** Decide whether or not to continue running.
     * @param timeStep the fixed timeStep.
     * @retval terminate if false, continue if true.
     */
    bool continueRun();
    bool writeToFile();
    void reset();
    double timeStep(){ return timeStep_; }
    double simTime() { return simTime_; }

    /** Create a directory. Platform independant.
     * @param directoryName The name of the directory to be created.
     */
    void createDirectory(std::string directoryName);

    //- Output messages

    void displayStartMessage();
    void displayUpdateMessage();
    void displayEndMessage();

    double residualNorm;

private:

    std::string terminationCondition_;
    int itrs_, maxIters_;
    double simTime_, timeStep_, maxSimTime_;

    int fileWriteInterval_;

    Time time_;
};

#endif
