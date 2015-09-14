/**
 * @file    Output.cpp
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
 * This file contains the implementations for class Output.
 */

#include <iostream>
#include <iomanip>

#include "Output.h"
#include "Parallel.h"

void Output::print(const std::string &message)
{
    if(Parallel::isMainProcessor())
        std::cout << message << std::endl;
}

void Output::print(const std::string &className, const std::string &message)
{
    if(Parallel::isMainProcessor())
        std::cout << className + ": " << message << std::endl;
}

void Output::print(const Array3D<int> &array3D)
{
    if(Parallel::isMainProcessor())
    {
        for(int k = 0; k < array3D.sizeK(); ++k)
        {
            for(int j = 0; j < array3D.sizeJ(); ++j)
            {
                for(int i = 0; i < array3D.sizeI(); ++i)
                {
                    std::cout << array3D(i, j, k) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
}

void Output::raiseException(std::string className, std::string methodName, std::string problemDescription)
{
    if(Parallel::isMainProcessor())
        throw ("in \"" + className + "::" + methodName + "\", " + problemDescription).c_str();
}

void Output::issueWarning(std::string className, std::string methodName, std::string warningDescription)
{
    if(Parallel::isMainProcessor())
        std::cout << "WARNING! in \"" << className << "::" << methodName << "\", " << warningDescription << std::endl;
}

void Output::pause()
{
    if(Parallel::isMainProcessor())
    {
        int i;
        print("Press Enter to continue...");
        std::cin >> i;
    }
    Parallel::barrier();
}

void Output::printCaffeHeader()
{
    using namespace std;

    if(Parallel::isMainProcessor())
    {
        printLine();

        cout << "|  || ___\n"
             << "|  ||/ _  \\ \\\n"				\
             << "|  ||\\__/ | | |\n"
             << "|  | \\___/  | | CAFFE\n"
             << "|   \\______/  /\n"
             << " \\___________/\n"
             << endl
             << "  Computational Algorithm Framework for Fluid Equations (CAFFE)\n"
             << endl
             << "                      Author: Adam O'Brien\n"
             << "	            E-mail: roni511@gmail.com\n"
             << endl;

        printLine();
    }
}

void Output::printLine()
{
    if(Parallel::isMainProcessor())
        std::cout << "------------------------------------------------------------------" << std::endl;
}
