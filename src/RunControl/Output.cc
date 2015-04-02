/**
 * @file    Output.c
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

#include "Output.h"

void Output::print(const std::string& message)
{
    std::cout << std::endl << message << std::endl;
}

void Output::print(const std::ostream& message)
{
    std::cout << std::endl << message << std::endl;
}

void Output::print(const std::vector<std::string>& vector)
{
    for(uint i = 0; i < vector.size(); ++i)
    {
        std::cout << "Element " << i << ": " << vector[i] << std::endl;
    }
}

void Output::print(std::string className, std::string message)
{
    std::cout << className + ": " << message << std::endl;
}

void Output::raiseException(std::string className, std::string methodName, std::string problemDescription)
{
    throw ("in \"" + className + "::" + methodName + "\", " + problemDescription).c_str();
}

void Output::raiseException(std::string className, std::string methodName, std::ostringstream &problemDescription)
{
    throw ("in \"" + className + "::" + methodName + "\", " + problemDescription.str()).c_str();
}

void Output::printCaffeHeader()
{
    using namespace std;


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

void Output::printLine()
{
    std::cout << "------------------------------------------------------------------" << std::endl;
}
