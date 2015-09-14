/**
 * @file    Output.h
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
 * This file contains the interface for class Output, which is used
 * for printing to the screen and raising exceptions.
 */

#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "Array3D.h"

class Output
{
public:

    static void print(const std::string &message);
    static void print(const std::string &className, const std::string &message);
    static void print(const Array3D<int> &array3D);

    /**
     * @brief Used to raise an exception, while also making sure relevant information about the method raising the exception.
     * @param className The name of the class raising the exception.
     * @param methodName The method of the class that is raising the exception.
     * @param problemDescription A description of the problem that occurred.
     */
    static void raiseException(std::string className, std::string methodName, std::string problemDescription);
    static void issueWarning(std::string className, std::string methodName, std::string warningDescription);

    static void pause();

    static void printCaffeHeader();
    static void printLine();
};

#endif
