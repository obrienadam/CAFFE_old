/**
 * @file    Input.h
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
 * This file contains the interface for class Input, which contains
 * fields for setting up and running simulations.
 */

#ifndef INPUT_H
#define INPUT_H

#include <string>

#include <boost/property_tree/ptree.hpp>

class Input
{
public:

    Input();
    Input(const std::string &caseFilename, const std::string &initialConditionsFilename);

    boost::property_tree::ptree caseParameters;
    boost::property_tree::ptree initialConditions;
};

#endif
