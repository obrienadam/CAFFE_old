/**
 * @file    Input.cc
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
 * This file contains all of the implementations for the methods of
 * class Input.
 */

#include <boost/property_tree/info_parser.hpp>

#include "Input.h"

Input::Input()
{
    using namespace boost::property_tree;

    read_info("case/case.info", caseParameters);
    read_info("case/initialConditions.info", initialConditions);
}

Input::Input(const std::string &caseFilename, const std::string &initialConditionsFilename)
{
    using namespace boost::property_tree;

    read_info(caseFilename, caseParameters);
    read_info(initialConditionsFilename, initialConditions);
}
