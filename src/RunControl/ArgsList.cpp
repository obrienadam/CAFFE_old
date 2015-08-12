/**
 * @file    ArgsList.cpp
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
 * This file contains the method implementations for class ArgsList.
 */

#include <iostream>
#include <string>

#include "ArgsList.h"

ArgsList::ArgsList()
    :
      optsDescription_("Supported options")
{
    using namespace boost::program_options;

    optsDescription_.add_options()
            ("help", "| shows this help message")
            ("version", "| show version info")
            ("file", value<std::string>(), "| load input file");
}

ArgsList::ArgsList(int argc, const char* argv[])
    :
      ArgsList()
{
    readArgs(argc, argv);
}

void ArgsList::readArgs(int argc, const char* argv[])
{
    using namespace std;
    using namespace boost::program_options;

    store(parse_command_line(argc, argv, optsDescription_), varsMap_);
    notify(varsMap_);
}
