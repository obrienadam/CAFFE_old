/**
 * @file    ArgsList.h
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
 * This file contains the interface for class Input, which is used
 * for reading and processing command line arguments, and passing
 * those arguments to the Input class.
 */

#ifndef ARGS_LIST_H
#define ARGS_LIST_H

#include <string>

#include <boost/program_options.hpp>

class ArgsList
{

public:

    ArgsList();
    ArgsList(int argc, const char* argv[]);

    void readArgs(int argc, const char* argv[]);

private:

    typedef boost::program_options::options_description OptsDescription;
    typedef boost::program_options::variables_map VarsMap;

    OptsDescription optsDescription_;
    VarsMap varsMap_;
};

#endif
