/**
 * @file    InputStringProcessing.cc
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
 * This file contains the method implementations for class InputStringProcessing.
 */

#include <boost/algorithm/string.hpp>

#include "InputStringProcessing.h"

std::string InputStringProcessing::processBuffer(std::string &buffer, bool removeAllWhitespace)
{
    using namespace boost::algorithm;

    // Trim the leading/lagging whitespace

    trim(buffer);

    // Remove all whitespace from buffer if desired (default)

    if(removeAllWhitespace)
        erase_all(buffer, " ");

    // Check if it is a comment line, if so discard the input

    if(buffer[0] == '#')
    {
        buffer.clear();
    }

    // Remove any comments on the line

    buffer = buffer.substr(0, buffer.find("#"));
    return buffer;
}

double InputStringProcessing::getNextElement(std::string &buffer)
{
    using namespace boost::algorithm;

    double element;

    // Ensure the buffer is trimmed

    trim_left_if(buffer, is_any_of("( "));

    // Extract a double element from the string

    element = stod(buffer.substr(0, buffer.find_first_of(" )")));

    // Remove the extracted element from the string

    buffer = buffer.substr(buffer.find_first_of(" )"), buffer.back());
    return element;
}

std::vector<std::string> InputStringProcessing::partition(std::string buffer, std::string delimiter)
{
    using namespace std;
    using namespace boost::algorithm;

    vector<string> partitionedString;

    split(partitionedString, buffer, is_any_of(delimiter));

    for(int i = 0; i < partitionedString.size(); ++i)
    {
        if(partitionedString[i].empty())
        {
            partitionedString.erase(partitionedString.begin() + i);
            --i;
        }
    }

    return partitionedString;
}
