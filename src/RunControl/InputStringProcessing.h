/**
 * @file    InputStringProcessing.h
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
 * This file contains the interface for class InputStringProcessing,
 * which is simply a collection of static methods for manipulating
 * strings.
 */

#ifndef INPUTSTRINGPROCESSING_H
#define INPUTSTRINGPROCESSING_H

#include <string>
#include <vector>

class InputStringProcessing
{
private:

public:

    /**
     * @brief Process a buffer containing a string that has been read in from a line of a file.
     * @param removeAllWhitespace A boolean that when true, removes all whitespace from the string and when false, only removes leading and trailing whitespace.
     * @return A string that is partititioned into sections.
     */
    static std::string processBuffer(std::string& buffer, bool removeAllWhitespace = true);

    /**
     * @brief Read in a string, and return the next floating point number found.
     * @return A floating point number found in the buffer string.
     */
    static double getNextElement(std::string &buffer);

    /**
     * @brief Partition a string into sections, which are delimited on a specified delimiter.
     * @param delimiter The string to be used for delimination.
     * @return A vector containing the partitioned string.
     */
    static std::vector<std::string> partition(std::string buffer, std::string delimiter);
};

#endif
