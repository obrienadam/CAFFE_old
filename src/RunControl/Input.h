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
 * The input class contains all input fields required for setting up
 * and running cases. Contains necessary methods for reading input
 * files.
 */

#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <fstream>
#include <map>

class Input
{
private:

    std::string filename_;
    std::ifstream fin_;                                 ///< The input file stream

public:

    std::map<std::string, int> inputInts;
    std::map<std::string, double> inputDoubles;
    std::map<std::string, std::string> inputStrings;

    /** Default constructor. Sets all default input values.
     */
    Input();

    /** Constructor that opens an input file.
     * @param filename the name of the input file.
     */
    Input(std::string filename);
    ~Input();

    /** Open an input file.
     * @param filename the name of the input file.
     */
    void openInputFile(std::string filename);

    /** Print all of the input data to the console.
     */
    void print();
};

#endif
