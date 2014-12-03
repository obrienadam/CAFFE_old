#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
#include <string>
#include <vector>

class Output
{

public:

    //- This static method can be replaced later to allow for parallel applications, GUIs etc.

    static void printToScreen(const std::string& message);

    static void printToScreen(const std::vector<std::string>& vector);

    static void printToScreen(const std::ostream& message);

    static void displayCaffeHeader();

    static void printLine();

};

#endif
