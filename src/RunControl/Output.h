#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

class Output
{

public:

    //- These functions are for arbitrary screen output

    static void printToScreen(const std::string& message);
    static void printToScreen(const std::vector<std::string>& vector);
    static void printToScreen(const std::ostream& message);

    //- These functions are for normalizing exception handling

    static void raiseException(std::string className, std::string methodName, std::string problemDescription);
    static void raiseException(std::string className, std::string methodName, std::ostringstream& problemDescription);

    static void displayCaffeHeader();
    static void printLine();

};

#endif
