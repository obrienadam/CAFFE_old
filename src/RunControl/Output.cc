#include <iostream>

#include "Output.h"

void Output::printToScreen(const std::string& message)
{
    std::cout << std::endl << message << std::endl;
}

void Output::printToScreen(const std::ostream& message)
{
    std::cout << std::endl << message << std::endl;
}

void Output::printToScreen(const std::vector<std::string>& vector)
{
    for(uint i = 0; i < vector.size(); ++i)
    {
        std::cout << "Element " << i << ": " << vector[i] << std::endl;
    }
}

void Output::printToScreen(std::string className, std::string message)
{
    std::cout << className + ": " << message << std::endl;
}

void Output::raiseException(std::string className, std::string methodName, std::string problemDescription)
{
    throw ("in \"" + className + "::" + methodName + "\", " + problemDescription).c_str();
}

void Output::raiseException(std::string className, std::string methodName, std::ostringstream &problemDescription)
{
    throw ("in \"" + className + "::" + methodName + "\", " + problemDescription.str()).c_str();
}

void Output::displayCaffeHeader()
{
    using namespace std;


    printLine();

    cout << "|  || ___\n"
         << "|  ||/ _  \\ \\\n"				\
         << "|  ||\\__/ | | |\n"
         << "|  | \\___/  | | CAFFE\n"
         << "|   \\______/  /\n"
         << " \\___________/\n"
         << endl
         << "  Computational Algorithm Framework for Fluid Equations (CAFFE)\n"
         << endl
         << "                      Author: Adam O'Brien\n"
         << "	            E-mail: roni511@gmail.com\n"
         << endl;

    printLine();
}

void Output::printLine()
{
    std::cout << "------------------------------------------------------------------" << std::endl;
}
