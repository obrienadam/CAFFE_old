#include <iostream>

#include "Output.h"

void Output::printToScreen(const std::string& message)
{

    std::cout << std::endl << message << std::endl;

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
