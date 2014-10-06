#ifndef ARGS_LIST_H
#define ARGS_LIST_H

#include <string>
#include "SmartPointer.h"

class ArgsList
{

private:

    SmartPointer<std::string> options_;
    SmartPointer<std::string> args_;

public:

    ArgsList();
    ArgsList(int argc, const char* argv[]);

    void readArgs(int argc, const char* argv[]);

};

#endif
