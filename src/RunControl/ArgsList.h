#ifndef ARGS_LIST_H
#define ARGS_LIST_H

#include <boost/program_options.hpp>

class ArgsList
{

private:

    boost::program_options::options_description optsDescription_;
    boost::program_options::variables_map varsMap_;

public:

    ArgsList();
    ArgsList(int argc, const char* argv[]);

    void readArgs(int argc, const char* argv[]);

};

#endif
