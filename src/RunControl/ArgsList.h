#ifndef ARGS_LIST_H
#define ARGS_LIST_H

#include <string>

#include <boost/program_options.hpp>

class ArgsList
{

private:

    typedef boost::program_options::options_description OptsDescription;
    typedef boost::program_options::variables_map VarsMap;

private:

    OptsDescription optsDescription_;
    VarsMap varsMap_;

    std::string inputFilename_;

    bool validOptionsSelected_;

public:

    ArgsList();
    ArgsList(int argc, const char* argv[]);

    void readArgs(int argc, const char* argv[]);

    friend class RunControl;
    friend class HexaMeshGen;

};

#endif
