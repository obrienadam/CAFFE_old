#ifndef TECPLOT_ASCII_H
#define TECPLOT_ASCII_H

#include <string>
#include <fstream>
#include <vector>

enum OutputFormat{POINT, BLOCK};

class TecplotAscii
{

    typedef std::vector< std::string > VariableList;

private:

    std::string filename_;

    std::ofstream fout_;
    
    OutputFormat outputFormat_;

    VariableList variableList_;

    int strandId_;

public:

    TecplotAscii(std::string filename = "Output.dat",
                 OutputFormat outputFormat = POINT,
                 int strandId = 1);

    ~TecplotAscii();

    void createHeader(VariableList variableList, int i, int j, int k);

};

#endif
