#infndef TECPLOT_ASCII_H
#define TECPLOT_ASCII_H

#include <string>
#include <fstream>

class TecplotAscii
{

 private:

  std::ofstream fout_;

 public:

  TecplotAscii();
  TecplotAscii(std::string filename, int nVariables, std::string variableNames[], int strandID = 1);

  
};

#endif
