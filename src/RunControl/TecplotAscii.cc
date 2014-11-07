#include "TecplotAscii.h"

TecplotAscii::TecplotAscii(std::string filename,
                           OutputFormat outputFormat,
                           int strandId)
    :
      filename_(filename),
      outputFormat_(outputFormat),
      strandId_(strandId)
{

    fout_.open(filename_.c_str());

}

TecplotAscii::~TecplotAscii()
{

    fout_.close();

}

void TecplotAscii::createHeader(VariableList variableList, int i, int j, int k)
{

    int itr;

    variableList_ = variableList;

    fout_ << "TITLE = \"" << filename_ << "\"\n"
          << "VARIABLES = ";

    for(itr = 0; itr < variableList_.size(); ++itr)
    {

        fout_ << "\"" << variableList_[i] << "\" ";

    }

    fout_ << "\n"
          << "ZONE i=" << i << ", j=" << j << ", k=" << k << ", f=";

    if(outputFormat_ == POINT)
    {

        fout_ << "POINT\n";

    }
    else if (outputFormat_ == BLOCK)
    {

        fout_ << "BLOCK\n";

    }
    else
    {

        throw "Unrecognized data format type in TecplotAscii.";

    }

}


