#ifndef SIMPLE_H
#define SIMPLE_H

#include "FvScheme.h"

class Simple : public FvScheme
{
public:

    int nConservedVariables();

    void discretize(std::vector<double>& timeDerivatives_);
    void copySolution(std::vector<double>& original);
    void updateSolution(std::vector<double>& timeDerivatives_, int method);
};

#endif
