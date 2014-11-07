#ifndef SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H

#include <string>
#include <sstream>
#include <vector>

#include "SmartPointer.h"
#include "DomainIncludes.h"
#include "SchemeIncludes.h"

#include "Input.h"

template <class DOMAIN_TYPE, class STATE_TYPE>
class SolverInterface
{

protected:

    //- Typedefs for the longer names

    typedef std::vector< std::vector<STATE_TYPE> > Array2D;
    typedef typename DOMAIN_TYPE::iterator Iterator;

    //- The base constructor (called by all derived classes)

    SolverInterface(std::string solverName = "Uninitialized Solver",
                    std::string timeUnits = "seconds")
        :
          solverName_(solverName),
          timeUnits_(timeUnits)
    {

    }

    std::string solverName_;

    //- Solver elements common to all sovlers

    int nSteps_, nElements_;
    double timeStep_;
    std::string timeUnits_;

    //- The domain and a data structure for storing the time derivatives of the domain

    DOMAIN_TYPE domain_;
    Array2D timeDerivatives_;

    //- Iterators for the domain

    Iterator itr_;
    Iterator begin_;
    Iterator end_;

    //- Helper function for initialization

    void initializeNumOfSteps(int nSteps, int nElements);

public:

    virtual void initialize(Input& input);

    virtual void solve() = 0;

    virtual double timeStep();

    std::string timeUnits();

};

#include "SolverInterfaceI.h"

#endif
