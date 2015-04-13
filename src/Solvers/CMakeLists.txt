# CMakeLists.txt for CAFFE/src/Solvers
# Author: Adam O'Brien

set(HEADERS SolverIncludes.h
            Solver.h
            Euler.h
            PredictorCorrector.h)

set(SOURCE Solver.cc
           Euler.cc
           PredictorCorrector.cc)

add_library(Solvers ${HEADERS} ${SOURCE})
