# CMakeLists.txt for CAFFE/src/Solvers
# Author: Adam O'Brien

set(HEADERS Solver.h
	    Euler.h)

set(SOURCE Solver.cc
	   Euler.cc)

add_library(Solvers ${HEADERS} ${SOURCE})
