# CMakeLists.txt for CAFFE/src/Solvers
# Author: Adam O'Brien

set(HEADERS Solver.h
            Diffusion.h
            Simple.h
            Piso.h
            IbPiso.h)

set(SOURCE Solver.cc
           Diffusion.cc
           Simple.cc
           Piso.cc
           IbPiso.cc)

add_library(Solvers ${HEADERS} ${SOURCE})
target_link_libraries(Solvers ${Boost_LIBRARIES})
