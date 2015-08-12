# CMakeLists.txt for CAFFE/src/Solvers
# Author: Adam O'Brien

set(HEADERS Solver.h
            Diffusion.h
            Simple.h
            Piso.h
            IbPiso.h
            Coupled.h)

set(SOURCE Solver.cpp
           Diffusion.cpp
           Simple.cpp
           Piso.cpp
           IbPiso.cpp
           Coupled.cpp)

add_library(Solvers ${HEADERS} ${SOURCE})
target_link_libraries(Solvers Field)
