# CMakeLists.txt for CAFFE/src/Domains/Field
# Author: Adam O'Brien

set(HEADERS Field.h
            PrimitiveBoundaryCondition.h
            FlowBoundaryConditions.h
            SimpleBoundaryCondition.h
            InitialConditions.h)

set(SOURCE Field.tpp
           PrimitiveBoundaryCondition.tpp
           PrimitiveBoundaryCondition.cpp
           FlowBoundaryConditions.cpp
           SimpleBoundaryCondition.cpp
           InitialConditions.tpp
           InitialConditions.cpp)

add_library(Field ${HEADERS} ${SOURCE})
target_link_libraries(Field)
