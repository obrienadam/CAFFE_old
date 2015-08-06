# CMakeLists.txt for CAFFE/src/Domains/Field
# Author: Adam O'Brien

set(HEADERS Field.h
            PrimitiveBoundaryCondition.h
            FlowBoundaryConditions.h
            SimpleBoundaryCondition.h
            InitialConditions.h)

set(SOURCE FieldI.h
           PrimitiveBoundaryConditionI.h
           PrimitiveBoundaryCondition.cc
           FlowBoundaryConditions.cc
           SimpleBoundaryCondition.cc
           InitialConditionsI.h
           InitialConditions.cc)

add_library(Field ${HEADERS} ${SOURCE})
target_link_libraries(Field)
