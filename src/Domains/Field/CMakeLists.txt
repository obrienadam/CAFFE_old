# CMakeLists.txt for CAFFE/src/Domains/Field
# Author: Adam O'Brien

set(HEADERS Field.h
            PrimitiveBoundaryCondition.h
            SimpleBoundaryCondition.h
            InitialConditions.h)

set(SOURCE FieldI.h
           PrimitiveBoundaryConditionI.h
           PrimitiveBoundaryCondition.cc
           SimpleBoundaryCondition.cc
           InitialConditionsI.h)

add_library(Field ${HEADERS} ${SOURCE})
target_link_libraries(Field)
