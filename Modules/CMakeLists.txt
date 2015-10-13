# CMakeLists.txt for CAFFE/Modules
# Author: Adam O'Brien

# Build the hexaMesh generator
add_executable(hexaMesh hexaMesh.cpp)
target_link_libraries(hexaMesh HexaFvmMesh)

# Build the diffusion solver
add_executable(caffeDiffusion caffeDiffusion.cpp)
target_link_libraries(caffeDiffusion RunControl
                                     Math
                                     Solvers
                                     HexaFvmMesh
                                     Schemes)

# Build the basic simple solver
add_executable(caffeSimple caffeSimple.cpp)
target_link_libraries(caffeSimple RunControl
                                  Math
                                  Solvers
                                  HexaFvmMesh
                                  Schemes)

# Build the basic piso solver
add_executable(caffePiso caffePiso.cpp)
target_link_libraries(caffePiso RunControl
                                Math
                                Solvers
                                HexaFvmMesh
                                Schemes)

# Build the immersed boundary piso solver
add_executable(caffeIbPiso caffeIbPiso.cpp)
target_link_libraries(caffeIbPiso RunControl
                                  Math
                                  Solvers
                                  HexaFvmMesh
                                  Schemes)

# Build the multi-phase piso solver
add_executable(caffeMultiphasePiso caffeMultiphasePiso.cpp)
target_link_libraries(caffeMultiphasePiso RunControl
                                          Math
                                          Solvers
                                          HexaFvmMesh
                                          Schemes)

# Build the coupled solver
add_executable(caffeCoupled caffeCoupled.cpp)
target_link_libraries(caffeCoupled RunControl
                                   Math
                                   Solvers
                                   HexaFvmMesh
                                   Schemes)
