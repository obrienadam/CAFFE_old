# CMakeLists.txt for CAFFE/Modules
# Author: Adam O'Brien

# Build the hexaMesh generator

add_executable(hexaMesh hexaMesh.cc)
target_link_libraries(hexaMesh HexaFvmMesh)

# Build the diffusion solver
add_executable(caffeDiffusion caffeDiffusion.cc)
target_link_libraries(caffeDiffusion RunControl
                                     Math
                                     Solvers
                                     HexaFvmMesh
                                     Schemes)

# Build the basic simple solver
add_executable(caffeSimple caffeSimple.cc)
target_link_libraries(caffeSimple RunControl
                                  Math
                                  Solvers
                                  HexaFvmMesh
                                  Schemes)

# Build the basic piso solver
add_executable(caffePiso caffePiso.cc)
target_link_libraries(caffePiso RunControl
                                Math
                                Solvers
                                HexaFvmMesh
                                Schemes)

# Build the immersed boundary piso solver
add_executable(caffeIbPiso caffeIbPiso.cc)
target_link_libraries(caffeIbPiso RunControl
                                  Math
                                  Solvers
                                  HexaFvmMesh
                                  Schemes)

# Build the multiphase simple solver
add_executable(caffeMultiphaseSimple caffeMultiphaseSimple.cc)
target_link_libraries(caffeMultiphaseSimple RunControl
                                            Math
                                            Solvers
                                            HexaFvmMesh
                                            Schemes)
