# CMakeLists.txt for CAFFE/Modules
# Author: Adam O'Brien

# Build the hexaMesh generator

add_executable(hexaMesh hexaMesh.cc)
target_link_libraries(hexaMesh HexaMesh)

# Build the advection diffusion solver

add_executable(caffeAdvectionDiffusion caffeAdvectionDiffusion.cc)
target_link_libraries(caffeAdvectionDiffusion RunControl
                                              Math
                                              Solvers
                                              HexaMesh
                                              Schemes)

add_executable(caffeSimple caffeSimple.cc)
target_link_libraries(caffeSimple RunControl
                                  Math
                                  Solvers
                                  HexaMesh
                                  Schemes)
