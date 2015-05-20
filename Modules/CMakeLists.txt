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

# Build the basic simple solver

add_executable(caffeSimple caffeSimple.cc)
target_link_libraries(caffeSimple RunControl
                                  Math
                                  Solvers
                                  HexaMesh
                                  Schemes)

# Build the immersed boundary simple solver

add_executable(caffeIbSimple caffeIbSimple.cc)
target_link_libraries(caffeIbSimple RunControl
                                    Math
                                    Solvers
                                    HexaMesh
                                    Schemes)

# Build the multiphase simple solver

add_executable(caffeMultiphaseSimple caffeMultiphaseSimple.cc)
target_link_libraries(caffeMultiphaseSimple RunControl
                                            Math
                                            Solvers
                                            HexaMesh
                                            Schemes)
