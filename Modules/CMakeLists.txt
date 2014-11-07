# CMakeLists.txt for CAFFE/Modules
# Author: Adam O'Brien

set(CAD_LIBS RunControl Math State)

add_executable(caffeAdvectionDiffusion caffeAdvectionDiffusion.cc)
target_link_libraries(caffeAdvectionDiffusion ${CAD_LIBS})
