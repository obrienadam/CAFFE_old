# CMakeLists.txt for CAFFE/src/RunControl
# Author: Adam O'Brien

set(HEADERS ArgsList.h
	    Input.h
	    Output.h
            RunControl.h
            Time.h)

set(SOURCE ArgsList.cc
 	   Input.cc
           Output.cc
           RunControl.cc
           Time.cc)

add_library(RunControl ${HEADERS} ${SOURCE})
target_link_libraries(RunControl Parallel ${Boost_LIBRARIES})
