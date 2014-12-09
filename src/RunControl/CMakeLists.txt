# CMakeLists.txt for CAFFE/src/RunControl
# Author: Adam O'Brien

set(HEADERS ArgsList.h
	    Version.h
            InputStringProcessing.h
	    Input.h
	    Output.h
            RunControl.h)

set(SOURCE ArgsList.cc
   	   Version.cc
           InputStringProcessing.cc
 	   Input.cc
           Output.cc
           RunControl.cc)

add_library(RunControl ${HEADERS} ${SOURCE})
target_link_libraries(RunControl ${Boost_LIBRARIES})
