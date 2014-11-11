# CMakeLists.txt for CAFFE/src/RunControl
# Author: Adam O'Brien

set(HEADERS ArgsList.h
	    Version.h
	    Input.h
	    Output.h
	    RunControl.h
	    TecplotAscii.h)

set(SOURCE ArgsList.cc
   	   Version.cc
 	   Input.cc
           Output.cc
  	   RunControl.cc
	   TecplotAscii.cc)

add_library(RunControl ${HEADERS} ${SOURCE})
target_link_libraries(RunControl boost_program_options
				 boost_system
				 boost_date_time
				 boost_chrono)
