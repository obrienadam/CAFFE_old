# CMakeLists.txt for CAFFE/src/RunControl
# Author: Adam O'Brien

set(SOURCE_FILES ArgsList.cc
		 Input.cc
		 Output.cc
		 RunControl.cc
		 TecplotAscii.cc)

add_library(RunControl ${SOURCE_FILES})
target_link_libraries(RunControl boost_program_options
				 boost_system
				 boost_date_time
				 boost_chrono)
