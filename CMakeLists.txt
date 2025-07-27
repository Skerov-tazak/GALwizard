cmake_minimum_required(VERSION 3.5)
project(GALwizard
		VERSION 0.0.1
		DESCRIPTION ""
		LANGUAGES CXX)
add_library(libGALwizard SHARED src/functions.cpp src/matrix.cpp)
target_include_directories(libGALwizard.so PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)


add_executable(GALtest 
				test/testing.cpp)
target_include_directories(GALtest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(GALtest libGALwizard)
