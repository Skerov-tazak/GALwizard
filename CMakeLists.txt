cmake_minimum_required(VERSION 3.5)
project(GALwizard
		VERSION 0.0.1
		DESCRIPTION ""
		LANGUAGES CXX)

#This Includes the Catch2 Framework
include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.8.1 # or a later release
)
FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_PATH}/extras)

#This Creates Targets for Complex Number Library
add_library(number STATIC src/complex/number.cpp)
target_include_directories(number PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/complex)
set_property(TARGET number PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(num_tests src/complex/number_test.cpp)
target_link_libraries(num_tests PRIVATE number Catch2::Catch2WithMain)

#This Creates Targets for The Linear Algebra Library
add_library(GALwizard SHARED src/functions.cpp src/matrix.cpp)
target_include_directories(GALwizard PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(GALwizard PUBLIC number)

add_executable(gal_tests test/test.cpp)
target_link_libraries(gal_tests PRIVATE GALwizard Catch2::Catch2WithMain)

add_executable(manual_test src/main.cpp)
target_link_libraries(manual_test PRIVATE GALwizard)


