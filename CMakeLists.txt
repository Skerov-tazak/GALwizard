cmake_minimum_required(VERSION 3.5)
project(GALwizard
		VERSION 0.0.1
		DESCRIPTION ""
		LANGUAGES CXX)

add_library(GALwizard SHARED src/functions.cpp src/matrix.cpp)
target_include_directories(GALwizard PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.8.1 # or a later release
)
FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_PATH}/extras)

add_executable(tests test/test.cpp)
target_link_libraries(tests PRIVATE GALwizard Catch2::Catch2WithMain)
