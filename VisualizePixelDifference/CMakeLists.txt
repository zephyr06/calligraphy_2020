cmake_minimum_required(VERSION 3.10)

# set the project name
project (PixelDifference)

# add the executables
add_executable(PixelDifference
    src/PixelDifference.cpp)
add_executable(PixelDifferenceTest
    test/TestPixelDifferenceMain.cpp
    test/TestPixelDifferenceTests.cpp)

# set include paths for Eigen
find_package(Eigen3 QUIET)
include_directories(${EIGEN3_INCLUDE_DIR})
set (USE_EIGEN=/usr/include/eigen3)

# set include paths for OpenCV
find_package( OpenCV REQUIRED )
include_directories(${OpenCV_INCLUDE_DIRS})
target_link_libraries( PixelDifferenceTest ${OpenCV_LIBS} )
target_link_libraries( PixelDifference ${OpenCV_LIBS} )

# specify the C++ standard
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED on)
set (CMAKE_CXX_FLAGS 
    "${CMAKE_CXX_FLAGS} -std=c++17")
set (VERBOSE on)

# set version numbers
set (PixelDifference_VERSION_MAJOR 1)
set (PixelDifference_VERSION_MINOR 2)

# configure a header file to pass some of the CMake settings to source code
configure_file (
    "${PROJECT_SOURCE_DIR}/src/PixelDifferenceMethods.h.in"
    "${PROJECT_SOURCE_DIR}/src/PixelDifferenceMethods.h"
)

target_include_directories(PixelDifference
    PUBLIC "${PROJECT_BINARY_DIR}")
