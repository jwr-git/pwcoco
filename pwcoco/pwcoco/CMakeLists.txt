cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 11)

# set the project name
project(PWCOCO VERSION 1.0)

# add the executable
add_executable(PWCOCO options.cpp coloc.cpp conditional.cpp data.cpp dcdflib.cpp helper_funcs.cpp)

set (EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")
if (NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Specify environment variable EIGEN3_INCLUDE_DIR to the location of Eigen3 package")
endif()
include_directories("${EIGEN3_INCLUDE_DIR}")
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()