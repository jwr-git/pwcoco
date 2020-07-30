cmake_minimum_required(VERSION 3.10)

# set the project name
project(PWCOCO VERSION 1.0)

# add the executable
add_executable(PWCOCO options.cpp coloc.cpp conditional.cpp data.cpp dcdflib.cpp helper_funcs.cpp)

SET(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")
IF(NOT EIGEN3_INCLUDE_DIR)
    MESSAGE(FATAL_ERROR "Specify environment variable EIGEN3_INCLUDE_DIR to the location of Eigen3 package")
ENDIF()
INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}")
