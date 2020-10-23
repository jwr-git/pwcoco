cmake_minimum_required(VERSION 3.10)
set (CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_C_STANDARD 11)

# set the project name
project(PWCOCO VERSION 1.0)

# Eigen library
if (DEFINED "$ENV{EIGEN3_INCLUDE_DIR}")
	set (EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")
	if (NOT EXISTS EIGEN3_INCLUDE_DIR)
    		message(FATAL_ERROR "Incorrect environment variable EIGEN3_INCLUDE_DIR for the location of Eigen3 package")
	endif()
endif()

# spdlog library
if (DEFINED "$ENV{SPDLOG_INCLUDE_DIR}")
	set (SPDLOG_INCLUDE_DIR "$ENV{SPDLOG_INCLUDE_DIR}")
	if (NOT EXISTS SPDLOG_INCLUDE_DIR)
		message (FATAL_ERROR "Incorrect environment variable SPDLOG_INCLUDE_DIR for the location of spdlog package")
	endif()
endif()

if (DEFINED "$ENV{EIGEN3_INCLUDE_DIR}" AND DEFINED "$ENV{SPDLOG_INCLUDE_DIR}")
	include_directories("${EIGEN3_INCLUDE_DIR}")
	include_directories("${SPDLOG_INCLUDE_DIR}")
else ()
	include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# add the executable
add_executable(PWCOCO src/options.cpp src/coloc.cpp src/conditional.cpp src/data.cpp src/dcdflib.cpp src/helper_funcs.cpp)