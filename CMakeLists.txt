cmake_minimum_required(VERSION 3.10)
set (CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_C_STANDARD 11)

# set the project name
project(pwcoco VERSION 1.0 LANGUAGES CXX)

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
add_executable(pwcoco src/options.cpp src/coloc.cpp src/conditional.cpp src/data.cpp src/dcdflib.cpp src/helper_funcs.cpp)
target_link_libraries(pwcoco PRIVATE OpenMP::OpenMP_CXX)
