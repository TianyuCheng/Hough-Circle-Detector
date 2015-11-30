# project name
PROJECT(HOUGH_CIRCLE_DETECTOR)

# 2.4.3 is actual minimum, but you have to add additional hacks
# uncomment this if you need 2.4.3 compatibility
# if(COMMAND cmake_policy)
#	 cmake_policy(SET CMP0003 NEW)
# endif(COMMAND cmake_policy)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# find and setup Qt4 for this project
FIND_PACKAGE(Qt4 REQUIRED)

# user-defined flags
# ADD_DEFINITIONS(-DTIMER=1)

# add cxx flags
SET(CMAKE_CXX_FLAGS "-O3 --std=c++11 ${CMAKE_CXX_FLAGS}")
# SET(CMAKE_CXX_FLAGS "-ggdb3 ${CMAKE_CXX_FLAGS}")
# SET(CMAKE_CXX_FLAGS "-pg ${CMAKE_CXX_FLAGS}")

# tell cmake to process CMakeLists.txt in these subdirs
SUBDIRS(original)
SUBDIRS(parallel)
