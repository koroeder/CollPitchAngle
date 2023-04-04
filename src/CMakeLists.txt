cmake_minimum_required(VERSION 3.0)

if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")
  message( FATAL_ERROR "CMake should not be built in the source directory" )
endif (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")

project(PITCHANGLE)

set(PROJECT_VERSION "1.0")
set(PROJECT_CONTACT "konstantin.roeder@kcl.ac.uk")

