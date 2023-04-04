cmake_minimum_required(VERSION 3.0)

if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")
  message( FATAL_ERROR "CMake should not be built in the source directory" )
endif (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")

project(PITCHANGLE)

set(PROJECT_VERSION "1.0")
set(PROJECT_CONTACT "konstantin.roeder@kcl.ac.uk")

get_filename_component(GIT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)
message("Root: ${CMAKE_CURRENT_SOURCE_DIR}") 

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

enable_language(Fortran)
set (CMAKE_Fortran_FLAGS "-ffixed-line-length-132 -ffree-line-length-0" CACHE TYPE STRING FORCE)
set (CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE TYPE STRING FORCE)
set (CMAKE_Fortran_FLAGS_DEBUG "-g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none" CACHE TYPE STRING FORCE)
set (FORTRAN_FREEFORM_FLAG "-ffree-form" CACHE TYPE STRING)

include_directories(${CMAKE_SOURCE_DIR})
# Glob all the sources
file(GLOB PITCH_LIB_SOURCES ${GIT_ROOT}/src/*.f90)
message("All sources: ${PITCH_LIB_SOURCES}") 
message("Root: ${SOURCE_DIR}") 
add_library(pitchlib ${PITCH_LIB_SOURCES})
set_target_properties(pitchlib PROPERTIES LINKER_LANGUAGE "Fortran")

add_executable(PitchAngles calc_pitch_angles.f90)
set_target_properties(PitchAngles PROPERTIES LINKER_LANGUAGE "Fortran") 
add_dependencies(PitchAngles pitchlib)
target_link_libraries(PitchAngles pitchlib)