cmake_minimum_required( VERSION 2.8 )
enable_language( CXX )
project( vtkPlot )

set( CMAKE_C_COMPILER "g++" )
add_definitions( -std=gnu++14 )

if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
  add_definitions(
    -O2
    -mtune=native
    -march=native
    )
endif ( CMAKE_BUILD_TYPE STREQUAL "Release" )
if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  add_definitions(
    -Og
    -ggdb
    -Wall
    -Wextra
    )
endif ( CMAKE_BUILD_TYPE STREQUAL "Debug" )

file(GLOB_RECURSE HDRS
  ${vtkPlot_SOURCE_DIR}/src/*.hpp
  ${vtkPlot_SOURCE_DIR}/src/*.h)
file(GLOB_RECURSE SRCS
  ${vtkPlot_SOURCE_DIR}/src/*.cpp
  ${vtkPlot_SOURCE_DIR}/src/*.c)

configure_file( configuration/root_directory.h.in configuration/root_directory.h )
include_directories( ${CMAKE_BINARY_DIR}/configuration )

add_executable( ${CMAKE_PROJECT_NAME} ${HDRS} ${SRCS} )
  
find_package( VTK 6.0 REQUIRED )
include( ${VTK_USE_FILE} )

target_link_libraries( ${CMAKE_PROJECT_NAME}
  ${VTK_LIBRARIES}
  )
