cmake_minimum_required(VERSION 2.8)
enable_language(CXX)
project(vtkPlot)

add_definitions( -std=gnu++11 -Wall )

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
    -Wextra
    -Werror
    -Wpedantic
    -Wno-error=unused-result
    -Wno-error=unused-parameter
    -Wno-error=unused-variable
    -Wno-error=unused-function
    )
endif ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
      


file(GLOB_RECURSE HDRS
  ${vtkPlot_SOURCE_DIR}/src/*.hpp
  ${vtkPlot_SOURCE_DIR}/src/*.h)
file(GLOB_RECURSE SRCS
  ${vtkPlot_SOURCE_DIR}/src/*.cpp
  ${vtkPlot_SOURCE_DIR}/src/*.c)

add_executable( vtkPlot ${HDRS} ${SRCS} )
  
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

target_link_libraries( vtkPlot ${VTK_LIBRARIES})
#if( VTK_LIBRARIES )
#  target_link_libraries( vtkPlot ${VTK_LIBRARIES})
#else()
#  target_link_libraries( vtkPlot vtkHybrid vtkWidgets)
#endif()
