
cmake_minimum_required(VERSION 2.8)

PROJECT(Point)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(Point MACOSX_BUNDLE Point)

if(VTK_LIBRARIES)
  target_link_libraries(Point ${VTK_LIBRARIES})
else()
  target_link_libraries(Point vtkHybrid vtkWidgets)
endif()
