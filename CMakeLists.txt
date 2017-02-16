cmake_minimum_required(VERSION 2.8)
 
PROJECT(create_atlas)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


add_executable(create_atlas create_atlas.cxx)

 
if(VTK_LIBRARIES)
  target_link_libraries(create_atlas ${VTK_LIBRARIES})
 
endif()