cmake_minimum_required(VERSION 2.6)
 
project(HaralickComputer)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
include_directories(BEFORE "${CMAKE_SOURCE_DIR}/include")
include_directories(BEFORE "${CMAKE_SOURCE_DIR}/test")
 
add_executable(haralick_iterations test/haralick_iterations.cpp)
 
if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(haralick_iterations ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(haralick_iterations ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )

#target_link_libraries(basic ${ITK_LIBRARIES})
