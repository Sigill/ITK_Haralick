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

add_executable(cooccurrence_matrix test/cooccurrence_matrix.cpp)
 
if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(cooccurrence_matrix ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(cooccurrence_matrix ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )

add_executable(filters test/filters.cpp)
 
if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(filters ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(filters ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )

add_executable(image_to_features test/image_to_features.cpp)
 
if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(image_to_features ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(image_to_features ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )

add_executable(image_to_local_features test/image_to_local_features.cpp)
 
if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(image_to_local_features ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(image_to_local_features ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )
