cmake_minimum_required(VERSION 2.6)
 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
project(ITK_Haralick)

SET(ITK_HARALICK_MAJOR_VERSION 0)
SET(ITK_HARALICK_MINOR_VERSION 0)
SET(ITK_HARALICK_PATCH_VERSION 0)
SET(ITK_HARALICK_VERSION "${ITK_HARALICK_MAJOR_VERSION}.${ITK_HARALICK_MINOR_VERSION}.${ITK_HARALICK_PATCH_VERSION}")

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
include_directories(BEFORE "${CMAKE_SOURCE_DIR}/include")
include_directories(BEFORE "${CMAKE_SOURCE_DIR}/test")

#FIND_PACKAGE(Valgrind REQUIRED)
#INCLUDE_DIRECTORIES(${VALGRIND_INCLUDE_DIR})

SET(INSTALL_DIR_MESSAGE "(If relative, will be appended to CMAKE_INSTALL_PREFIX)" CACHE INTERNAL "Internal message")
SET(INSTALL_DIR_INCLUDE "include" CACHE PATH "Installation directory for header files ${INSTALL_DIR_MESSAGE}")
SET(INSTALL_DIR_CMAKE "lib/cmake" CACHE PATH "Installation directory for cmake files ${INSTALL_DIR_MESSAGE}")
foreach(p INCLUDE CMAKE)
  set(var INSTALL_DIR_${p})
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
  MARK_AS_ADVANCED(${var})
endforeach()

add_executable(whole_image_to_local_glcm_algorithm test/whole_image_to_local_glcm_algorithm.cpp)
add_executable(GreyLevelCooccurrenceMatrix test/GreyLevelCooccurrenceMatrix.cpp)
add_executable(ScalarImageToHaralickTextureFeaturesImage test/ScalarImageToHaralickTextureFeaturesImage.cpp)
add_executable(ScalarImageToHaralickTextureFeaturesImage_statistics test/ScalarImageToHaralickTextureFeaturesImage_statistics.cpp)
add_executable(GLCMImageCalculator test/GLCMImageCalculator.cpp)
add_executable(HaralickFeaturesGLCMCalculator test/HaralickFeaturesGLCMCalculator.cpp)

if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(whole_image_to_local_glcm_algorithm ITKReview ${ITK_LIBRARIES})
  target_link_libraries(GreyLevelCooccurrenceMatrix ITKReview ${ITK_LIBRARIES})
  target_link_libraries(ScalarImageToHaralickTextureFeaturesImage ITKReview ${ITK_LIBRARIES})
  target_link_libraries(ScalarImageToHaralickTextureFeaturesImage_statistics ITKReview ${ITK_LIBRARIES})
  target_link_libraries(GLCMImageCalculator ITKReview ${ITK_LIBRARIES})
  target_link_libraries(HaralickFeaturesGLCMCalculator ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(whole_image_to_local_glcm_algorithm ${ITK_LIBRARIES})
  target_link_libraries(GreyLevelCooccurrenceMatrix ${ITK_LIBRARIES})
  target_link_libraries(ScalarImageToHaralickTextureFeaturesImage ${ITK_LIBRARIES})
  target_link_libraries(ScalarImageToHaralickTextureFeaturesImage_statistics ${ITK_LIBRARIES})
  target_link_libraries(GLCMImageCalculator ${ITK_LIBRARIES})
  target_link_libraries(HaralickFeaturesGLCMCalculator ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )

INSTALL( FILES
  include/IncrementalWeightedMeanVarianceComputer.h
  include/itkGLCMImageCalculator.h
  include/itkGLCMImageCalculator.hxx
  include/itkGreyLevelCooccurrenceMatrix.h
  include/itkGreyLevelCooccurrenceMatrix.hxx
  include/itkHaralickFeaturesGLCMCalculator.h
  include/itkHaralickFeaturesGLCMCalculator.hxx
  include/itkScalarImageToHaralickTextureFeaturesImageFilter.h
  DESTINATION "${INSTALL_DIR_INCLUDE}")


# Write the file used to match a specific version of the library
CONFIGURE_FILE(cmake/ITK_HaralickConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/ITK_HaralickConfigVersion.cmake" @ONLY)
INSTALL(FILES "${PROJECT_BINARY_DIR}/ITK_HaralickConfigVersion.cmake" DESTINATION ${INSTALL_DIR_CMAKE})

configure_file(cmake/ITK_HaralickBuildTreeSettings.cmake.in "${PROJECT_BINARY_DIR}/ITK_HaralickBuildTreeSettings.cmake" @ONLY)

file(RELATIVE_PATH CONF_REL_INCLUDE_DIR "${INSTALL_DIR_CMAKE}" "${INSTALL_DIR_INCLUDE}")
CONFIGURE_FILE(cmake/ITK_HaralickConfig.cmake.in "${PROJECT_BINARY_DIR}/ITK_HaralickConfig.cmake" @ONLY)
INSTALL(FILES "${PROJECT_BINARY_DIR}/ITK_HaralickConfig.cmake" DESTINATION ${INSTALL_DIR_CMAKE})
