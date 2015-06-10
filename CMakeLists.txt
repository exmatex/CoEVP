cmake_minimum_required(VERSION 2.6)

project(CoEVP)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)

enable_language(CXX)
include(CMakeDetermineCXXCompiler)

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++0x" HAVE_C0X)
if(NOT HAVE_C0X)
    message(FATAL_ERROR "flag '-std=c++0x' does not work")
endif(NOT HAVE_C0X)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# Apple Clang needs to be directed to correct stdlib
if (APPLE)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
  endif()
endif()

########################################################################
# User input options                                                   #
########################################################################
option(BUILD_SHARED_LIBS "Build shared libs" OFF)

########################################################################
#Find external packages
########################################################################
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFile)
foreach(HEADER math.h stdio.h stdlib.h string.h)
  check_include_file(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

include(CheckIncludeFileCXX)
foreach(HEADER algorithm cassert cfloat cmath cstdio cstdlib cstring ctime iostream list map vector)
  check_include_file_cxx(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

set(MATH_LIBRARIES "m" CACHE STRING "math library")
mark_as_advanced( MATH_LIBRARIES )
include(CheckLibraryExists)
foreach(FUNC sqrt)
  check_library_exists(${MATH_LIBRARIES} ${FUNC} "" FOUND_${FUNC}_${MATH_LIBRARIES})
  if(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
    message(FATAL_ERROR "Could not find needed math function - ${FUNC}")
  endif(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
endforeach(FUNC)

######################################
# actual targets                     # 
######################################
include_directories(${CMAKE_BINARY_DIR})
file(GLOB_RECURSE CM_SOURCES CM/*.cc)
file(GLOB_RECURSE CM_HEADERS CM/*.h)
set(CM_INCLUDE_DIRS)
foreach(_HEADER ${CM_HEADERS})
	get_filename_component(_DIR ${_HEADER} DIRECTORY)
	list(APPEND CM_INCLUDE_DIRS ${_DIR})
endforeach()
list(REMOVE_DUPLICATES CM_INCLUDE_DIRS)
include_directories(${CM_INCLUDE_DIRS}
	CM/src/adaptive_sampling/interpolation_database
	CM/src/adaptive_sampling/interpolation
	CM/src/adaptive_sampling/mtl_headers
	CM/src/adaptive_sampling/database
	CM/src/adaptive_sampling/utils)

add_definitions(-DDBL_SNAN_IS_BROKEN -DFLT_SNAN_IS_BROKEN)
add_library(cm ${CM_SOURCES})

file(GLOB LULESH_SOURCES LULESH/*_main.cc)
add_executable(lulesh ${LULESH_SOURCES})

file(GLOB LULESH_CORE_SOURCES LULESH/*.cc)
list(REMOVE_ITEM LULESH_CORE_SOURCES ${LULESH_SOURCES})
add_library(lulesh_core ${LULESH_CORE_SOURCES})

set_target_properties(lulesh PROPERTIES RUNTIME_OUTPUT_DIRECTORY LULESH)
target_link_libraries(lulesh lulesh_core cm ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
