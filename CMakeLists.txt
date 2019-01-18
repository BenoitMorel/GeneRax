cmake_minimum_required(VERSION 3.0.1)

project(JointSearch
        VERSION 0.0)

set(JOINTSEARCH_BUILD_AS_LIBRARY OFF CACHE BOOL "Build JointSearch as shared library (instead of stand-alone executable)")

if(JOINTSEARCH_BUILD_AS_LIBRARY)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DJOINTSEARCH_BUILD_AS_LIB=ON -fPIC")
endif()

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++14 CXX_STANDARD_14_SUPPORT -g)
if(CXX_STANDARD_14_SUPPORT)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++14 support.")
endif()


include(FindMPI)
if(MPI_FOUND)
  set(WITH_MPI TRUE)
  add_definitions(-DWITH_MPI)
endif()


add_subdirectory(ext)

set(JOINTSEARCH_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src 
   ${CMAKE_CURRENT_SOURCE_DIR}/ext
   ${PLLMOD_EXPORTED_INCLUDE_PATH})

add_subdirectory(src)
