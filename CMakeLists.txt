cmake_minimum_required(VERSION 3.0.1)

project(JointSearchProject
        VERSION 0.0)

set(JOINTSEARCH_BUILD_AS_LIBRARY OFF CACHE BOOL "Build JointSearch as shared library (instead of stand-alone executable)")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

if(JOINTSEARCH_BUILD_AS_LIBRARY)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DJOINTSEARCH_BUILD_AS_LIB=ON")
endif()

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++14 CXX_STANDARD_14_SUPPORT -g)
if(CXX_STANDARD_14_SUPPORT)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++14 support.")
endif()


add_subdirectory(ext)

if (NOT DISABLE_MPI)
  include(FindMPI)
  if(MPI_FOUND)
    set(WITH_MPI TRUE)
    add_definitions(-DWITH_MPI)
    message(STATUS "WITH MPI ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()


set(JOINTSEARCH_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src/core
   ${CMAKE_CURRENT_SOURCE_DIR}/ext
   ${PLLMOD_EXPORTED_INCLUDE_PATH}
   ${MPI_INCLUDE_PATH})
  

add_subdirectory(src)

