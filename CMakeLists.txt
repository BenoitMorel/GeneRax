cmake_minimum_required(VERSION 3.0.1)

project(JointTreeSearch
        VERSION 0.0)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++14 CXX_STANDARD_14_SUPPORT)
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

set(JOINTSEARCH_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src 
   ${CMAKE_CURRENT_SOURCE_DIR}/ext)

add_subdirectory(ext)
add_subdirectory(src)
