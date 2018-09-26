# CMake script for Bio++ and Treerecs Library
# Authors:
#   Nicolas COMTE
# Created: 27/07/2017

# File list
set(${PROJECT_NAME}_SOURCE_FILES
  Bpp/Seq/Io/NexusTools.cpp
  )

# Build the static lib
add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCE_FILES})
target_include_directories(${PROJECT_NAME}
  PUBLIC ${BPP-CORE_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}
  PUBLIC bpp-core)

# Export variables as CACHE entries
string(TOUPPER ${PROJECT_NAME} UPPER_PROJECT_NAME)
set(${UPPER_PROJECT_NAME}_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include dirs" FORCE)
