
cmake_minimum_required(VERSION 3.0.1)

project(boost)

# Export variables as CACHE entries
set(BOOST_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include dirs" FORCE)
