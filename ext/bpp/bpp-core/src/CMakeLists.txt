# CMake script for Bio++ and Treerecs Library
# Authors:
#   Nicolas COMTE
# Created: 27/07/2017

# File list
set(${PROJECT_NAME}_SOURCE_FILES
  Bpp/App/ApplicationTools.cpp
  Bpp/BppString.cpp
  Bpp/Exceptions.cpp
  Bpp/Graph/GlobalGraph.cpp
  Bpp/Io/BppOParametrizableFormat.cpp
  Bpp/Io/FileTools.cpp
  Bpp/Numeric/AbstractParametrizable.cpp
  Bpp/Numeric/AdaptiveKernelDensityEstimation.cpp
  Bpp/Numeric/AutoParameter.cpp
  Bpp/Numeric/DataTable.cpp
  Bpp/Numeric/NumTools.cpp
  Bpp/Numeric/Parameter.cpp
  Bpp/Numeric/ParameterExceptions.cpp
  Bpp/Numeric/ParameterList.cpp
  Bpp/Numeric/VectorTools.cpp
  Bpp/Numeric/Random/RandomFactory.h
  Bpp/Numeric/Random/RandomTools.cpp
  Bpp/Numeric/Random/Uniform01K.cpp
  Bpp/Text/KeyvalTools.cpp
  Bpp/Text/NestedStringTokenizer.cpp
  Bpp/Text/StringTokenizer.cpp
  Bpp/Text/TextTools.cpp
  Bpp/Utils/AttributesTools.cpp
  )

# Build the static lib
add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCE_FILES})

# Export variables as CACHE entries
string(TOUPPER ${PROJECT_NAME} UPPER_PROJECT_NAME)
set(${UPPER_PROJECT_NAME}_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include dirs" FORCE)
