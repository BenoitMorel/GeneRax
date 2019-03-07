


add_subdirectory(core)


# SOURCES
set(JointSearchLight_SOURCES JointSearchLight.cpp)

# MAIN LIBRARY OR EXECUTABLE
if (JOINTSEARCH_BUILD_AS_LIBRARY)
  add_library(JointSearchLight MODULE ${JointSearchLight_SOURCES})
else()
  add_executable(JointSearchLight ${JointSearchLight_SOURCES})
endif()


# LINK LIBRARIES
set(PLL_LIBRARIES 
    pllmodalgorithm_static
    pllmodbinary_static
    pllmodmsa_static
    pllmodoptimize_static
    pllmodtree_static
    pllmodutil_static
   )
target_link_libraries(JointSearchLight 
  jointsearch-core
  ${MPI_LIBRARIES}
  ${PLL_LIBRARIES}
  )
if (JOINTSEARCH_BUILD_AS_LIBRARY AND NOT APPLE)
  target_link_libraries(JointSearchLight gcc_s gcc)
endif()

# PROPERTIES
set_target_properties(JointSearchLight
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  )
set_target_properties(JointSearchLight 
  PROPERTIES
  COMPILE_FLAGS "-g ${MPI_COMPILE_FLAGS}"
  )
set_target_properties(JointSearchLight 
  PROPERTIES
  LINK_FLAGS "${MPI_LINK_FLAGS}"
  )

# INCLUDE
target_include_directories(JointSearchLight
  PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
  )
