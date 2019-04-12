


add_subdirectory(core)
  
set(PLL_LIBRARIES 
    pllmodalgorithm_static
    pllmodbinary_static
    pllmodmsa_static
    pllmodoptimize_static
    pllmodtree_static
    pllmodutil_static
   )




macro(add_program program_name program_sources)
  if (JOINTSEARCH_BUILD_AS_LIBRARY)
    add_library(${program_name} SHARED ${program_sources})
  else()
    add_executable(${program_name} ${program_sources})
  endif()
  target_link_libraries(${program_name} 
    jointsearch-core
    ${MPI_LIBRARIES}
    ${PLL_LIBRARIES}
    )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  if (JOINTSEARCH_BUILD_AS_LIBRARY AND NOT APPLE)
    target_link_libraries(${program_name} gcc_s gcc)
  endif()

  set_target_properties(${program_name}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
  set_target_properties(${program_name} 
    PROPERTIES
    COMPILE_FLAGS "-g ${MPI_COMPILE_FLAGS}"
    )
  set_target_properties(${program_name} 
    PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}"
    )

  target_include_directories(${program_name}
    PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
    )
endmacro()



add_subdirectory(JointSearch)
add_subdirectory(generax)
add_subdirectory(SpeciesRax)


