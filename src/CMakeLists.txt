# Add executable JointSearch


add_executable(${PROJECT_NAME} JointSearch.cpp
  Arguments.cpp
  ParallelContext.cpp
  likelihoods/LibpllEvaluation.cpp
  treeSearch/JointTree.cpp
  treeSearch/Moves.cpp
  treeSearch/Rollbacks.cpp
  treeSearch/NNISearch.cpp
  treeSearch/SPRSearch.cpp
  )

target_link_libraries(${PROJECT_NAME}
  ale)
set_target_properties(${PROJECT_NAME}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  )

if (WITH_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
  if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
      COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()

  if(MPI_LINK_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif()


target_include_directories(${PROJECT_NAME}
  PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
  )

target_compile_definitions(${PROJECT_NAME}
  PRIVATE "-DVERSION_NUMBER=${PROJECT_VERSION}")

