# Add executable JointSearch
add_executable(${PROJECT_NAME} JointSearch.cpp
  Arguments.cpp
  likelihoods/LibpllEvaluation.cpp
  treeSearch/JointTree.cpp
  treeSearch/Moves.cpp
  treeSearch/Rollbacks.cpp
  treeSearch/NNISearch.cpp
  treeSearch/SPRSearch.cpp
  )

target_link_libraries(${PROJECT_NAME}
  PRIVATE ale)
set_target_properties(${PROJECT_NAME}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  )

target_include_directories(${PROJECT_NAME}
  PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
  )

target_compile_definitions(${PROJECT_NAME}
  PRIVATE "-DVERSION_NUMBER=${PROJECT_VERSION}")

