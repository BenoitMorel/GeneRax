# Add executable JointSearch

#add_executable(${PROJECT_NAME} JointSearch.cpp


set(JOINTSEARCH_SOURCES JointSearch.cpp
  Arguments.cpp
  ParallelContext.cpp
  Logger.cpp
  likelihoods/LibpllEvaluation.cpp
  likelihoods/ALEEvaluation.cpp
  likelihoods/ale/UndatedDLModel.cpp
  likelihoods/dated_ale/exODT.cpp
  likelihoods/dated_ale/ALE.cpp
  likelihoods/dated_ale/model.cpp
  treeSearch/JointTree.cpp
  treeSearch/Moves.cpp
  treeSearch/Rollbacks.cpp
  treeSearch/SearchUtils.cpp
  treeSearch/NNISearch.cpp
  treeSearch/SPRSearch.cpp
  )

if (JOINTSEARCH_BUILD_AS_LIBRARY)
  message(STATUS "${PROJECT_NAME} will be built as a library.")
  add_library(${PROJECT_NAME} MODULE ${JOINTSEARCH_SOURCES})
else()
  add_executable(${PROJECT_NAME} ${JOINTSEARCH_SOURCES})
endif()


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
      COMPILE_FLAGS "-g ${MPI_COMPILE_FLAGS}")
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

if (BUILD_AS_LIBRARY AND NOT APPLE)
  target_link_libraries(${PROJECT_NAME} gcc_s gcc)
endif()
