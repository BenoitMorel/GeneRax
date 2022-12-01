
set(speciessplitsearch_SOURCES speciessplitsearch.cpp 
  SpeciesSplitArguments.cpp
  SpeciesSplitSpeciesTreeOptimizer.cpp
  )
add_program(speciessplitsearch "${speciessplitsearch_SOURCES}")
target_link_libraries(speciessplitsearch PRIVATE mpi-scheduler) 


