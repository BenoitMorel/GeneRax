
set(speciesrax_SOURCES speciesrax.cpp 
  SpeciesRaxArguments.cpp
  )
add_program(speciesrax "${speciesrax_SOURCES}" False)
target_link_libraries(speciesrax mpi-scheduler) 



