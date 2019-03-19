
set(generax_SOURCES generax.cpp 
  GeneRaxArguments.cpp
  GeneRaxSlaves.cpp
  )
add_program(generax "${generax_SOURCES}" False)
target_link_libraries(generax mpi-scheduler) 



