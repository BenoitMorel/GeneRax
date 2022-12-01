
set(generax_SOURCES generax.cpp 
  GeneRaxCore.cpp
  GeneRaxInstance.cpp
  GeneRaxArguments.cpp
  )
add_program(generax "${generax_SOURCES}")
target_link_libraries(generax PRIVATE mpi-scheduler) 

