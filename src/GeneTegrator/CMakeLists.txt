
set(genetegrator_SOURCES genetegrator.cpp 
  GeneTegratorArguments.cpp
  GTSpeciesTreeOptimizer.cpp
  TrimFamilies.cpp
  )
add_program(genetegrator "${genetegrator_SOURCES}")
target_link_libraries(genetegrator PRIVATE mpi-scheduler) 


