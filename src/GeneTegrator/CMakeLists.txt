
set(genetegrator_SOURCES genetegrator.cpp 
  GeneTegratorArguments.cpp
  GTSpeciesTreeOptimizer.cpp
  )
add_program(genetegrator "${genetegrator_SOURCES}")
target_link_libraries(genetegrator mpi-scheduler) 


