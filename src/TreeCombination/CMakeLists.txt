
set(treecombination_SOURCES treecombination.cpp 
  )
add_program(treecombination "${treecombination_SOURCES}")
target_link_libraries(treecombination mpi-scheduler) 

