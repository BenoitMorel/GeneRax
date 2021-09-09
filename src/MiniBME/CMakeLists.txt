
set(minibme_SOURCES minibme.cpp 
  MiniBMEArguments.cpp
  MiniBMEOptimizer.cpp
  )
add_program(minibme "${minibme_SOURCES}")
target_link_libraries(minibme mpi-scheduler) 


