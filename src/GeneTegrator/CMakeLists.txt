
set(genetegrator_SOURCES genetegrator.cpp 
  GeneTegratorArguments.cpp
  )
add_program(genetegrator "${genetegrator_SOURCES}")
#target_link_libraries(genetegrator mpi-scheduler) 


