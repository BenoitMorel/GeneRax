
set(njrax_SOURCES njrax.cpp 
  )
add_program(njrax "${njrax_SOURCES}")
target_link_libraries(njrax mpi-scheduler) 


