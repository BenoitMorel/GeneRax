
set(lightsearch_SOURCES lightsearch.cpp 
  )
add_program(lightsearch "${lightsearch_SOURCES}")
target_link_libraries(lightsearch mpi-scheduler) 

