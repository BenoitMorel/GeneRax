
set(asteroid_SOURCES asteroid.cpp 
  AsteroidArguments.cpp
  AsteroidOptimizer.cpp
  )
add_program(asteroid "${asteroid_SOURCES}")
target_link_libraries(asteroid mpi-scheduler) 


