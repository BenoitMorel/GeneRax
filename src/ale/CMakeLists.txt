
set(alegenerax_SOURCES ale.cpp 
  AleArguments.cpp
  AleEvaluator.cpp
  AleOptimizer.cpp
  TrimFamilies.cpp
  )
add_program(alegenerax "${alegenerax_SOURCES}")
target_link_libraries(alegenerax PRIVATE mpi-scheduler) 


