
set(generax_SOURCES generax.cpp 
  GeneRaxArguments.cpp
  GeneRaxSlaves.cpp
  )
add_program(generax "${generax_SOURCES}")
target_link_libraries(generax mpi-scheduler) 

set(generaxslaves_SOURCES GeneRaxSlaves.cpp)
add_program(generaxslaves "${generaxslaves_SOURCES}")
target_compile_definitions(generaxslaves PRIVATE SLAVE_EXEC=1)
