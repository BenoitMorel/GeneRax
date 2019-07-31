
set(bench_rec_likelihood_SOURCES bench_rec_likelihood.cpp 
  )
add_program(bench_rec_likelihood "${bench_rec_likelihood_SOURCES}")
target_link_libraries(bench_rec_likelihood mpi-scheduler) 

