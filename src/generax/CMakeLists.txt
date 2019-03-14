
set(generax_SOURCES generax.cpp GeneRaxArguments.cpp)
set(generax_optimize_gene_trees_SOURCES generax_optimize_gene_trees.cpp GeneRaxArguments.cpp)
add_program(generax "${generax_SOURCES}" False)
add_program(generax_optimize_gene_trees "${generax_optimize_gene_trees_SOURCES}" True) 
target_link_libraries(generax mpi-scheduler) 



