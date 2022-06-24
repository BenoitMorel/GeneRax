
set(species_tree_tests_SOURCES species_tree_tests.cpp 
  )
#add_program(rooted_newick_parser_tests "rooted_newick_parser_tests.cpp")
add_program(species_tree_tests "species_tree_tests.cpp")
add_program(pllrooted_tree_tests "pllrooted_tree_tests.cpp")
add_program(pllunrooted_tree_tests "pllunrooted_tree_tests.cpp")
add_program(polytomy_solver_tests "polytomy_solver_tests.cpp")
add_program(polytree_tests "polytree_tests.cpp")
add_program(conditional_clade_tests "conditional_clade_tests.cpp")
add_program(species_splits_tests "species_splits_tests.cpp")
add_program(split_score_tests "split_score_tests.cpp")
add_program_corax(test_coraxlib "test_coraxlib.cpp")
