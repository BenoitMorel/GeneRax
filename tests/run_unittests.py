import os
import sys
import subprocess

script_dir = os.path.dirname(os.path.realpath(__file__))
repo_dir = os.path.realpath(os.path.join(script_dir, os.pardir))
species_tree_test = os.path.join(repo_dir, "build", "bin", "species_tree_tests")



subprocess.check_call([species_tree_test])


