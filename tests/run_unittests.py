import os
import sys
import subprocess

script_dir = os.path.dirname(os.path.realpath(__file__))
repo_dir = os.path.realpath(os.path.join(script_dir, os.pardir))
species_tree_test = os.path.join(repo_dir, "build", "bin", "species_tree_tests")


def call_test(test_name):
  test_path = os.path.join(repo_dir, "build", "bin", test_name)
  subprocess.check_call([test_path])

call_test("species_tree_tests")
call_test("pllrooted_tree_tests")



