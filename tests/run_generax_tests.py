import os
import sys
import subprocess
import shutil
from  distutils.spawn import find_executable

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
REPO_DIR = os.path.realpath(os.path.join(SCRIPT_DIR, os.pardir))
GENERAX = os.path.join(REPO_DIR, "build", "bin", "generax")
FAMILIES_SCRIPT = os.path.join(REPO_DIR, "scripts", "build_families_file.py")
DATA_DIR = os.path.join(REPO_DIR, "data")
OUTPUT = os.path.join(REPO_DIR, "tests", "outputs")

try:
  os.makedirs(OUTPUT)
except:
  pass

def is_mpi_installed():
  return find_executable("mpiexec") is not None

def get_test_name(dataset, with_starting_tree, strategy, model, cores):
  test_name = dataset
  test_name += "_" + str(with_starting_tree)
  test_name += "_" + strategy
  test_name += "_" + model
  test_name += "_" + str(cores)
  return test_name

def reset_dir(directory):
  shutil.rmtree(directory, ignore_errors=True)
  os.makedirs(directory)

def generate_families_file(test_data, with_starting_tree, test_output):
  families_file = os.path.join(test_output, "families.txt")
  command = []
  command.append("python")
  command.append(FAMILIES_SCRIPT)
  command.append(os.path.join(test_data, "alignments"))
  if (with_starting_tree):
    command.append(os.path.join(test_data, "raxml_trees"))
  else:
    command.append("NONE")
  command.append(os.path.join(test_data, "mappings"))
  command.append("GTR")
  command.append(families_file)
  logs_file_path = os.path.join(test_data, "families_script_logs.txt")
  with open(logs_file_path, "w") as writer:
    subprocess.check_call(command, stdout = writer, stderr = writer)
  return families_file

def run_generax(test_data, test_output, families_file, strategy, model, cores):
  command = []
  if (cores > 1):
    command.append("mpiexec")
    command.append("-np")
    command.append(str(cores))
  command.append(GENERAX)
  command.append("-f")
  command.append(families_file)
  command.append("-s")
  command.append(os.path.join(test_data, "speciesTree.newick"))
  command.append("--strategy")
  command.append(strategy)
  command.append("--rec-model")
  command.append(model)
  command.append("-p")
  command.append(os.path.join(test_output, "generax"))
  logs_file_path = os.path.join(test_output, "tests_logs.txt")
  with open(logs_file_path, "w") as writer:
    subprocess.check_call(command, stdout = writer, stderr = writer)


def run_test(dataset, with_starting_tree, strategy, model, cores):
  test_name = get_test_name(dataset, with_starting_tree, strategy, model, cores)
  test_output = os.path.join(OUTPUT, test_name)
  reset_dir(test_output)
  test_data = os.path.join(DATA_DIR, dataset)
  try:
    families_file = generate_families_file(test_data, with_starting_tree, test_output)
    run_generax(test_data, test_output, families_file, strategy, model, cores)
    print("Test " + test_name + ": ok") 
  except:
    print("Test " + test_name + ": FAILED") 
    return False
  return True

dataset_set = ["simulated_2"]
with_starting_tree_set = [True, False]
strategy_set = ["SPR", "EVAL"]
model_set = ["UndatedDL", "UndatedDTL"]
cores_set = [1]
if (is_mpi_installed()):
  cores_set.append(3)

ok = True
for dataset in dataset_set:
  for with_starting_tree in with_starting_tree_set:
    for strategy in strategy_set:
      for model in model_set:
        for cores in cores_set:
          ok = ok and run_test(dataset, with_starting_tree, strategy, model, cores)
if (not ok):
  print("[Error] Some tests failed, please fix them!")
  exit(1)
else:
  print("All tests succeeded!")





