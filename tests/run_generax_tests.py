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

def generate_families_file_data(test_data, with_starting_tree, test_output):
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

def generate_families_file(test_output, alignments = "NONE", starting_trees = "NONE", mappings = "NONE", subst_model = "NONE"):
  families_file = os.path.join(test_output, "families.txt")
  command = []
  command.append("python")
  command.append(FAMILIES_SCRIPT)
  command.append(alignments)
  command.append(starting_trees)
  command.append(mappings)
  command.append(subst_model)
  command.append(families_file)
  logs_file_path = os.path.join(test_output, "families_script_logs.txt")
  with open(logs_file_path, "w") as writer:
    subprocess.check_call(command, stdout = writer, stderr = writer)
  return families_file


def run_reconciliation(species_tree, families_file, model, test_output, cores):
  command = []
  if (cores > 1):
    command.append("mpiexec")
    command.append("-np")
    command.append(str(cores))
  command.append(GENERAX)
  command.append("-f")
  command.append(families_file)
  command.append("-s")
  command.append(species_tree)
  command.append("--rec-model")
  command.append(model)
  command.append("--do-not-optimize-gene-trees")
  command.append("--reconcile")
  command.append("--dup-rate")
  command.append("0.2")
  command.append("--loss-rate")
  command.append("0.2")
  command.append("--transfer-rate")
  command.append("0.1")
  command.append("-p")
  command.append(os.path.join(test_output, "generax"))
  logs_file_path = os.path.join(test_output, "tests_logs.txt")
  with open(logs_file_path, "w") as writer:
    subprocess.check_call(command, stdout = writer, stderr = writer)

def is_string_in_file(string, file_name):
  return string in open(file_name).read()

def count_string_in_file(string, file_name):
  return open(file_name).read().count(string)

def check_reconciliation(test_output, model):
  reconciliations_path = os.path.join(test_output, "generax", "reconciliations")
  nhx_dup_A = os.path.join(reconciliations_path, "gene_dup_A_reconciliated.nhx")
  nhx_dup_AB = os.path.join(reconciliations_path, "gene_dup_AB_reconciliated.nhx")
  nhx_transfer_A_D = os.path.join(reconciliations_path, "gene_transfer_A_D_reconciliated.nhx")
  if (not is_string_in_file("[&&NHX:S=A:D=Y:H=N:B=0]", nhx_dup_A)):
    print("Failed to infer a duplication in species A (" + nhx_dup_A + ")")
    return False
  if (count_string_in_file("=Y", nhx_dup_A) != 1):
    print("Inferred to many events in " + nhx_dup_A)
    return False;
  if (not is_string_in_file("[&&NHX:S=AB:D=Y:H=N:B=0]", nhx_dup_AB)):
    print("Failed to infer a duplication in species AB (" + nhx_dup_AB + ")")
    return False
  if (count_string_in_file("=Y", nhx_dup_AB) != 1):
    print("Inferred to many events in " + nhx_dup_AB)
    return False;
  if (model == "UndatedDTL"):
    if (not is_string_in_file("[&&NHX:S=A:D=N:H=Y@A@D:B=0]", nhx_transfer_A_D)):
      print("Failed to infer a transfer from A to D (" + nhx_transfer_A_D + ")")
      return False
    if (count_string_in_file("=Y", nhx_transfer_A_D) != 1):
      print("Inferred to many events in " + nhx_transfer_A_D)
      return False;
  return True

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
    families_file = generate_families_file_data(test_data, with_starting_tree, test_output)
    run_generax(test_data, test_output, families_file, strategy, model, cores)
    print("Test " + test_name + ": ok") 
  except:
    print("Test " + test_name + ": FAILED") 
    return False
  return True

def run_reconciliation_test(cores, model):
  test_name = "reconciliation_" + model
  test_output = os.path.join(OUTPUT, test_name)
  reset_dir(test_output)
  reconciliation_dir = os.path.join(DATA_DIR, "reconciliation")
  gene_trees = os.path.join(reconciliation_dir, "gene_trees")
  species_tree = os.path.join(reconciliation_dir, "ABCD_species.newick")
  families_file = generate_families_file(test_output, starting_trees = gene_trees)
  ok = False
  run_reconciliation(species_tree, families_file, model, test_output, cores)
  ok = check_reconciliation(test_output, model)
  if (ok):
    print("Test " + test_name + ": ok")
  else:
    print("Test " + test_name + ": FAILED")
  return ok

dataset_set = ["simulated_2"]
with_starting_tree_set = [True, False]
strategy_set = ["SPR", "EVAL"]
model_set = ["UndatedDL", "UndatedDTL"]
cores_set = [1]
if (is_mpi_installed()):
  cores_set.append(3)

all_ok = True
all_ok = all_ok and run_reconciliation_test(1, "UndatedDTL")
all_ok = all_ok and run_reconciliation_test(1, "UndatedDL")
for dataset in dataset_set:
  for with_starting_tree in with_starting_tree_set:
    for strategy in strategy_set:
      for model in model_set:
        for cores in cores_set:
          all_ok = all_ok and run_test(dataset, with_starting_tree, strategy, model, cores)

if (not all_ok):
  print("[Error] Some tests failed, please fix them!")
  exit(1)
else:
  print("All tests succeeded!")





