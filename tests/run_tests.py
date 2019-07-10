import os
import sys
import subprocess

script_dir = os.path.dirname(os.path.realpath(__file__))
repo_dir = os.path.realpath(os.path.join(script_dir, os.pardir))
executable = os.path.join(repo_dir, "build", "bin", "JointSearch")
data_dir = os.path.join(repo_dir, "data")
output_root = os.path.join(repo_dir, "tests", "outputs")
executable = os.path.join(repo_dir, "build", "bin", "JointSearch")
try:
  os.makedirs(output_root)
except:
  pass

def get_true_hash(family_dir):
  return int(open(os.path.join(family_dir, "trueHash.txt")).readlines()[0])

def get_infered_hash(output_prefix):
  lines = open((output_prefix +".stats")).readlines()
  for line in lines:
    if line.startswith("hash"):
      return int(line.split(" ")[1][:-1])

def test_family_with_parameters(family_dir, with_gene_tree):
  if (with_gene_tree):
    print("  with raxml starting gene tree")
  else:
    print ("  with random starting gene tree")
  alignment = os.path.join(family_dir, "alignment.msa")
  speciesTree = os.path.join(family_dir, "speciesTree.newick")
  mapping = os.path.join(family_dir, "mapping.link")
  true_hash = get_true_hash(family_dir)
  output_prefix = os.path.join(output_root, "test")
  gene_tree = os.path.join(family_dir, "raxmlGeneTree.newick")
  command = []
  command.append(executable)
  command.append("-s")
  command.append(speciesTree)
  command.append("-a")
  command.append(alignment)
  command.append("-m")
  command.append(mapping)
  command.append("-p")
  command.append(output_prefix)
  command.append("--strategy")
  command.append("SPR")
  if (with_gene_tree):
    command.append("-g")
    command.append(gene_tree)
  print("Executing " + " ".join(command))
  FNULL = open(os.devnull, 'w')
  subprocess.check_call(command, stdout = FNULL)
  infered_hash = get_infered_hash(output_prefix)
  print(infered_hash)
  print(get_infered_hash(output_prefix))
  print(true_hash)
  assert(infered_hash == true_hash)
  print("ok!")

def test_family(family_dir):
  print("Test family " + os.path.basename(family_dir) + ":")
  test_family_with_parameters(family_dir, True)
  test_family_with_parameters(family_dir, False)

family_dir = os.path.join(data_dir, "simulated_1")
test_family(family_dir)

