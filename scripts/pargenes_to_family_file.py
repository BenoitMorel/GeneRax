import os
import sys

def get_model_name(model_path):
  full_model = open(model_path).read().split(",")[0]
  model_name = full_model.split("{")[0]
  return model_name

def get_family_to_alignments(pargenes_dir):
  logs = os.path.join(pargenes_dir, "pargenes_logs.txt")
  lines = open(logs).readlines()
  found = False
  pargenes_args = None
  for line in lines:
    if (found):
      pargenes_args = line.split()
      break
    if (line.startswith("ParGenes was called as follow:")):
      found = True
  found = False
  alignments_dir = ""
  for arg in pargenes_args:
    if (found):
      alignments_dir = arg
      break
    if (arg == "-a" or arg == "--alignments-dir"):
      found = True
  print(alignments_dir)
  fam_to_ali = {}
  for alignment in os.listdir(alignments_dir):
    family = alignment.replace(".", "_")
    ali = os.path.abspath(os.path.join(alignments_dir, alignment))
    fam_to_ali[family] = ali
  return fam_to_ali

def export(pargenes_dir, output):
  results_dir = os.path.join(pargenes_dir, "mlsearch_run", "results")
  if (not os.path.isdir(results_dir)):
    print("Error: " + results_dir + " is not a directory, please check your pargenes directory")
    sys.exit(1)
  fam_to_ali = get_family_to_alignments(pargenes_dir)
  writer = open(output, "w")
  writer.write("[FAMILIES]\n")
  for family in os.listdir(results_dir):
    writer.write("- " + family + "\n")
    tree_path = os.path.join(results_dir, family, family + ".raxml.bestTree")
    model_path = os.path.join(results_dir, family, family + ".raxml.bestModel")
    model = get_model_name(model_path)
    # here one could adapt this script to add the path
    # to the mapping file, either from the family name (family)
    # or from the path to the alignment (fam_to_ali[family])
    writer.write("starting_gene_tree = " + tree_path + "\n")
    writer.write("subst_model = " + model + "\n")
    writer.write("alignment = " + fam_to_ali[family] + "\n")

if (__name__ == "__main__"): 
  if (len(sys.argv) != 3): 
    print(" [Error] Invalid syntax")
    print(" Usage: python " + os.path.basename(__file__) + " pargenes_run_dir output_family_file")
    sys.exit(1)
  pargenes_dir = sys.argv[1]
  output = sys.argv[2]
  export(pargenes_dir, output)
