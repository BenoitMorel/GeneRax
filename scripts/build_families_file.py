import os
import sys

"""
Script that generates a family file for GeneRax.
"""

def get_family(filename):
  return filename.split(".")[0]

def join_abs(directory, filename):
  return os.path.abspath(os.path.join(directory, filename))

def build_families(alignments_dir, trees_dir, mappings_dir, model, output_file):
  families = {}
  trees = {}
  mappings = {}
  if (trees_dir != "NONE"):
    for tree in os.listdir(trees_dir):
      trees[get_family(tree)] = join_abs(trees_dir, tree)
  if (mappings_dir != "NONE"):
    for mapping in os.listdir(mappings_dir):
      mappings[get_family(mapping)] = join_abs(mappings_dir, mapping)
  for alignment in os.listdir(alignments_dir):
    family = get_family(alignment) 
    cell = {}
    cell["alignment"] = join_abs(alignments_dir, alignment)
    if (family in trees):
      cell["starting_gene_tree"] = trees[family]
    if (family in mappings):
      cell["mapping"] = mappings[family]
    cell["subst_model"] = model
    if (family in families):
      print("[Error]: family " + family + " was found twice in the alignments directory.")
      exit(1)
    families[family] = cell
  
  with open(output_file, "w") as writer:
    writer.write("[FAMILIES]\n")
    for family in families:
      cell = families[family]
      writer.write("- " + family + "\n")
      for param in cell:
        writer.write(param + " = " + cell[param] + "\n")
  print("Output families file" + output_file)

if (__name__ == "__main__"): 
  if (len(sys.argv) != 6): 
     print(" [Error] Invalid syntax")
     print(" Usage: python " + os.path.basename(__file__) + " alignments_dir trees_dir mappings_dir subst_model output_file")
     print(" set trees_dir and/or mappings_dir to NONE if you don't have any.")
     print(" Example: python "  + os.path.basename(__file__) + " /home/myalignments/ NONE NONE GTR /home/families.txt")
     exit(1)
  alignments_dir = sys.argv[1]
  trees_dir = sys.argv[2]
  mappings_dir = sys.argv[3]
  model = sys.argv[4]
  output_file = sys.argv[5]
  build_families(alignments_dir, trees_dir, mappings_dir, model, output_file)




