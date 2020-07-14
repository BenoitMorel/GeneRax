import os
import sys
import shutil
import random 
import string
import subprocess

runs_number = 50

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
REPO_DIR = os.path.realpath(os.path.join(SCRIPT_DIR, os.pardir))
GENERAX = os.path.join(REPO_DIR, "build", "bin", "generax")
FAMILIES_SCRIPT = os.path.join(REPO_DIR, "scripts", "build_families_file.py")
DATA_DIR = os.path.join(REPO_DIR, "data")
OUTPUT = os.path.join(REPO_DIR, "tests", "outputs")

MIN_LABEL_SIZE = 2
MAX_LABEL_SIZE = 20
MIN_TAXA_NUMBER_S = 4
MAX_TAXA_NUMBER_S = 20
MIN_TAXA_NUMBER_G = 4
MIN_TAXA_NUMBRE_G = 20
MIN_SITES = 5
MAX_SITES = 50
MIN_FAMILIES = 1
MAX_FAMILIES = 5

def mkdir(p):
  try:
    os.mkdir(p)
  except:
    pass


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

class Node(object):
  def __init__(self):
    self.children = []
    self.label = ""
    self.parent = None

  def is_leaf(self):
    return len(self.children) == 0

  def is_root(self):
    return self.parent == None

  def add_child(self, child):
    self.children.append(child)
    child.parent = self

  def remove_child(self, child):
    assert(child in self.children)
    self.children.remove(child)
    child.parent = None

  def to_string(self):
    res = ""
    if (not self.is_leaf()):
      children_strings = []
      for child in self.children:
        children_strings.append(child.to_string())
      res += "(" + ",".join(children_strings) + ")"
    res += self.label
    return res

class Tree(object):
  def __init__(self, taxa_set, rooted):
    self.root = None
    self.nodes = []
    taxa_list = list(taxa_set)
    random.shuffle(taxa_list)
    for taxon in taxa_list:
      self.add_leaf_rand(taxon)
    if (not rooted):
      left = self.root.children[0]
      if (left.is_leaf()):
        self.root.children[0] = self.root.children[1]
        self.root.children[1] = left
        left = self.root.children[0]
      assert(not left.is_leaf())
      left_left = left.children[0]
      left_right = left.children[1]
      self.root.remove_child(left)
      self.root.add_child(left_left)
      self.root.add_child(left_right)

  def create_node(self):
    node = Node()
    self.nodes.append(node)
    return node

  def get_rand_node(self):
    if (len(self.nodes) != 0):
      return random.choice(self.nodes)
    else:
      return None

  def add_leaf_rand(self, label):
    rand_node = self.get_rand_node()
    leaf = self.create_node()
    leaf.label = label
    if (self.root == None):
      self.root = leaf
      return
    parent_rand_node = rand_node.parent
    new_parent = self.create_node()
    if (parent_rand_node != None):
      parent_rand_node.remove_child(rand_node)
      parent_rand_node.add_child(new_parent)
    else:
      self.root = new_parent
    new_parent.add_child(rand_node)
    new_parent.add_child(leaf)


  def to_string(self):
    return self.root.to_string() + ";"

def get_rand_bool():
  return random.choice([True, False])

def get_rand_label(min_size, max_size):
  size = random.randint(min_size, max_size)
  return ''.join(random.choice(string.ascii_lowercase + string.ascii_uppercase + string.digits) for _ in range(size))

def get_rand_taxa_set(taxa_number, prefix = ""):
  taxa_set = set()
  for i in range(0, taxa_number):
    taxa_set.add(prefix + get_rand_label(MIN_LABEL_SIZE, MAX_LABEL_SIZE))
  return taxa_set

def get_rand_taxa_number_s():
  return random.randint(MIN_TAXA_NUMBER_S, MAX_TAXA_NUMBER_S) 

def get_rand_taxa_number_g():
  return random.randint(MIN_TAXA_NUMBER_S, MAX_TAXA_NUMBER_S) 

def get_rand_families_number():
  return random.randint(MIN_FAMILIES, MAX_FAMILIES)

def get_rand_rtree(taxa_set):
  tree = Tree(taxa_set, True)
  return tree.to_string()

def get_rand_utree(taxa_set):
  tree = Tree(taxa_set, False)
  return tree.to_string()

def get_rand_rec_model():
  models = ["UndatedDL", "UndatedDTL"]#, "ParsimonyDL"]
  return random.choice(models)

def get_rand_subst_model_dna():
  models = ["GTR", "JC", "K81uf", "TIM3uf"]
  res = random.choice(models)
  if (get_rand_bool()):
    res += "+G"
  return res

def write_rand_map(species_taxa_set, gene_taxa_set,  output):
  species_list = list(species_taxa_set)
  with open(output, "w") as writer:
    for gene in gene_taxa_set:
      species = random.choice(species_list)
      writer.write(gene + " " + species + "\n")
  
def write_rand_msa(gene_taxa_set, sites, output):
  alphabet = "ACGT-"
  with open(output, "w") as writer:
    for gene in gene_taxa_set:
      writer.write(">" + gene + "\n")
      writer.write(''.join(random.choice(alphabet) for _ in range(sites)))
      writer.write("\n")


def get_rand_per_family_gene_names(family_names):
  res = {}
  for family in family_names:
    genes_number = get_rand_taxa_number_g()
    res[family] = get_rand_taxa_set(genes_number, "gene_")
  return res

class Runner():
  def __init__(self, seed):
    self.ok = False
    self.test_dir = os.path.join(OUTPUT, "randomized_" + str(seed))
    self.destroy_dir()
    mkdir(self.test_dir)
    self.seed = seed
    random.seed(seed)
    self.taxa = get_rand_taxa_set(get_rand_taxa_number_s())
    self.family_names = get_rand_taxa_set(get_rand_families_number(), "fam_") 
    self.per_family_gene_names = get_rand_per_family_gene_names(self.family_names)
    self.cores = random.randint(1, 4)
    # create species tree
    self.species_tree = self._get_path("species_tree.newick")
    open(self.species_tree, "w").write(get_rand_rtree(self.taxa))
     
    # create alignments, mappings and gene trees
    self.input_alignments = self._get_path("alignments")
    self.input_starting_trees = self._get_path("starting_trees")
    self.input_mappings = self._get_path("mappings")
    mkdir(self.input_alignments)
    mkdir(self.input_mappings)
    mkdir(self.input_starting_trees)
    subst_model = get_rand_subst_model_dna()
    for family in self.per_family_gene_names:
      genes = self.per_family_gene_names[family]
      sites = random.randint(MIN_SITES, MAX_SITES)
      ali_out = os.path.join(self.input_alignments, family + ".fasta")
      write_rand_msa(genes, sites, ali_out)
      mapping_out = os.path.join(self.input_mappings, family + ".txt")
      write_rand_map(self.taxa, genes, mapping_out)
      tree_out = os.path.join(self.input_starting_trees, family + ".newick")
      open(tree_out, "w").write(get_rand_utree(genes))
    self.families_file = generate_families_file(self.test_dir, alignments = self.input_alignments, starting_trees = self.input_starting_trees, mappings = self.input_mappings, subst_model = subst_model)
    
    self.rec_model = get_rand_rec_model()
    self.reconcile = get_rand_bool()
    self.reconcile_samples = 0
    if (self.reconcile):
      if (get_rand_bool()):
        self.reconcile_samples = random.randint(1, 5)

  def _get_path(self, p):
    return os.path.join(self.test_dir, p)

  def run(self):
    command = []
    if (self.cores > 1):
      command.append("mpiexec")
      command.append("-np")
      command.append(str(self.cores))
    command.append(GENERAX)
    command.append("-f")
    command.append(self.families_file)
    command.append("-s")
    command.append(self.species_tree)
    command.append("--rec-model")
    command.append(self.rec_model)
    if (self.reconcile):
      command.append("--reconcile")
    else:
      command.append("--do-not-reconcile")
    if (self.reconcile_samples > 0):
      command.append("--reconciliation-samples")
      command.append(str(self.reconcile_samples))
    logs_file_path = os.path.join(self.test_dir, "tests_logs.txt")
    command.append("-p")
    command.append(os.path.join(self.test_dir, "generax"))
    try:
      with open(logs_file_path, "w") as writer:
        subprocess.check_call(command, stdout = writer, stderr = writer)
    except:
      return
    self.ok = True

  def check(self):
    return self.ok
  
  def destroy_dir(self):
    try:
      shutil.rmtree(self.test_dir)
    except:
      pass

print("Running " + str(runs_number) + " randomized tests")
all_ok = True
mkdir(OUTPUT)

for i in range(0, runs_number):
  seed = i + 1
  runner = Runner(seed)
  runner.run()
  if (not runner.check()):
    all_ok = False
    print("Randomized test with seed " + str(seed) + " failed. Outputdir: " + runner.test_dir)
  else:
    print("Test with seed " + str(seed) + " ok")
    #runner.destroy_dir()

if (all_ok):
  print("All randomized tests successfully finished")
else:
  print("Some of the randomized tests failed")
  sys.exit(1)


