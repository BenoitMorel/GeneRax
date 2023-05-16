#!/usr/bin/env python
import os
import sys
import shutil
import subprocess


def extract(llpath, treespath, conseldir, test, cutoff, outputdir):
  try:
    os.mkdir(outputdir)
  except:
    print("Could not create directory " + outputdir + ".")
    print("Please check that this directory does not already exists")
    sys.exit(1)
  llcopy = os.path.join(outputdir, "input_likelihoods.txt")
  run_name = os.path.join(outputdir, "input_likelihoods")
  treescopy = os.path.join(outputdir, "input_trees.txt")
  shutil.copy(llpath, llcopy)
  shutil.copy(treespath, treescopy)
  makermt = os.path.join(conseldir, "makermt")
  consel = os.path.join(conseldir, "consel")
  catpv = os.path.join(conseldir, "catpv")
  cmd1 = [makermt, "--puzzle", llcopy]
  cmd2 = [consel, run_name]
  cmd3 = [catpv, run_name]
  subprocess.check_call(cmd1)
  subprocess.check_call(cmd2)
  res = subprocess.check_output(cmd3)
  print(res)
  lines = res.split("\n")
  labels = lines[2].split()
  lines = lines[3:]
  item_index = labels.index("item")
  metric_index = labels.index(test)
  trees = open(treespath).read().split("\n")

  accepted_trees = os.path.join(outputdir, "accepted_trees.txt")
  accepted_count = 0
  with open(accepted_trees, "w") as writer:
    for line in lines:
      values = line.split()
      if (len(values) == 0):
        continue
      metric = float(values[metric_index]) 
      if (metric < 0.05):
        continue
      item = int(values[item_index])
      tree = trees[item - 1]
      writer.write(str(metric) + "\n")
      writer.write(tree + "\n")
      accepted_count += 1
  print(str(accepted_count) + " trees passed the " + test + "test with cutoff " + str(cutoff) + ". See output in " + accepted_trees)

if (__name__ == "__main__"):
  if (len(sys.argv) != 7):
    print("Syntax python " + os.path.basename(__file__) + " llpath treespath conselbindir test cutoff outputdir")
    print("  llpath: path to file containing per tree and per-family likelihoods")
    print("  treespath: path to file containing the list of newick trees")
    print("  conselbindir: path to directory containing consel binaries")
    print("  test: statistical test to use. Possible values: au np bp pp kh sh wkh wsh")
    print("  cutoff: keep only trees with pvalues above this cutoff (typically 0.05)")
    print("  outputdir: output directory. Will be created")

    sys.exit(1)
  llpath = sys.argv[1]
  treespath = sys.argv[2]
  conseldir = sys.argv[3]
  test = sys.argv[4]
  cutoff = float(sys.argv[5])
  outputdir = sys.argv[6]
  extract(llpath, treespath, conseldir, test, cutoff, outputdir)
  








