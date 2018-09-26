dataset=simuls_998
gene_tree=${dataset}/raxmlGeneTree.newick
alignment=${dataset}/alignment.msa
species_tree=${dataset}/speciesTree.newick
strategy=NNI
threads=4
output=output


../build/bin/JointTreeSearch ${gene_tree} ${alignment} ${species_tree} ${strategy} ${threads} ${output}
