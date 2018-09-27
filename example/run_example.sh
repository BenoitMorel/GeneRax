dataset=simuls_998
gene_tree=${dataset}/raxmlGeneTree.newick
alignment=${dataset}/alignment.msa
species_tree=${dataset}/speciesTree.newick
strategy=NNI
threads=4
prefix=output


../build/bin/JointTreeSearch -g ${gene_tree} -a ${alignment} -s ${species_tree} --strategy ${strategy} -t ${threads} -p ${prefix}
