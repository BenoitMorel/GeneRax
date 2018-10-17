dataset=simuls_998
gene_tree=${dataset}/raxmlGeneTree.newick
alignment=${dataset}/alignment.msa
species_tree=${dataset}/speciesTree.newick
threads=4
prefix=output

echo "Running NNI example"
strategy=NNI
../build/bin/JointSearch -g ${gene_tree} -a ${alignment} -s ${species_tree} --strategy ${strategy} -t ${threads} -p ${prefix}_${strategy} --verbose

echo "Running SPR example"
strategy=SPR
../build/bin/JointSearch -g ${gene_tree} -a ${alignment} -s ${species_tree} --strategy ${strategy} -t ${threads} -p ${prefix}_${strategy} --verbose

