generax=build/bin/generax
example=examples/compare_two_species_trees
output1=$example/output1
output2=$example/output2
families=$example/families_joint_likelihood.txt
speciestree1=data/plants/species_trees/speciesTree.newick
speciestree2=data/plants/species_trees/otherSpeciesTree.newick
cores=2
sprradius=2

if [ ! -f "$generax" ]; then
  echo "Can't find executable $generax. Please check that you compiled generax and that you ran this script from the repository root directory."
  exit 1
fi

rm -rf $output1
rm -rf $output2

mpiexec -np $cores $generax --families $families --species-tree $speciestree1 --rec-model UndatedDL --per-family-rates  --prefix $output1 --max-spr-radius $sprradius
mpiexec -np $cores $generax --families $families --species-tree $speciestree2 --rec-model UndatedDL --per-family-rates  --prefix $output2 --max-spr-radius $sprradius

