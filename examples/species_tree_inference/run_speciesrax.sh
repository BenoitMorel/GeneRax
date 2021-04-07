generax=build/bin/generax
example=examples/species_tree_inference
output=$example/output
families=$example/families_speciesrax.txt
speciestree=MiniNJ
cores=2
sprradius=2

if [ ! -f "$generax" ]; then
  echo "Can't find executable $generax. Please check that you compiled generax and that you ran this script from the repository root directory."
  exit 1
fi

rm -rf $output

mpiexec -np $cores $generax --families $families --species-tree $speciestree --strategy SKIP --rec-model UndatedDTL --per-family-rates --prune-species-tree --si-estimate-bl --si-quartet-support --prefix $output --si-strategy HYBRID

