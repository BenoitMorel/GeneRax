![alt text](doc/generax_logo.png?raw=true "GeneRax")


# GeneRax 

GeneRax is a parallel tool for species tree-aware maximum likelihood based gene family tree inference under gene duplication, transfer, and loss.

It infers gene family trees from their aligned sequences, the mapping between genes and species, and a rooted undated species tree. In addition, it infers the duplication, transfer and loss events that best (in terms of maximum likelihood) reconcile the gene family trees with the species trees.

It accounts for sequence substitutions, gene duplication, gene loss and horizontal gene transfer.

When using GeneRax, please cite: [https://academic.oup.com/mbe/article/doi/10.1093/molbev/msaa141/5851843](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msaa141/5851843)

GeneRax is also available on [`bioconda`](https://anaconda.org/bioconda/generax) 

# SpeciesRax

SpeciesRax is part of the GeneRax tool and is available since GeneRax v2.0.0. SpeciesRax infers a rooted species tree from a set of unrooted gene trees. When using SpeciesRax, please cite [https://academic.oup.com/mbe/article/39/2/msab365/6503503](https://academic.oup.com/mbe/article/39/2/msab365/6503503)

Update: we have released a new tool [AleRax](https://github.com/BenoitMorel/AleRax) ([paper](https://academic.oup.com/bioinformatics/article/40/4/btae162/7633408)) which an independent study ([preprint](https://www.biorxiv.org/content/10.1101/2024.11.20.624597v4.full.pdf))  has shown to be more accurate than SpeciesRax and other competing methods for species tree inference.

## Requirement

(If you are not installing with bioconda)

* A Linux or MacOS environnement
* gcc 5.0 or > 
* CMake 3.6 or >
* MPI

## Installation 

(Please note that you can also install through [`bioconda`](https://anaconda.org/bioconda/generax))


 To download GeneRax, please use git,  and clone with --recursive!!!

```
git clone --recursive https://github.com/BenoitMorel/GeneRax
```

To build the sources:
```
./install.sh
```
## Running

See the wiki (https://github.com/BenoitMorel/GeneRax/wiki)

### Running SpeciesRax with example data

Download and extract the test datasets:

```
wget https://cme.h-its.org/exelixis/material/generax_data.tar.gz
tar xzf generax_data.tar.gz
```

Prepare a families file (e.g. `families.txt`) listing gene trees and mappings:

```
[FAMILIES]
- HBG284008_sim
starting_gene_tree = generax_data/cyano_simulated/families/HBG284008_sim/gene_trees/raxml-ng.LG+G+I.geneTree.newick
mapping = generax_data/cyano_simulated/families/HBG284008_sim/mappings/mapping.link
- HBG284199_sim
starting_gene_tree = generax_data/cyano_simulated/families/HBG284199_sim/gene_trees/raxml-ng.LG+G+I.geneTree.newick
mapping = generax_data/cyano_simulated/families/HBG284199_sim/mappings/mapping.link
```

Run SpeciesRax to infer a species tree:

```
./build/bin/generax --families families.txt --species-tree MiniNJ --strategy SKIP --rec-model UndatedDTL --per-family-rates --prune-species-tree --si-strategy HYBRID --prefix output
```

### Checkpointing

For long-running SpeciesRax analyses, enable checkpointing with the `--checkpoint` flag. This saves progress after each search iteration so the run can be resumed if interrupted:

```
./build/bin/generax --families families.txt --species-tree MiniNJ --strategy SKIP --rec-model UndatedDTL --per-family-rates --prune-species-tree --si-strategy HYBRID --prefix output --checkpoint
```

If the run is interrupted, simply re-run the same command. It will detect the checkpoint in the output directory and resume from where it left off.

A checkpoint test script is provided in `tests/test_speciesrax_checkpoint.sh`. It requires the test data above and should be run from the repository root.

## Issues and questions

For questions, issues or feedback, please post on our google group: https://groups.google.com/g/generaxusers
When reporting an issue, please send us at least the command line you ran, the logs file and the families file. The more information we get, the quicker we can solve the problems :-)

