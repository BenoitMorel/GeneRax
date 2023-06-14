![alt text](doc/generax_logo.png?raw=true "GeneRax")

News: 
* 03.2021: GeneRax 2 is released, allowing rooted species tree inference from gene trees (I am still working on updating the bioconda package) 
* 11.2020: we just created a new GeneRax google groups for asking questions and reporting issues:  https://groups.google.com/g/generaxusers

# GeneRax 

GeneRax is a parallel tool for species tree-aware maximum likelihood based gene family tree inference under gene duplication, transfer, and loss.

It infers gene family trees from their aligned sequences, the mapping between genes and species, and a rooted undated species tree. In addition, it infers the duplication, transfer and loss events that best (in terms of maximum likelihood) reconcile the gene family trees with the species trees.

It accounts for sequence substitutions, gene duplication, gene loss and horizontal gene transfer.

When using GeneRax, please cite: [https://academic.oup.com/mbe/article/doi/10.1093/molbev/msaa141/5851843](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msaa141/5851843)

GeneRax is also available on [`bioconda`](https://anaconda.org/bioconda/generax) 

# SpeciesRax

SpeciesRax is part of the GeneRax tool and is available since GeneRax v2.0.0. SpeciesRax infers a rooted species tree from a set of unrooted gene trees. When using SpeciesRax, please cite [https://academic.oup.com/mbe/article/39/2/msab365/6503503](https://academic.oup.com/mbe/article/39/2/msab365/6503503)

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

## Issues and questions

For questions, issues or feedback, please post on our google group: https://groups.google.com/g/generaxusers
When reporting an issue, please send us at least the command line you ran, the logs file and the families file. The more information we get, the quicker we can solve the problems :-)

