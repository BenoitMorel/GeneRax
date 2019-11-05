# GeneRax 

GeneRax is a parallel tool for species tree-aware maximum likelihood based gene tree inference under gene duplication, transfer, and loss.

It infers gene trees from their aligned sequences, the mapping between genes and species, and a rooted undated species tree.

It accounts for sequence substitutions, gene duplication, gene loss and horizontal gene transfer.

Preprint: https://www.biorxiv.org/content/biorxiv/early/2019/09/26/779066.full.pdf

GeneRax is also available on [`bioconda`](https://anaconda.org/bioconda/generax) 

## Requirement

(If you are not installing with bioconda)

* A Linux or MacOS environnement
* gcc 5.0 or > 
* CMake 3.6 or >
* MPI
* bison and flex parsers


## Installation 

(Please note that you can also install through [`bioconda`](https://anaconda.org/bioconda/generax))

Installing the dependencies on Ubuntu and other Debian-based systems:
```
sudo apt-get install flex bison libgmp3-dev
```

On other systems: [`GNU Bison`](http://www.gnu.org/software/bison/) [`Flex`](http://flex.sourceforge.net/) [`GMP`](https://gmplib.org/)


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

