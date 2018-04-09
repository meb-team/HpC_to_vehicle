# HpC_to_vehicle
[This project is in development. Documentation is fairly light. You are welcomed to use this software, but please expect it to change in non-trivial ways. Feel free to return your impressions.]

## Contents 

* [Introduction](https://github.com/meb-team/HpC_to_vehicle/blob/master/README.md#introduction)
* [Installation](https://github.com/meb-team/HpC_to_vehicle/blob/master/README.md#installation)
* [How to use ?](https://github.com/meb-team/HpC_to_vehicle/blob/master/README.md#how-to-use)
* [Bugs](https://github.com/meb-team/HpC_to_vehicle/blob/master/README.md#bugs)
* [Citation](https://github.com/meb-team/HpC_to_vehicle/blob/master/README.md#citation)

## Introduction

The perl script HpC_to_vehicle retrieves the results of the perl script AccNetPhylA (files to visualize the networks, multiple alignment of Muscle, phylogenetic results of PHYLIP, ...) to create HpC (Homologous protein Cluster) and GU (Genomic Unit), to thus make different networks representing the links between "vehicles" of genes or proteins.

For each vehicle, there are two types of graphs:
		-network between vehicles that can belong to the same GU
		-network between different GU vehicles

For gene vehicle analysis, a patristic distance (calculated from tree branch lengths describe the amount of genetic change represented by a tree) will be used to build two more graphs.

In the end, there will be 4 networks resulting from an analysis with the genes and 2 from an analysis with the proteins.

## Installation 

### Prerequisites

 * [Perl](https://www.perl.org/) - scripting language
 * [Python 2.7](https://www.python.org/download/releases/2.7/) or [Python 3+](https://www.python.org/download/releases/3.0/) - scripting language
 * [MUSCLE](https://www.drive5.com/muscle/) - Tool for MUltiple Sequence Comparison by Log-Expectation
 * [trimAl](http://trimal.cgenomics.org/) - Tool for removal of spurious sequences or poorly aligned regions
 * [PHYLIP](http://evolution.genetics.washington.edu/phylip.html) - PHYlogeny Inference Programs
 * [PhyML](http://www.atgc-montpellier.fr/phyml/) - Phylogeny software

### GitHub

Choose somewhere to put it, for example in your home directory (no root access required):

```bash
cd $HOME
```

Clone the latest version of the repository:

```bash
git clone https://github.com/cleasiguret/HpC_to_vehicle.git
```

## How to use ?

```
    Usage : 
      Vehicle analysis with genes :
	      perl HpC_to_vehicle.pl -equi equivalence.txt -netin Network.csv -tabin Table.csv -seq nt 
	      -ppA 'Analyse/muscle/*' -pg 'directoryGene/*'
	
	    Vehicle analysis with proteins :
	      perl HpC_to_vehicle.pl -equi equivalence.txt -netin Network.csv -tabin Table.csv -seq aa 
	      -ppM 'Analyse/phylip/'


    Global options:

    -Help|help|h, produces this help file. 
	
    -Verbose|verbose|v, Verbose mode. Can be used multiple times for increased verbosity.
	
    -equi, Equivalence file of AccNetPhylA.pl.
	
    -netin, Network filename of AccNetPhylA.pl.
	
    -tabin, Table filename of AccNetPhylA.pl.

    -seq, Proteic or nucleic sequence (nt|aa).
	
    -ppA, Path of the protein alignment files to analyze (with -seq nt).
	
    -ppM, Path of the protein matrix files to analyze (with -seq aa).
	
    -pg, Path of the fasta files of genes to to analyze (with -seq nt).
	
    -clean, Remove and class files. Default(yes).

    -dir, Directory's name. Default(Analyse_vehicle).
  
```


## Bugs

* Submit problems or requests here: https://github.com/meb-team/HpC_to_vehicle/issues


## Citation

### Authors
* Siguret Clea - [cleasiguret](https://github.com/cleasiguret)
