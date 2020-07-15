# MuteBaSS
MuteBaSS (**Mu**l**t**i-sp**e**cies **Ba**lancing **S**election **S**ummary statistics) is a set of Python scripts to perform scans with summary statitics HKA, NCD, NCD<sub>opt</sub>, NCD<sub>sub</sub>, and NCD<sub>mid</sub>, for detecting footprints of long-term balancing selection affecting one or more species (Cheng & DeGiorgio, 2019). Operation of this package requires a UNIX environment with Python 2.7 and above.

Please cite us as:

> [X. Cheng and M. DeGiorgio. Detection of shared balancing selection in the absence of trans-species polymorphism. 2019. *Mol. Biol. Evo.* 36(1):177--199](https://academic.oup.com/mbe/article/36/1/177/5150441)

## Installation
 We distribute MuteBaSS in compressed (tar.gz) format. In addition to the `MuteBaSS.py` script, we also included the user manual and example data. The scripts included are designed to perform on a UNIX system. To unpack `MuteBaSS` from the command line, go to the directory where it is stored, and enter

    tar -xzvf MuteBaSS.tar.gz
    cd MuteBaSS/

## Quick Guide
### Performing scans
To run <code>MuteBaSS.py</code>, first make sure other python scripts are in the same directory, then format the dash commands and their arguments in a single line as

	python MuteBaSS.py -i <input file>  -c <p,x> [--check] [--tree <tree>] [--fixSize] 
        -w <window size> -s <step size> -o <output file> 
        [--NCD] [--tf <tf>] [--NCDopt] [--NCDsub] [--NCDmid]
        [--HKA] [--config <config file>] [--getConfig]

Alternatively, running command <code>python MuteBaSS.py</code> or <code>python MuteBaSS.py -h</code> would return the help page with detailed instructions.

### Input format
#### Input file
The input file should be *tab-delimited*, and must at least include, in addition to the physical position of each informative site, the number of ancestral alleles (denoted as *x*) and the total number of alleles sampled (denoted as *n*). To be an informative site, a site should either be polymorphic in only one of the species examined, or be monomorphically different across all species, with the pattern agreeing with the species tree. Sites that fit neither of these two types should be discarded. All sites should be bi-allelic. All input files should include one-line headers, otherwise the first line will be automatically excluded from analyses. To examine *K* species, the input file should at least contain the following columns:
    
| position | x1 | n1 | ... | xk | nk |
|----------|----|----|-----|----|----|
|  |  |   |  |  |  |

where **position** is for physical positions, and **xj** and **nj** denote the ancestral and total allele counts, respectively, at this position in species *j*, *j* = 1, 2, .., *K*. For this particular input above, the `-c <p,x>` argument should be `-c 1,2`. 
  
Here's another example input:
| Chr | pos | anc | drv | x1 | n1 | ... | xk | nk |
|-----|-----|-----|-----|----|----|-----|----|----|
|  |  |  |  |  |  |  | | |

for which the <code>-c <p,x></code> argument should be `-c 2,5`. 
  
  
#### Configuration file for HKA
When choosing to perform scans with the HKA statistic, users must provide a corresponding configuration file with `--config` argument. This file records the proportions (conditional on informative sites) of within-species polymorphisms and cross-species substitutions for each set of sample sizes. Configuration files should **not** have headers, should be *tab-delimited*, and each line should present the needed information in the following order in one line:

| \<n1\> | \<n2\> | ... | \<nk\> | \<poly\> | \<sub\> |
|----------|----|----|-----|----|----|

You can find an example file for 4-species configuration file in the folder <code>test/</code>, which reads:

| 50 | 50 | 50 | 50 | 0.283164 | 0.716836 |
|----|----|----|----|----------|----------|

