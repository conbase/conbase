# conbase

## conbase algorithm description
 
conbase is a method for calling clonal somatic mutations in single cell data exhibiting high error and allelic dropout rates. Increased confidence in variants is achieved by requiring that the mutation is present on either the maternal or paternal allele. By taking advantage of read phasing, absence of a mutation can be determined despite high rates of allelic dropout. The somatic mutation site prediction is based on a set of parameters associated with observations in reads and read pairs covering heterozygous SNVs and sites containing an alternative allele present in a subset of samples. 

conbase consists of two separate subprograms, conbase_stats and conbase_analyze. conbase_stats iterates over positions upstream and downstream heterozygous SNVs in a genomic window, and prints statistics about sites where an alternative base is present. The genomic window size is determined by the maximal distance from a SNV covered by read pairs in any direction in any of the samples (maximum library fragment size). 

Algorithm developed by [Ezeddin Al Hakim](https://github.com/ezeddin), [Marie Kindblom](https://github.com/mkindblom) & [Joanna Hård](https://github.com/joannahard).

## How to run

Stats.py:   

    python3 bin/Main.py stats <snp path> <bam path> <reference path> <number of nodes> <output name>


Analyze.py: 

    python3 bin/Main.py analyze <json path> <output name>
    

input format bam file

    NAME	PATH
    cell1	data/cell1.bam
    .
    .
    .

input format snp file

    CHROM	POS	REF	ALT
    1	23142496	C	A
    .
    .
    .

## Authors

* **Marie Kindblom** - (https://github.com/mkindblom)
* **Ezeddin Al Hakim** - (https://github.com/ezeddin)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
