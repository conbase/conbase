![conbase](https://raw.githubusercontent.com/conbase/conbase/master/conbase_logo.png)

## Conbase: a software for unsupervised discovery of clonal somatic mutations in single cells through read phasing
 
Conbase is a software for identification of somatic mutations in single cell DNA sequencing data exhibiting high rates of allelic dropout and at low read depth. Conbase leverages data from multiple samples in a dataset, and utilizes read phasing to call somatic single nucleotide variants and to accurately predict genotypes in whole genome amplified single cells.

Algorithm developed by [Ezeddin Al Hakim](https://github.com/ezeddin), [Marie Kindblom](https://github.com/mkindblom) & [Joanna Hård](https://github.com/joannahard).

## How to run

Stats.py:   

    python3 bin/Main.py stats <snp path> <bam path> <reference path> <number of threads <output name>


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

## Implementation by

* **Marie Kindblom** - (https://github.com/mkindblom)
* **Ezeddin Al Hakim** - (https://github.com/ezeddin)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
