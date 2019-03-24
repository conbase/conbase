![conbase](https://raw.githubusercontent.com/conbase/conbase/master/conbase_logo.png)

## Conbase: a software for unsupervised discovery of clonal somatic mutations in single cells through read phasing
 
Conbase is a software for unsupervised discovery of clonal somatic mutations in single DNA sequencing data from whole genome amplified single cells. Conbase leverages phased read data from multiple samples in a dataset to achieve increased confidence in somatic variant calls and genotype predictions.

Algorithm developed by [Ezeddin Al Hakim](https://github.com/ezeddin), [Marie Kindblom](https://github.com/mkindblom) & [Joanna Hård](https://github.com/joannahard).


## Preprocessing

Conbase takes three inputs: bams files from single cells and an unamplified bulk sample, heterozygous germline SNV coordinates, a human reference genome in fasta format. A pipeline for read processesing, bulk variant calling and bulk variant filtering is available *here*

It is recommended that heterozygous germline SNVs called in the unamplified bulk sample are conservatively filtered, since false positive SNVs result in higher false discovery rate in somatic variant calling with Conbase. To reduce the number of false positive SNVs it is recommended to run Misc.py, included in the Conbase package.

Misc.py:
	python3 Misc.py --duplicate_regions <prefiltered snp file> <bulk bam path> <reference path> <number of threads> <filtered snp file>


Before running Misc.py, adjust parameter settings in Params.py

Params.py:
	# Defining variant supporting sites in the unamplified bulk sample
        # Bulk
	"bulk_ref_limit" : 0.9,  #the fraction of reads supporting the reference base in the unamplified bulk sample. The remaining fraction supports a non-reference base

	# The number of sites supporting a variant allowed per genomic window in the unamplified bulk sample
        "snp_nr_limit" : 10, #the number of variant supporting sites
        "snp_dist_limit" : 1000, #the size of the genomic window (bp)

Input format prefiltered snp file
Each row correspond to the coordinates for filtered germline heterozygous SNVs, the reference base (REF) and the alternative base (ALT) observed in the unamplified bulk sample

    CHROM       POS     REF     ALT
    1   23142496        C       A
    .
    .
    .



## How to run Conbase

Stats.py:   

    python3 bin/Main.py stats <filtered snp file path> <bam file path> <reference path> <number of threads> <output>


Analyze.py: 

    python3 bin/Main.py analyze <output.json path> <output>
    

Input format bam file
Each row correspond to the sample name and path to a bam file

    NAME	PATH
    BULK	data/bulk.bam
    cell_1	data/cell_1.bam
    cell_2	data/cell_2.bam
    .
    .


Input format filtered snp file
Each row correspond to the coordinates for filtered germline heterozygous SNVs and the reference (REF) and alternative (ALT) bases observed in the unamplified bulk sample

    CHROM	POS	REF	ALT
    1	23142496	C	A
    .
    .
    .


Before running Stats.py and Analyze.py, adjust parameter settings in Params.py
Conbase variant calling is based on assumptions associated with expected observations in true clonal somatic variant sites. These assumptions are coded as adjustable parameters defined in the Params.py file.



	Stats #################
	stats_params = {
	# Read
	"fragment_length" : 650, #insert size of sequencing library
	"mapping_quality" : 20, #required mapping quality of analyzed reads
	"base_quality" : 20, #required base quality of analyzed bases
	# ALT
	"dp_limit" : 10, #the number of reads required to cover a putative variant site
	"alt_ratio_limit" : 0.2, #the required variant allele fraction in putative variant sites
	"sample_vote_limit" : 2, #the number of samples required to support a putative variant
	"vote_ratio_limit" : 0.9, #the fraction of samples required to support the same base in a putative variant site
	"snp_read_limit" : 1, #what?
	# Indel
	"indel_ratio" : 0.05, #the fraction of reads (across the dataset) allowed to contain an indel in a putative variant site
	# Bulk
	"bulk_ref_limit" : 1,  #the fraction of reads required to support the reference base in the unamplified bulk sample
	"acceptable_bases" : ['A', 'C', 'G', 'T'],
	}
	############### Analyze #################
	analyze_params = {
	
	# Genotype prediction in samples with sufficient data amount
	"dp_tuple_limit" : 10, #the number of phased reads required for predicting a genotype
	"snp_total_vote" : 0.9, #the fraction of available gSNVs required to contribute to genotype prediction
	"snp_vote_ratio" : 0.9, #the fraction of available gSNVs required to agree on the same genotype prediction
        

	"tuples_ratio" : 0.9, #the minimal ratio required of the total that is made out of a given tuple pair
	"tuples_internal_ratio" : 0.1, #the minimal internal ratio required for a tuple pair
	"tuples_c2_external_error_ratio" : 0.1, # 1-tuples_ratio

	"tuple_group_ratio" : 0.01, #the max "error" allowed between two groups of tuplepairs for a sample to vote for a winning tuplepair
	"win_internal_group_ratio" : 0.1, #minimal internal ratio for a winning tuplepair in a voting sample
	"sample_tuple_vote_limit" : 2, #the minimum number of samples that needs to vote for a winning tuple pair
	"vote_tuple_ratio_limit" : 0.9, #the minimum percentage of votes that the winning tuple pair needs to acquire


	# Inferring genotypes in samples with low read depth
        "a1_lower_limit" : 2, #minimum depth of alternative alleles (site olsy passes if conflicting_upper_limit and c3_conflicting_upper_limit also holds)
        "c3_homo_limit" : 2, #min depth of tuple type "RA" required when the winning tuplepair is {AA, RR} for the site to be regarded unmutated
        "c3_a1_limit" : 2, #min depth for an alternative base in the case of no support of a mutated site


	# Filtering variant sites
	"conflicting_upper_limit" : 0, #the number of samples with predicted genotypes allowed to display support for discordant alleles
	"c3_conflicting_upper_limit" : 0, #the number of samples with inferred genotypes allowed to display support for discordant alleles
	"homo_error_allowed" : 0, #the number of reads allowed to support a non-reference base in a sample predicted to be unmutated
	"bulk_dp_interval" : (15,65), #the minimal and maximal read depth required in the unamplified bulk sample


	}

	################ Misc ##################
	misc_params = {
	
	# The number of somatic variant calls expected per genomic window
	"mut_nr_limit" : 1, #the number of somatic variant calls
	"mut_dist_limit" : 1000, #the size of the genomic window (bp)



	# Defining output to display
	# MIN_HET : the minimum number of samples required to be predicted mutated
	# MIN_HOM : the minimum number of samples required to be predicted unmutated
	# MAX_HET : the maximum number of samples required to be predicted mutated
        # MAX_HOM : the maximum number of samples required to be predicted unmutated

	}
	trees = {
	    'tree1':{
	        'samples':['cell_1','cell_2','cell_3','cell_4','cell_5','cell_6','cell_7','cell_8','cell_9','cell_10'],
	        'params': { 'MIN_HET': 2, 'MIN_HOM': 1, 'MAX_HET': 9, 'MAX_HOM': 8}
	    },
	    'all':{
	        'samples':['cell_1','cell_2','cell_3','cell_4','cell_5','cell_6','cell_7','cell_8','cell_9','cell_10'],
	        'params': { 'MIN_HET': 2, 'MIN_HOM': 1, 'MAX_HET': 9, 'MAX_HOM': 8}
	    }
	}



	# Conbase can also be run in a supervised mode, to parse variant calls enriched in a defined population


	}
	trees = {
	    'tree1':{
	        'samples':['cell_1','cell_2','cell_3','cell_4','cell_5'],
	        'params': { 'MIN_HET': 2, 'MIN_HOM': 1, 'MAX_HET': 4, 'MAX_HOM': 3}
	    },
	    'tree2':{
	        'samples':['cell_6','cell_7','cell_8','cell_9','cell_10'],
	        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 0, 'MAX_HOM': 5}
	    },
	    'all':{
	        'samples':['cell_1','cell_2','cell_3','cell_4','cell_5','cell_6','cell_7','cell_8','cell_9','cell_10'],
	        'params': { 'MIN_HET': 2, 'MIN_HOM': 1, 'MAX_HET': 4, 'MAX_HOM': 8}
	    }
	}




## Implementation by
Marie Kindblom (https://github.com/mkindblom)
Ezeddin Al Hakim (https://github.com/ezeddin)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
