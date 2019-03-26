![conbase](https://raw.githubusercontent.com/conbase/conbase/master/conbase_logo.png)

## Conbase: a software for unsupervised discovery of clonal somatic mutations in single cells through read phasing

Conbase is a software for unsupervised discovery of clonal somatic mutations in single DNA sequencing data from whole genome amplified single cells. Conbase leverages phased read data from multiple samples in a dataset to achieve increased confidence in somatic variant calls and genotype predictions.
https://doi.org/10.1186/s13059-019-1673-8

Algorithm developed by [Ezeddin Al Hakim](https://github.com/ezeddin), [Marie Kindblom](https://github.com/mkindblom) & [Joanna Hård](https://github.com/joannahard).


## Optional Preprocessing

Misc.py

    python3 Misc.py --duplicate_regions <prefiltered germline SNV file> <bulk bam> <reference> <number of threads> <output>

In Params.py, adjust:



	"bulk_ref_limit" : 0.9,  #Definition of a variant site, where the remaining fraction supports a non-reference base
    "snp_nr_limit" : 10, #The maximum number of variant sites per genomic window
    "snp_dist_limit" : 1000, #Genomic window size (bp)

Input format prefiltered germline SNV file

    CHROM       POS     REF     ALT
    1   23142496        C       A
    .
    .



## How to run Conbase

Stats.py:   

    python3 bin/Main.py stats <germline SNV file> <bam path file> <reference> <number of threads> <output>


Analyze.py:

    python3 bin/Main.py analyze <output.json> <output>


Input format bam path file


    NAME	PATH
    BULK	data/bulk.bam
    cell_1	data/cell_1.bam
    cell_2	data/cell_2.bam
    .
    .


Input format germline SNV file


    CHROM	POS	REF	ALT
    1	23142496	C	A
    .
    .


In Params.py, optionally adjust:

	################# Stats #################
	stats_params = {  

	"fragment_length" : 650, #Insert size of sequencing library  

	"mapping_quality" : 20, #Required mapping quality of analyzed reads  

	"base_quality" : 20, #Required base quality of analyzed bases

	"dp_limit" : 5, #Minimum read depth  

	"alt_ratio_limit" : 0.2, #Minimum variant allele fraction  

	"sample_vote_limit" : 2, #Minimum number of samples required to support a variant

	"vote_ratio_limit" : 0.9, #Minimum fraction of samples required to support the same base

	"snp_read_limit" : 1, #Minimun number of reads covering an SNP for a given sample to regard its reads (and voting) for the SNP

	"indel_ratio" : 0.05, #The fraction of reads in the entire dataset allowed to contain an indel  

	"bulk_ref_limit" : 1,  #The fraction of reads required to support the reference base in the unamplified bulk sample  

	"acceptable_bases" : ['A', 'C', 'G', 'T'],
	}  




    ############### Analyze #################  
	analyze_params = {

	# Genotyping is based on phased read concordance. saved as tuples containing base observations in somatic variant sites and germline SNV sites present in the same read or read pair
   
    # Please note the following annotations: a heterozygous site can be of type = {RR, AA} or {RA, AR} and a homozygous site is always {RR, RA} where R is the reference base and A the alternative. Any of these three combinations is a combination of two tuples, and will be referred to as a "tuple pair" below.

	"dp_tuple_limit" : 5, #Minimum number of phased reads  

	"snp_total_vote" : 0.9, #Fraction of phased germline SNVs required to contribute to genotype prediction  

	"snp_vote_ratio" : 0.9, #Fraction of phased germline SNVs required to support the same genotype prediction  

	"tuples_ratio" : 0.9, #The minimal ratio required of the total that is made out of a given tuple pair  

	"tuples_internal_ratio" : 0.1, #The minimal internal ratio required for a tuple pair  

	"tuples_c2_external_error_ratio" : 0.1, # 1-tuples_ratio  

    # Tuple
	"tuple_group_ratio" : 0.01, #The max "error" allowed between two groups of tuplepairs for a sample to vote for a winning tuplepair  

	"win_internal_group_ratio" : 0.1, #Minimal internal ratio for a winning tuplepair in a voting sample  

	"sample_tuple_vote_limit" : 2, #the minimum number of samples that needs to vote for a winning tuple pair  

	"vote_tuple_ratio_limit" : 0.9, #the minimum percentage of votes that the winning tuple pair needs to acquire


	# Genotype prediction in samples with low read depth
    "a1_lower_limit" : 2, #minimum depth of alternative alleles (site only passes if conflicting_upper_limit and c3_conflicting_upper_limit also holds)  

    "c3_homo_limit" : 2, #min depth of tuple type "RA" required when the winning tuplepair is {AA, RR} for the site to be regarded unmutated  

    "c3_a1_limit" : 2, #min depth for an alternative base in the case of no support of a mutated site


	# Filtering variant sites
	"conflicting_upper_limit" : 0, #The number of samples with predicted genotypes allowed to display support for discordant alleles  

	"c3_conflicting_upper_limit" : 0, #The number of samples with inferred genotypes allowed to display support for discordant alleles  

	"homo_error_allowed" : 0, #The number of reads allowed to support a non-reference base in a sample predicted to be unmutated  

	"bulk_dp_interval" : (15,65), #The minimal and maximal read depth required in the unamplified bulk sample in a predicted somatic variant site 


	}

	################ Misc ##################
	misc_params = {

	"mut_nr_limit" : 1, #Maximum number of somatic variants per genomic window
	"mut_dist_limit" : 1000, #Genomic window size (bp)



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
