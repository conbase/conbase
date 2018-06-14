############### Stats #################
stats_params = {
# Read
"fragment_length" : 200,
"mapping_quality" : 20,
"base_quality" : 20,
# ALT
"dp_limit" : 10,
"alt_ratio_limit" : 0.2,
"sample_vote_limit" : 2,  
"vote_ratio_limit" : 1,
"snp_read_limit" : 1,
# Indel
"indel_ratio" : 0.01,
# Bulk
"bulk_ref_limit" : 1,  #this requires the remaining percent to be only a1 
}
############### Analyze #################
analyze_params = {
"dp_ms_limit" : 2,
"snp_total_vote" : 0.9,
"snp_vote_ratio" : 0.9,

"msp_ratio" : 0.9,
"msp_internal_ratio" : 0.1,
"msp_c2_external_error_ratio" : 0.1, # 1-msp_ratio 

"c3_a1_limit" : 2,
"c3_homo_limit" : 2,

"homo_error_allowed" : 0,

"ms_group_ratio" : 0.01,
"win_internal_group_ratio" : 0.1,

"sample_ms_vote_limit" : 2,
"vote_ms_ratio_limit" : 0.9,

"conflicting_upper_limit" : 0,
"c3_conflicting_upper_limit" : 0,
"a1_lower_limit" : 2,
"bulk_dp_interval" : (15,55), # Must have the following format (min, max)
}

################ Misc ##################
misc_params = {
"mut_nr_limit" : 1,
"mut_dist_limit" : 1000,

"snp_nr_limit" : 10,
"snp_dist_limit" : 1000,
}
trees = {
    'tree1':{
        'samples':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'],
        'params': { 'MIN_HET': 2, 'MIN_HOM': 1, 'MAX_HET': 100, 'MAX_HOM': 100}
    },
    'all':{
        'samples':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'],
        'params': { 'MIN_HET': 2, 'MIN_HOM': 1, 'MAX_HET': 100, 'MAX_HOM': 100}
    }
}

