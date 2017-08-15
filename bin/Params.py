############### Stats #################
stats_params = {
# Read
"fragment_length" : 600,
"mapping_quality" : 20,
"base_quality" : 20,
# ALT
"dp_limit" : 15,
"alt_ratio_limit" : 0.3,
"sample_vote_limit" : 2,  
"vote_ratio_limit" : 0.7,
"snp_read_limit" : 1,
# Indel
"indel_ratio" : 0.05,
# Bulk
"bulk_ref_limit" : 0.91,  #this requires the remaining percent to be only a1 
}
############### Analyze #################
analyze_params = {
"dp_ms_limit" : 8,
"snp_total_vote" : 0.9,
"snp_vote_ratio" : 0.9,

"msp_ratio" : 0.9,
"msp_internal_ratio" : 0.1,
"msp_c2_external_error_ratio" : 0.1, # 1-msp_ratio 

"c3_a1_limit" : 3,
"c3_homo_limit" : 3,

"homo_error_allowed" : 0,

"ms_group_ratio" : 0.01,
"win_internal_group_ratio" : 0.1,

"sample_ms_vote_limit" : 3,
"vote_ms_ratio_limit" : 0.9,

"conflicting_upper_limit" : 0,
"c3_conflicting_upper_limit" : 0,
"a1_lower_limit" : 2,
"bulk_dp_interval" : (15,45), # Must have the following format (min, max)
}

################ Misc ##################
misc_params = {
"mut_nr_limit" : 1000,
"mut_dist_limit" : 1000,

"snp_nr_limit" : 50,
"snp_dist_limit" : 1000,
}

trees = {
    'tree1':{
        'samples':['fibroblast_1','fibroblast_2','fibroblast_3','fibroblast_5','fibroblast_6'],
        'params': { 'HET': 0, 'HOM': 1}
    },
    'tree2':{
        'samples':['fibroblast_22','fibroblast_24','fibroblast_27','fibroblast_30','fibroblast_33','fibroblast_34','fibroblast_36','fibroblast_38','fibroblast_39','fibroblast_40','fibroblast_41'],
        'params': { 'HET': 0, 'HOM': 1}
    },
    'all':{
        'samples':['fibroblast_1','fibroblast_2','fibroblast_3','fibroblast_5','fibroblast_6','fibroblast_22','fibroblast_24','fibroblast_27','fibroblast_30','fibroblast_33','fibroblast_34','fibroblast_36','fibroblast_38','fibroblast_39','fibroblast_40','fibroblast_41'],
        'params': { 'HET': 0, 'HOM': 2}
    }
}