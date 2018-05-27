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
"dp_tuple_limit" : 10,
"snp_total_vote" : 0.9,
"snp_vote_ratio" : 0.9,

"tuples_ratio" : 0.9,
"tuples_internal_ratio" : 0.1,
"tuples_c2_external_error_ratio" : 0.1, # 1-tuples_ratio 

"c3_a1_limit" : 2,
"c3_homo_limit" : 2,

"homo_error_allowed" : 0,

"tuple_group_ratio" : 0.01,
"win_internal_group_ratio" : 0.1,

"sample_tuple_vote_limit" : 2,
"vote_tuple_ratio_limit" : 0.9,

"conflicting_upper_limit" : 100,
"c3_conflicting_upper_limit" : 100,
"a1_lower_limit" : 2,
"bulk_dp_interval" : (15,55), # Must have the following format (min, max)
}

################ Misc ##################
misc_params = {
"mut_nr_limit" : 100,
"mut_dist_limit" : 1000,

"snp_nr_limit" : 10,
"snp_dist_limit" : 1000,
}



#trees = {   
#    'tree1':{
#        'samples':["fib22", "fib24", "fib27", "fib30", "fib33", "fib34", "fib36", "fib38", "fib39", "fib40", "fib41"],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 100}
#    },
#    'tree2':{
#        'samples':["fib1", "fib2", "fib4", "fib5", "fib6"],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 100}
#    },
#    'all':{
#        'samples':["fib22", "fib24", "fib27", "fib30", "fib33", "fib34", "fib36", "fib38", "fib39", "fib40", "fib41", "fib1", "fib2", "fib4", "fib5", "fib6"],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 1, 'MAX_HET': 100, 'MAX_HOM': 100}
#    }
#}




#trees = {
#    'tree1':{
#        'samples':["fibroblast_22", "fibroblast_24", "fibroblast_27", "fibroblast_30", "fibroblast_33", "fibroblast_34", "fibroblast_36", "fibroblast_38", "fibroblast_39", "fibroblast_40", "fibroblast_41"],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 100}
#    },
#    'tree2':{
#        'samples':["fibroblast_1", "fibroblast_2", "fibroblast_3", "fibroblast_4", "fibroblast_5", "fibroblast_6"],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 100}
#    },
#    'all':{
#        'samples':["fibroblast_22", "fibroblast_24", "fibroblast_27", "fibroblast_30", "fibroblast_33", "fibroblast_34", "fibroblast_36", "fibroblast_38", "fibroblast_39", "fibroblast_40", "fibroblast_41", "fibroblast_1", "fibroblast_2", "fibroblast_3", "fibroblast_4", "fibroblast_5", "fibroblast_6"],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 1, 'MAX_HET': 100, 'MAX_HOM': 100}
#    }
#}











trees = {   
    'tree1':{
        'samples':['110','111','112','113','114','118','119','120','121'],
        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 100}
    },
    'tree2':{
        'samples':['134','135','136','137','138','142','143','144','145'],
        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 100}
    },
    'all':{
        'samples':['110','111','112','113','114','118','119','120','121','134','135','136','137','138','142','143','144','145'],
        'params': { 'MIN_HET': 1, 'MIN_HOM': 1, 'MAX_HET': 100, 'MAX_HOM': 100}
    }
}


#trees = {
#    'tree1':{
#        'samples':['malbac_4','malbac_4_long','malbac_5','malbac_6'],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 100}
#    },
#    'tree2':{
#        'samples':['sigma_1','sigma_1_long','sigma_2','sigma_3'],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 100}
#    },
#    'all':{
#        'samples':['malbac_4','malbac_4_long','malbac_5','malbac_6','sigma_1','sigma_1_long','sigma_2','sigma_3'],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 100}
#    }
#}





#trees = {   
#    'tree1':{
#        'samples':['Cardio_1_D7','Cardio_1_F7','Cardio_1_F9','Cardio_1_G2','Cardio_1_G9','Cardio_1_H4','Cardio_2_C2','Cardio_2_G4'],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 0, 'MAX_HOM': 100}
#    },
#    'tree2':{
#        'samples':['NEURONS','NeuN_2_D8','NeuN_2_E8','NeuN_3_A8','NeuN_3_C6','NeuN_3_C9','NeuN_3_D8','NeuN_3_E2','NeuN_3_E8'],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 0, 'MAX_HOM': 100}
#    },
#    'tree3':{
#        'samples':['Tcell_1_D6','Tcell_1_F5', 'Tcell_1_F6', 'Tcell_1_G7'],
#        'params': {'MIN_HET':2, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 0}
#    },
#    'tree4':{
#        'samples':['Monocyte_1_A3','Monocyte_1_A5', 'Monocyte_1_E3', 'Monocyte_1_H6'],
#        'params': {'MIN_HET':0, 'MIN_HOM': 0, 'MAX_HET': 0, 'MAX_HOM': 100}
#    },
#    'tree5':{
#        'samples':['MICROGLIA','Microglia_1_D6','Microglia_1_E8','Microglia_1_F3','Microglia_1_F6','Microglia_1_E5','Microglia_1_G2'],
#        'params': {'MIN_HET':0, 'MIN_HOM': 0, 'MAX_HET': 0, 'MAX_HOM': 100}
#
#   },
#   'all':{
#        'samples':['Cardio_1_D7','Cardio_1_F7','Cardio_1_F9','Cardio_1_G2','Cardio_1_G9','Cardio_1_H4','Cardio_2_C2','Cardio_2_G4','NEURONS','NeuN_2_D8','NeuN_2_E8','NeuN_3_A8','NeuN_3_C6','NeuN_3_C9','NeuN_3_D8','NeuN_3_E2','NeuN_3_E8','Tcell_1_D6','Tcell_1_F5', 'Tcell_1_F6', 'Tcell_1_G7','Monocyte_1_A3','Monocyte_1_A5', 'Monocyte_1_E3', 'Monocyte_1_H6','MICROGLIA','Microglia_1_D6','Microglia_1_E8','Microglia_1_F3','Microglia_1_F6','Microglia_1_E5','Microglia_1_G2'],
#       'params': { 'MIN_HET': 2, 'MIN_HOM': 1, 'MAX_HET': 100, 'MAX_HOM': 100}
#    }
#}



#trees = {   
#    'tree1':{
#        'samples':['Cardio_1_D7','Cardio_1_F7','Cardio_1_F9','Cardio_1_G2','Cardio_1_G9','Cardio_1_H4','Cardio_2_C2','Cardio_2_G4','NEURONS','NeuN_2_D8','NeuN_2_E8','NeuN_3_A8','NeuN_3_C6','NeuN_3_C9','NeuN_3_D8','NeuN_3_E2','NeuN_3_E8','Monocyte_1_A3','Monocyte_1_A5', 'Monocyte_1_E3', 'Monocyte_1_H6','MICROGLIA','Microglia_1_D6','Microglia_1_E8','Microglia_1_F3','Microglia_1_F6','Microglia_1_E5','Microglia_1_G2'],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 0, 'MAX_HOM': 100}
#    },
#    'tree2':{
#        'samples':['Tcell_1_D6','Tcell_1_F5', 'Tcell_1_F6', 'Tcell_1_G7'],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 100}
#    },
#    'all':{
#        'samples':['Cardio_1_D7','Cardio_1_F7','Cardio_1_F9','Cardio_1_G2','Cardio_1_G9','Cardio_1_H4','Cardio_2_C2','Cardio_2_G4','NEURONS','NeuN_2_D8','NeuN_2_E8','NeuN_3_A8','NeuN_3_C6','NeuN_3_C9','NeuN_3_D8','NeuN_3_E2','NeuN_3_E8','Monocyte_1_A3','Monocyte_1_A5', 'Monocyte_1_E3', 'Monocyte_1_H6','MICROGLIA','Microglia_1_D6','Microglia_1_E8','Microglia_1_F3','Microglia_1_F6','Microglia_1_E5','Microglia_1_G2','Tcell_1_D6','Tcell_1_F5', 'Tcell_1_F6', 'Tcell_1_G7'],
#        'params': { 'MIN_HET': 0, 'MIN_HOM': 0, 'MAX_HET': 100, 'MAX_HOM': 100}
#    }
#}







