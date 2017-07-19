############### Stats #################
# Read
fragment_length = 600
mapping_quality = 20
base_quality = 20
# ALT
dp_limit = 15
alt_ratio_limit = 0.3
sample_vote_limit = 2   
vote_ratio_limit = 0.7
snp_read_limit = 1
# Indel
indel_ratio = 0.05
# Bulk
bulk_ref_limit = 0.91 #this requires the remaining percent to be only a1 
############### Analyze #################
dp_ms_limit = 10
snp_total_vote = 0.9
snp_vote_ratio = 0.9

msp_ratio = 0.9
msp_internal_ratio = 0.05
msp_c2_external_error_ratio = 1-msp_ratio

c2_a1_limit = 3
c2_homo_limit = 3

homo_error_allowed = 0

ms_group_ratio = 0.01
win_internal_group_ratio = 0.1

sample_ms_vote_limit = 2
vote_ms_ratio_limit = 0.9

conflicting_upper_limit = 20
a1_lower_limit = 2
bulk_dp_interval = (15,45) # Must have the following format (min, max)
################ Misc ##################
snp_nr_limit = 10
snp_dist_limit = 1000


