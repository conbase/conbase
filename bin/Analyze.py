import json
from Stats import *
import pdb
from Params import analyze_params


    
def count_total_tuples(site_tuples):
    return sum([site_tuples['RR'], site_tuples['RA'], site_tuples['AR'], site_tuples['AA']])

def get_tuples_ratio(site_tuples, tuple_star):
    """ 
        get_tuples_ratio: calculates the ratio between tuples for a given tuple
        Args:
            site_tuples (dict): tuples for a Site object 
            tuple_star (str): tuple with the most votes (tp*) 
        Returns:
            het, homo_R, homo_A1 (str): ratios for het, homo_R, homo_A1 (may be None if count_total_tuples <= 0)
    """
    T = count_total_tuples(site_tuples)

    if T > 0:
        RR = float(site_tuples['RR'])/T
        RA = float(site_tuples['RA'])/T
        AR = float(site_tuples['AR'])/T
        AA = float(site_tuples['AA'])/T

        het1 = [RR, AA]
        het2 = [RA, AR]
        
        if tuple_star is None:
            het = ("HET-C1", 0, 0)
        else:
            if tuple_star == 'AA':
                het_ = (sum(het1), het1, "AA")
            elif tuple_star == 'AR':
                het_ = (sum(het2), het2, "AR")
                
            het = ("HET-C1", sum(het_[1]), ratio(het_[1][0], het_[1][1]), het_[2])

        homo_R = ("HOMO-C1", (RR+RA), ratio(RR,RA))
        homo_A1 = ("HOMO-A1", (AR+AA), ratio(AR,AA))
        site_tuples['stats'] = "H0:{homo_R}/{homo_R_ratio}, H:{het}/{het_ratio}, A1:{homo_A1}/{homo_A1_ratio}".format(homo_R= round(homo_R[1],2), homo_R_ratio=round(homo_R[2],2), het=round(het[1],2), het_ratio=round(het[2],2), homo_A1=round(homo_A1[1],2), homo_A1_ratio = round(homo_A1[2],2))
        
        site_tuples['het'] = "H0:{homo_R}/{homo_R_ratio}".format(homo_R= round(homo_R[1],2), homo_R_ratio=round(homo_R[2],2))
        site_tuples['homo_R'] ="H:{het}/{het_ratio}".format(het=round(het[1],2), het_ratio=round(het[2],2))
        site_tuples['homo_A1'] = "A1:{homo_A1}/{homo_A1_ratio}".format(homo_A1=round(homo_A1[1],2), homo_A1_ratio = round(homo_A1[2],2))
        return het, homo_R, homo_A1
    else:
        return None, None, None


def json_to_site(f):
    """ 
        json_to_site: parse json file to Site objects
        Args:
            f (file pointer): file pointer to json file  
        Returns:
            site (Site object): yields Site objects
    """
    for line in f:
        s = json.loads(line.strip('\n'))    
        samples = dict()
        sample_names = []
        for sample in s['samples'].values():
            snp_pos_list = list(sample['tuples'].keys())
            sample_names.append(sample['name'])
            site_tuples_dict = dict()
            for pos in snp_pos_list:
                site_tuples = dict()
                site_tuples['RR'] = sample['tuples'][pos]['RR']
                site_tuples['RA'] = sample['tuples'][pos]['RA']
                site_tuples['AR'] = sample['tuples'][pos]['AR']
                site_tuples['AA'] = sample['tuples'][pos]['AA']
                site_tuples['stats'] = sample['tuples'][pos]['stats']
                site_tuples_dict[int(pos)] = site_tuples

            samples[sample['name']] = Sample(sample['name'],sample['AD'], site_tuples_dict, sample['info'])
        site = Site(s['chrom'], s['pos'], s['ref'], s['alts'], s['kind'], samples, sample_names)
        site.bulk = s['bulk']
        yield site

def ratio(num1, num2):
    if (num1 + num2) > 0:
        return min(num1, num2)/(num1 + num2)
    else:
        return 0

def define_tuple_star(snp_pos, site):
    """ 
        define_tuple_star: defines tp* (see conbase paper)
        Args:
            snp_pos (int): position for snp
            site (Site object): Site object
        Returns:
            (str): tp* (may be None)
    """
    pairs = {'AR' : 0, 'AA' : 0}
    for sample in site.samples.values():
        if snp_pos in sample.tuples.keys() and count_total_tuples(sample.tuples[snp_pos]) >= analyze_params["dp_tuple_limit"]:
            site_tuples = sample.tuples[snp_pos]
            RR, RA, AR, AA = site_tuples['RR'], site_tuples['RA'], site_tuples['AR'], site_tuples['AA']
            
            min_, min_name = min(([AR, RA], 'AR'), ([RR, AA], 'AA'), key=lambda x: x[0][0]+x[0][1])
            max_, max_name = max(([AR, RA], 'AR'), ([RR, AA], 'AA'), key=lambda x: x[0][0]+x[0][1])
            
            if(ratio(sum(min_),sum(max_)) < analyze_params["tuple_group_ratio"] and ratio(max_[0], max_[1]) > analyze_params["win_internal_group_ratio"] ):
                #het gets to vote
                pairs[max_name] += 1
            elif(ratio(sum(min_),sum(max_)) < analyze_params["tuple_group_ratio"] and max((RR, 'RR'), (RA, 'RA'), (AR, 'AR'), (AA, 'AR'))[1] in {'AR', 'AA'}):
                #HET-C2 gets to vote
                pairs[max_name] += 1

    max_pair = sorted(list(pairs.items()), reverse=True, key=lambda x: x[1])
    if max_pair[0][1] >= analyze_params["sample_tuple_vote_limit"] and ratio(max_pair[0][1], max_pair[1][1]) <= (1-analyze_params["vote_tuple_ratio_limit"]):
        return max_pair[0][0]
    else:
        return None

def analyze(site):
    """ 
        analyze: run analyze (described in detail in the conbase paper)
        Args:
            site (Site object): Site object
    """
    for sample in site.samples.values():
        votes = {'HET-C1':0, 'HET-C2':0, 'HOMO-C1':0,'HOMO-C2':0, 'HOMO-A1':0, 'CONFLICT':0}
        nr_snp_allowed_voting = 0
        tuple_stars = dict()
        for snp_pos in sample.tuples.keys():    
            site_tuples = sample.tuples[snp_pos]
            tuple_star = define_tuple_star(snp_pos, site)
            tuple_stars[snp_pos] = tuple_star
            het, homo_R, homo_A1 = get_tuples_ratio(site_tuples, tuple_star)
            
            if count_total_tuples(site_tuples) >= analyze_params["dp_tuple_limit"]:
                tuple['voted'] = ''
                if tuple_star != None:
                    site.snp_tuple_star[snp_pos] = tuple_star              
                    nr_snp_allowed_voting += 1
                    ##CASE: HOMO-C2 or HET-C2
                    if het[2] < analyze_params["tuples_internal_ratio"] and homo_R[2] < analyze_params["tuples_internal_ratio"]:
                        if homo_R[1] >= analyze_params["tuples_ratio"] and het[1] < analyze_params["tuples_c2_external_error_ratio"]:
                            hr = max((tuple['RR'], 'RR'), (tuple['RA'], 'RA'))
                            if (hr[1] == 'RR' and tuple_star == 'AR') or (hr[1] =='RA' and tuple_star == 'AA'):
                                tuple['voted'] = 'HOMO-C2' #we know for certain this is not a mutation
                            else:
                                tuple['voted'] = 'UNKNOWN'
                        elif homo_R[1] < analyze_params["tuples_c2_external_error_ratio"] and het[1] >= analyze_params["tuples_ratio"]:
                            ha = max((tuple['AA'], 'AA'), (tuple['AR'], 'AR'))
                            if (ha[1] == 'AA' and tuple_star == 'AR') or (ha[1] =='AR' and tuple_star == 'AA'):
                                tuple['voted'] = 'CONFLICT'#this is contradictory
                            else:
                                tuple['voted'] = 'HET-C2'
                        else:
                            tuple['voted'] = 'UNKNOWN'

                    ##CASE: HET or HOMO-C1 or HOMO-A1 if it's only one true tuples
                    else:
                        tuples_list = [m for m in [het, homo_R, homo_A1] if m[1] >= analyze_params["tuples_ratio"] and m[2] >= analyze_params["tuples_internal_ratio"]]
                        if len(tuples_list) == 1:
                            if tuples_list[0][0] == 'HET-C1':
                                if tuples_list[0][3] == tuple_star:
                                    tuple['voted'] = 'HET-C1'   #mut-snp pair match!
                                else:                   
                                    tuple['voted'] = 'CONFLICT'              #match is not right
                            else:
                                tuple['voted'] = tuples_list[0][0]       #homo-C1 or homo-A1
                        elif len(tuples_list) > 1:
                            tuple['voted'] = 'CONFLICT'
                        else:
                            error_list = [m for m in [het, homo_R, homo_A1] if m[1] >= analyze_params["tuples_c2_external_error_ratio"] and m[2] >= analyze_params["tuples_internal_ratio"]]
                            if len(error_list) > 1:
                                tuple['voted'] = 'CONFLICT'
                            else:
                                tuple['voted'] = 'UNKNOWN'

                    if tuple['voted'] != '' and tuple['voted'] != 'UNKNOWN':
                        votes[tuple['voted']] += 1

        for k,v in list(votes.items()):
            if k == 'HET-C2' and votes['HOMO-A1'] == 0 and votes['HET-C1'] > 0:
                votes['HET-C1'] += v
                votes[k] = 0
            elif k == 'HET-C2' and votes['HET-C1'] == 0 and votes['HOMO-A1'] > 0 :
                votes['HOMO-A1'] += v
                votes[k] = 0
            elif k == 'HOMO-C2' and votes['HOMO-C1'] > 0:
                votes['HOMO-C1'] += v
                votes[k] = 0
        

        total_votes = sum(list(votes.values()))
        if total_votes >= 1:
            if votes['HOMO-C2'] != 0 and votes['HET-C2'] != 0:
                sample.info = 'CONFLICT'
            else:
                max_vote = sorted(votes.items(), reverse = True, key = lambda t: t[1])
                vote_limit = max_vote[0][1]/total_votes
                if vote_limit >= analyze_params["snp_total_vote"]:
                    if max_vote[0][0] == 'HOMO-C1' or max_vote[0][0] == 'HOMO-C2':
                        alts_dp = sum([v for k, v in sample.AD.items() if k != site.ref])

                        if sum(sample.AD.values()) > 0 and float(alts_dp)/sum(sample.AD.values()) <= analyze_params["homo_error_allowed"]:
                            if float(total_votes)/nr_snp_allowed_voting >= analyze_params["snp_vote_ratio"]:
                                sample.info = max_vote[0][0]
                            else:
                                sample.info = 'UNKNOWN'
                        else:
                            sample.info = 'CONFLICT'
                    else:
                        if float(total_votes)/nr_snp_allowed_voting >= analyze_params["snp_vote_ratio"]:
                            sample.info = max_vote[0][0]
                        else:
                            sample.info = 'UNKNOWN'
                else:
                    sample.info = 'CONFLICT'
        else:
            support_unmutated = 0
            for snp_pos in sample.tuples.keys():    
                site_tuples = sample.tuples[snp_pos]
                tuple_star = tuple_stars[snp_pos]
                if (tuple_star == 'AA' and site_tuples['RA'] >= analyze_params["c3_homo_limit"]) or (tuple_star == 'AR' and site_tuples['RR'] >= analyze_params["c3_homo_limit"]):
                    support_unmutated += 1
            if support_unmutated > 0:
                if sample.AD[site.alts['A1']] == 0:
                    sample.info = 'HOMO-C3'
                elif sample.AD[site.alts['A1']] >= analyze_params["c3_a1_limit"]:
                    sample.info = 'C3-CONFLICT'
                else:
                    sample.info = 'UNKNOWN'
            elif sample.AD[site.alts['A1']] >= analyze_params["c3_a1_limit"]:
                sample.info = 'HET-C3'
            else:
                if sum(sample.AD.values()) > 0:
                    sample.info = 'NOT-INFORMATIVE'
                else:
                    sample.info = 'ZERO-READS'

    
            