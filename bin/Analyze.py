import json
from Stats import *
import pdb
from Params import analyze_params


def json_to_site(f):
    for line in f:
        s = json.loads(line.strip('\n'))    
        samples = dict()
        sample_names = []
        for sample in s['samples'].values():
            snp_pos_list = list(sample['MSP'].keys())
            sample_names.append(sample['name'])
            MSP_dict = dict()
            for pos in snp_pos_list:
                msp_object = MSP()
                msp_object.RR = sample['MSP'][pos]['RR']
                msp_object.RA = sample['MSP'][pos]['RA']
                msp_object.AR = sample['MSP'][pos]['AR']
                msp_object.AA = sample['MSP'][pos]['AA']
                msp_object.stats = sample['MSP'][pos]['stats']
                MSP_dict[int(pos)] = msp_object

            samples[sample['name']] = Sample(sample['name'],sample['AD'], sample['GT'], None, MSP_dict, sample['info'])
        site = Site(s['CHROM'], s['POS'], s['REF'], s['ALTS'], s['TYPE'], samples, sample_names)
        site.BULK_INFO = s['BULK_INFO']
        yield site

def ratio(num1, num2):
    if (num1 + num2) > 0:
        return min(num1, num2)/(num1 + num2)
    else:
        return 0

def define_ms_pair(snp_pos, site):
    pairs = {'AR' : 0, 'AA' : 0}
    for sample in site.samples.values():
        if snp_pos in sample.MSP.keys() and sample.MSP[snp_pos].get_ms_total() >= analyze_params["dp_ms_limit"]:
            msp = sample.MSP[snp_pos]
            RR, RA, AR, AA = msp.RR, msp.RA, msp.AR, msp.AA
            
            min_, min_name = min(([AR, RA], 'AR'), ([RR, AA], 'AA'), key=lambda x: x[0][0]+x[0][1])
            max_, max_name = max(([AR, RA], 'AR'), ([RR, AA], 'AA'), key=lambda x: x[0][0]+x[0][1])
            
            if(ratio(sum(min_),sum(max_)) < analyze_params["ms_group_ratio"] and ratio(max_[0], max_[1]) > analyze_params["win_internal_group_ratio"] ):
                #het gets to vote
                pairs[max_name] += 1
            elif(ratio(sum(min_),sum(max_)) < analyze_params["ms_group_ratio"] and max((RR, 'RR'), (RA, 'RA'), (AR, 'AR'), (AA, 'AR'))[1] in {'AR', 'AA'}):
                #ado-a1 gets to vote
                pairs[max_name] += 1

    max_pair = sorted(list(pairs.items()), reverse=True, key=lambda x: x[1])
    if max_pair[0][1] >= analyze_params["sample_ms_vote_limit"] and ratio(max_pair[0][1], max_pair[1][1]) <= (1-analyze_params["vote_ms_ratio_limit"]):
        return max_pair[0][0]
    else:
        return None

def gt_ratio(site):
    for sample in site.samples.values():
        votes = {'HET':0, 'HOMO-R':0, 'HOMO-A1':0, 'ADO-R':0, 'ADO-A1':0, 'X':0}
        nr_snp_allowed_voting = 0
        ms_pairs = dict()
        for snp_pos in sample.MSP.keys():    
            msp = sample.MSP[snp_pos]
            ms_pair = define_ms_pair(snp_pos, site)
            ms_pairs[snp_pos] = ms_pair
            het, homo_R, homo_A1 = msp.get_msp_ratio(ms_pair)
            
            if msp.get_ms_total() >= analyze_params["dp_ms_limit"]:
                msp.voted = ''
                # if site.POS == 178756156 and sample.name == 'fibroblast_41' and snp_pos == 178756123:
                #     import pdb; pdb.set_trace()
                if ms_pair != None:
                    site.snp_ms_win[snp_pos] = ms_pair              
                    nr_snp_allowed_voting += 1
                    ##CASE: ADO-R or ADO-A1
                    if het[2] < analyze_params["msp_internal_ratio"] and homo_R[2] < analyze_params["msp_internal_ratio"]:
                        if homo_R[1] >= analyze_params["msp_ratio"] and het[1] < analyze_params["msp_c2_external_error_ratio"]:
                            hr = max((msp.RR, 'RR'), (msp.RA, 'RA'))
                            if (hr[1] == 'RR' and ms_pair == 'AR') or (hr[1] =='RA' and ms_pair == 'AA'):
                                msp.voted = 'ADO-R' #we know for certain this is not a mutation
                            else:
                                msp.voted = 'unknown'
                        elif homo_R[1] < analyze_params["msp_c2_external_error_ratio"] and het[1] >= analyze_params["msp_ratio"]:
                            ha = max((msp.AA, 'AA'), (msp.AR, 'AR'))
                            if (ha[1] == 'AA' and ms_pair == 'AR') or (ha[1] =='AR' and ms_pair == 'AA'):
                                msp.voted = 'X'#this is contradictory
                            else:
                                msp.voted = 'ADO-A1'
                        else:
                            msp.voted = 'unknown'

                    ##CASE: HET or HOMO-R or HOMO-A1 if it's only one true msp
                    else:
                        msp_list = [m for m in [het, homo_R, homo_A1] if m[1] >= analyze_params["msp_ratio"] and m[2] >= analyze_params["msp_internal_ratio"]]
                        if len(msp_list) == 1:
                            if msp_list[0][0] == 'HET':
                                if msp_list[0][3] == ms_pair:
                                    msp.voted = 'HET'   #mut-snp pair match!
                                else:                   
                                    msp.voted = 'X'              #match is not right
                            else:
                                msp.voted = msp_list[0][0]       #homo-R or homo-A1
                        elif len(msp_list) > 1:
                            msp.voted = 'X'
                        else:
                            error_list = [m for m in [het, homo_R, homo_A1] if m[1] >= analyze_params["msp_c2_external_error_ratio"] and m[2] >= analyze_params["msp_internal_ratio"]]
                            if len(error_list) > 1:
                                msp.voted = 'X'
                            else:
                                msp.voted = 'unknown'

                    if msp.voted != '' and msp.voted != 'unknown':
                        votes[msp.voted] += 1

        for k,v in list(votes.items()):
            if k == 'ADO-A1' and votes['HOMO-A1'] == 0 and votes['HET'] > 0:
                votes['HET'] += v
                votes[k] = 0
            elif k == 'ADO-A1' and votes['HET'] == 0 and votes['HOMO-A1'] > 0 :
                votes['HOMO-A1'] += v
                votes[k] = 0
            elif k == 'ADO-R' and votes['HOMO-R'] > 0:
                votes['HOMO-R'] += v
                votes[k] = 0
        

        total_votes = sum(list(votes.values()))
        if total_votes >= 1:
            if votes['ADO-R'] != 0 and votes['ADO-A1'] != 0:
                sample.info = 'X'
            else:
                max_vote = sorted(votes.items(), reverse = True, key = lambda t: t[1])
                vote_limit = max_vote[0][1]/total_votes
                if vote_limit >= analyze_params["snp_total_vote"]:
                    if max_vote[0][0] == 'HOMO-R' or max_vote[0][0] == 'ADO-R':
                        alts_dp = sum([v for k, v in sample.AD.items() if k != site.REF])

                        if sum(sample.AD.values()) > 0 and float(alts_dp)/sum(sample.AD.values()) <= analyze_params["homo_error_allowed"]:
                            if float(total_votes)/nr_snp_allowed_voting >= analyze_params["snp_vote_ratio"]:
                                sample.info = max_vote[0][0]
                            else:
                                sample.info = 'unknown'
                        else:
                            sample.info = 'X'
                    else:
                        if float(total_votes)/nr_snp_allowed_voting >= analyze_params["snp_vote_ratio"]:
                            sample.info = max_vote[0][0]
                        else:
                            sample.info = 'unknown'
                else:
                    #TODO Homo-R = ADO-R
                    sample.info = 'X'
        else:
            if sample.AD[site.ALTS['A1']] >= analyze_params["c2_a1_limit"]:
                sample.info = 'HET-C3'
            else:
                found = 0
                for snp_pos in sample.MSP.keys():    
                    msp = sample.MSP[snp_pos]
                    ms_pair = ms_pairs[snp_pos]
                    if (ms_pair == 'AA' and msp.RA >= analyze_params["c2_homo_limit"]) or (ms_pair == 'AR' and msp.RR >= analyze_params["c2_homo_limit"]):
                        found += 1 
                if found > 0:
                    if sample.AD[site.ALTS['A1']] == 0:
                        sample.info = 'HOMO-C3'
                    else:
                        sample.info = 'unknown'

