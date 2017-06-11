import json
from Stats import *
import pdb
import Params as params

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
        if snp_pos in sample.MSP.keys() and sample.MSP[snp_pos].get_ms_total() >= params.dp_ms_limit:
            msp = sample.MSP[snp_pos]
            RR, RA, AR, AA = msp.RR, msp.RA, msp.AR, msp.AA
            
            min_, min_name = min(([AR, RA], 'AR'), ([RR, AA], 'AA'), key=lambda x: x[0][0]+x[0][1])
            max_, max_name = max(([AR, RA], 'AR'), ([RR, AA], 'AA'), key=lambda x: x[0][0]+x[0][1])
            
            if(ratio(sum(min_),sum(max_)) < params.ms_group_ratio and ratio(max_[0], max_[1]) > params.win_internal_group_ratio ):
                #het gets to vote
                pairs[max_name] += 1
            elif(ratio(sum(min_),sum(max_)) < params.ms_group_ratio and max((RR, 'RR'), (RA, 'RA'), (AR, 'AR'), (AA, 'AR'))[1] in {'AR', 'AA'}):
                #ado-a1 gets to vote
                pairs[max_name] += 1

    max_pair = sorted(list(pairs.items()), reverse=True, key=lambda x: x[1])
    if max_pair[0][1] > params.sample_ms_vote_limit and ratio(max_pair[0][1], max_pair[1][1]) < params.vote_ms_ratio_limit:
        return max_pair[0][0]
    else:
        return None

def gt_ratio(site):
    for sample in site.samples.values():
        vote = {'HET':0, 'HOMO-R':0, 'HOMO-A1':0, 'ADO-R':0, 'ADO-A1':0, 'X':0}
        for snp_pos in sample.MSP.keys():    
            msp = sample.MSP[snp_pos]
            
            ms_pair = define_ms_pair(snp_pos, site)
            het, homo_R, homo_A1 = msp.get_msp_ratio(ms_pair)
            
            if msp.get_ms_total() >= params.dp_ms_limit:
                msp.voted = True #TODO
                site.snp_ms_win[snp_pos] = ms_pair                    
                if ms_pair != None:
                    ##CASE: ADO-R or ADO-A1        
                    if homo_R[2] <= 0.01 and het[2] <= 0.01 and homo_A1[2] <= 0.01:
                        if homo_R[1] >= 0.99 and homo_A1[1] <= 0.01:
                            hr = max((msp.RR, 'RR'), (msp.RA, 'RA'))
                            if (hr[1] == 'RR' and ms_pair == 'AR') or (hr[1] =='RA' and ms_pair == 'AA'):
                                vote['ADO-R'] += 1 #we know for certain this is not a mutation
        
                        elif homo_A1[1] >= 0.99 and homo_R[1] <= 0.01:
                            ha = max((msp.AA, 'AA'), (msp.AR, 'AR'))
                            if (ha[1] == 'AA' and ms_pair == 'AR') or (ha[1] =='AR' and ms_pair == 'AA'):
                                vote['X'] += 1 #this is contradictory
                            else:
                                vote['ADO-A1'] += 1
                        else:
                            msp.voted = False
                    ##CASE: HET or HOMO-R or HOMO-A1 if it's only one true msp
                    else:   
                        msp_list = [m for m in [het, homo_R, homo_A1] if m[1] >= params.msp_ratio and m[2] >= params.msp_internal_ratio]
                        if len(msp_list) == 1:
                            if msp_list[0][0] == "HET":
                                if msp_list[0][3] == ms_pair:
                                    vote[msp_list[0][0]] += 1   #mut-snp pair match!
                                else:                   
                                    vote['X'] += 1              #match is not right
                            else:
                                vote[msp_list[0][0]] += 1       #homo-R or homo-A1
                        elif len(msp_list) > 1:
                            vote['X'] += 1
                        else:
                            msp.voted = False
                    #if msp.voted:
                    #    print('new snp, sample:',sample.name , vote)
                else:
                    msp_list = [m for m in [het, homo_R, homo_A1] if m[1] >= params.msp_ratio and m[2] >= params.msp_internal_ratio]
                    if len(msp_list) == 1:
                        if msp_list[0][0] != "HET" and msp_list[0][0] != "HOMO-A1":
                            vote[msp_list[0][0]] += 1       #homo-R or homo-A1
                    elif len(msp_list) > 1:
                        vote['X'] += 1
                    else:
                        msp.voted = False                 
        
        total_votes = sum(list(vote.values()))
        if total_votes >= 1:
            if vote['ADO-R'] != 0 and vote['ADO-A1'] != 0:
                sample.info = 'X'

            max_vote = sorted(vote.items(), reverse = True, key = lambda t: t[1])
            max_list = [max_vote[0][0], max_vote[1][0]]
            vote_limit = max_vote[0][1]/total_votes
            if vote_limit > params.snp_total_vote:
                sample.info = max_vote[0][0]  
            else:
                sample.info = 'X'
