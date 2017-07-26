import sys, csv, pysam, os
from File_Output import *
import json
import multiprocessing as mp
from Params import stats_params

acceptable_bases = {'A','C','G','T'}
skipped_mate = 0

class JSONSerializable(object):
    def __repr__(self):
        return json.dumps(self.__dict__, default = lambda o: o.__dict__)

class Site(JSONSerializable):
    def __init__(self, CHROM, POS, REF, ALTS, TYPE, samples, sample_names):
        self.CHROM = CHROM
        self.POS = POS  #other index
        self.REF = REF
        self.ALTS = ALTS
        self.BULK_INFO = dict()
        self.TYPE = TYPE
        self.snp_ms_win = dict()
        self.samples = samples
        if not samples:
            self.init_samples(sample_names)

    def init_samples(self, sample_names):
        for sample_name in sample_names:
            self.samples[sample_name] = Sample(sample_name, {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}, dict(), None, dict(), 'None')

    def real_POS(self):
        return self.POS + 1

    def get_ALTS(self):
        alts = list(self.ALTS.keys())
        alts.sort()
        return ",".join([self.ALTS[alt] for alt in alts if self.ALTS[alt] != None])
    
    def get_bulk_info(self):
        if len(self.BULK_INFO) == 0:
            return ''
            
        bases = {'A', 'C', 'G', 'T'}
        info_str = "DP:{DP}, ".format(DP = self.BULK_INFO['SUM'])
        for b in bases:
            r = round(100*self.BULK_INFO[b],3)
            if r != 0:
                info_str += "{base}:{ratio}, ".format(base = b, ratio = r)
        return info_str

class Sample(JSONSerializable):
    def __init__(self, name, AD, GT, ST, MSP, info):
        self.name       = name
        self.AD         = AD
        self.GT         = GT
        self.ST         = ST
        self.snp_reads  = 0
        self.MSP        = MSP
        self.info       = info
        self.indels     = 0

    # TODO
    def get_AD(self, site):
        str_AD = ""
        alts = [k for k, v in site.ALTS.items() if v != None and self.AD[v] > 0]
        alts.sort()
        alts_AD = [k + ":" + str(self.AD[site.ALTS[k]]) for k in alts]

        if self.AD[site.REF] != 0 and len(alts) != 0:
            str_AD = 'R:{R}/{ALTS}'.format(R = self.AD[site.REF], ALTS = ",".join(alts_AD))
        elif self.AD[site.REF] == 0 and len(alts) != 0:
            str_AD = "./{ALTS}".format(ALTS = ",".join(alts_AD))
        elif self.AD[site.REF] != 0 and len(alts) == 0:
            str_AD = 'R:{R}/.'.format(R = self.AD[site.REF])
        elif self.AD[site.REF] == 0 and len(alts) == 0:
            str_AD = "./."
        # else:
        #     print("FEEEEEEEEEEELLLLLLLLLLLL")
        
        alts_value = [v for k, v in site.ALTS.items() if v != None and self.AD[v] > 0]
        other = sum([v for k, v in self.AD.items() if k != site.REF and k not in alts_value and v > 0])


        str_AD += ", O: {other}".format(other = other)
        return str_AD

    def get_GT(self):
        str_GT = ""
        for pos, msp in self.MSP.items():
            str_GT += "{pos}: {msp} | ".format(pos = int(pos)+1, msp = msp.get_stats())
        return str_GT

    def get_MSP(self):
        str_MSP = ""
        for pos, msp in self.MSP.items():
            str_MSP += "{pos}: RR:{RR}, RA:{RA}, AR:{AR}, AA:{AA} | ".format(pos = int(pos)+1, RR = str(msp.RR), RA = str(msp.RA), AR = str(msp.AR), AA = str(msp.AA))
        return str_MSP
    
class Read(object):
    def __init__(self, id, mate, sequence, tuples, start, end, base_quality, mapping_quality, has_snp):
        self.id = id
        self.mate = mate
        self.start = start
        self.end = end
        self.bases = self.init_bases(tuples, sequence)
        self.base_quality = self.init_base_quality(tuples, base_quality)
        self.mapping_quality = mapping_quality
        self.has_snp = has_snp

    def __str__(self):
        return "ID: {ID}, MATE: {MATE}, SEQ: {SEQ}, START: {START}, END: {END}".format(ID = self.id, MATE = self.mate.id, SEQ = self.bases, START = self.start, END = self.end)
    
    # D: AA<C>G -> [(0,100), (1,101), (None,102), (2,103)]
    # I: AA<C>G -> [(0,100), (1,101), (2,None), (3,102)]
    def init_bases(self,tuples,sequence):
        bases = dict()
        base_quality = dict()
        for i, p in tuples:
            if p != None and i != None:
                bases[p] = sequence[i]
            elif p != None and i == None: # D
                bases[p] = 'D'
            # elif p == None and i != None: #I
            #     bases[p] == 'I'

            # if self.start <= 76111775 and self.end >= 76111775 and p==None:
        return bases 

    def init_base_quality(self, tuples, base_quality_list):
        base_quality = dict()
        for i, p in tuples:
            if p != None and i != None:
                base_quality[p] = base_quality_list[i] 
            elif p != None and i == None:
                base_quality[p] = 0      
        return base_quality 

class MSP(JSONSerializable):
    def __init__(self):
        self.RR = 0
        self.RA = 0
        self.AR = 0
        self.AA = 0
        self.stats = ''
        self.voted = ''
        self.het = None
        self.homo_R = None
        self.homo_A1 = None
        
    def get_ms_total(self):
        return sum([self.RR,self.RA,self.AR,self.AA])

    def get_msp_ratio(self, ms_pair):
        T = self.get_ms_total()

        if T > 0:
            RR = float(self.RR)/T
            RA = float(self.RA)/T
            AR = float(self.AR)/T
            AA = float(self.AA)/T

            het1 = [RR, AA]
            het2 = [RA, AR]
            
            if ms_pair is None:
                het = ("HET", 0, 0)
            else:
                if ms_pair == 'AA':
                    het_ = (sum(het1), het1, "AA")
                elif ms_pair == 'AR':
                    het_ = (sum(het2), het2, "AR")
                    
                het = ("HET", sum(het_[1]), ratio(het_[1][0], het_[1][1]), het_[2])

            homo_R = ("HOMO-R", (RR+RA), ratio(RR,RA))
            homo_A1 = ("HOMO-A1", (AR+AA), ratio(AR,AA))
            self.stats = "H0:{homo_R}/{homo_R_ratio}, H:{het}/{het_ratio}, A1:{homo_A1}/{homo_A1_ratio}".format(homo_R= round(homo_R[1],2), homo_R_ratio=round(homo_R[2],2), het=round(het[1],2), het_ratio=round(het[2],2), homo_A1=round(homo_A1[1],2), homo_A1_ratio = round(homo_A1[2],2))
            
            self.het = "H0:{homo_R}/{homo_R_ratio}".format(homo_R= round(homo_R[1],2), homo_R_ratio=round(homo_R[2],2))
            self.homo_R ="H:{het}/{het_ratio}".format(het=round(het[1],2), het_ratio=round(het[2],2))
            self.homo_A1 = "A1:{homo_A1}/{homo_A1_ratio}".format(homo_A1=round(homo_A1[1],2), homo_A1_ratio = round(homo_A1[2],2))
            return het, homo_R, homo_A1
        else:
            return None, None, None
            
    def get_stats(self):
        return self.stats

def get_reads(snp, bams):
    sample_reads = dict()
    for bam_name, bam_file in bams:
        reads = list()
        mates = dict()
        left = snp.POS-stats_params["fragment_length"]
        right = snp.POS+stats_params["fragment_length"]

        if left < 0:
            left = 0
#        print(bam_name)

        for read in bam_file.fetch(snp.CHROM,  left, right):
            if read.mapping_quality >= stats_params["mapping_quality"] and read.is_paired and read.is_proper_pair:
                r = Read(read.query_name, None, read.query_sequence, read.get_aligned_pairs(),
                    read.reference_start, read.reference_end-1, read.query_qualities, read.mapping_quality, False)
                if snp.POS in r.bases.keys():
                    r.has_snp = True
                try:
                    m = mates[r.id]
                    if m.mate != None:
                        global skipped_mate
                        skipped_mate += 1
                        continue
                    r.mate = m
                    m.mate = r
                    if r.has_snp or r.mate.has_snp:
                        snp.samples[bam_name].snp_reads += 1
                        r.has_snp = True
                        r.mate.has_snp = True
                except KeyError:
                    mates[r.id] = r

                reads.append(r)
        sample_reads[bam_name] = reads
    return sample_reads

def get_references(chrom, start, end, ref_file):
    try:
        sequence = ref_file.fetch(chrom, start, end+1)
    except KeyError:
        sequence = ref_file.fetch("chr"+chrom, start, end+1)
    aligned_pairs = dict()
    for i in range(0, end-start+1):
        pos = i + start
        base = sequence[i]
        if base not in acceptable_bases:
            base = None
        aligned_pairs[pos] = base
    return aligned_pairs

def init_site(snp, sample_names, reference, pos):
    ref = reference[pos]
    if pos == snp.POS:
        site = snp
    else:
        site = Site(snp.CHROM, pos, ref, None, '', dict(), sample_names)
    return site

# snp_limit returns min and max position for all reads belonging to a snp, otherwise none
def snp_limits(snp, reads):
    start = list()
    end = list()
    found_snp = False
    for sample_name, sample_reads in reads.items():
        if snp.samples[sample_name].snp_reads > stats_params["snp_read_limit"]:
            found_snp = True
            start_sample = list()
            end_sample = list()
            for read in sample_reads:
                if read.has_snp:
                    start_sample.append(read.start)
                    end_sample.append(read.end)
            start.append(min(start_sample))
            end.append(max(end_sample))
    if found_snp:
        return min(start), max(end)
    else:
        return None, None

def allele_counter(reads, site, pos):
    for sample_name, sample_reads in reads.items():
        for read in sample_reads:
            if pos in read.bases.keys():
                base = read.bases[pos].upper()
                if base in acceptable_bases and read.bases[pos] != None and read.base_quality[pos] > stats_params["base_quality"]:
                    site.samples[sample_name].AD[base] += 1
                elif base in {'D', 'I'}:
                     site.samples[sample_name].indels += 1

def is_indel(site):
    tot_indel_ratio = 0.0
    for sample in site.samples.values():
        tot_indel_ratio += float(sample.indels)/sum(sample.AD.values()) if sum(sample.AD.values()) > 0 else 0
    return (tot_indel_ratio/len(site.samples)) > stats_params["indel_ratio"]

def ratio(num1, num2):
    if (num1 + num2) > 0:
        return min(num1, num2)/(num1 + num2)
    else:
        return 0.0

def define_altenative(site):
    if site.ALTS == None:
        allele_vote = {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}
        for sample in site.samples.values():
            lst = []
            for b in sample.AD.keys(): 
                if b != site.REF:
                    if ratio(sample.AD[b], sample.AD[site.REF]) >= stats_params["alt_ratio_limit"] and \
                     (sample.AD[site.REF]+sample.AD[b]) >= stats_params["dp_limit"]:
                        lst.append((b, sample.AD[b]))
                    elif ratio(sample.AD[b], sample.AD[site.REF]) <= stats_params["alt_ratio_limit"] and sample.AD[b] >= stats_params["dp_limit"]:
                        lst.append((b, sample.AD[b]))

            if lst != []:
                lst.sort(key = lambda t: t[1], reverse=True)
                #R:10/A1:5, A2:5 <-- will not be allowed to vote
                if len(lst) == 1 or (lst[0][1] > lst[1][1]):
                    allele_vote[lst[0][0]] += 1

        max_vote_allele = sorted([b for b in allele_vote.items() if b[0] != site.REF], reverse = True, key = lambda t: t[1])
        sum_votes = sum([b[1] for b in max_vote_allele])
        vote_ratio = 0
        if sum_votes > 0:
            vote_ratio = max_vote_allele[0][1]/sum_votes

        site.ALTS = {'A1': None, 'A2': None, 'A3':None}

        if vote_ratio >= stats_params["vote_ratio_limit"] and len(max_vote_allele) >= stats_params["sample_vote_limit"]:
            i = 1
            for b in max_vote_allele:
                if b[1] != 0:
                    alt_str = "A" + str(i)
                    site.ALTS.update({alt_str:b[0]})
                    i += 1
        else:
            site.TYPE = 'UNDEF'
        alts = [alt for alt in site.ALTS.values() if alt != None]
        if len(alts) == 0:
            site.TYPE = 'HOMO'

def mut_snp(snp, sites, reads):
    for sample_name, sample_reads in reads.items():
        for read in sample_reads:
            if read.has_snp:
                snp_read = read
                if snp.POS not in snp_read.bases.keys():
                    snp_read = read.mate
                    if snp.POS not in snp_read.bases.keys():
                        print(snp_read.id)
                        print("Base not in snp mate pos", snp.POS)

                for pos in read.bases.keys():
                    if pos in sites.keys() and sites[pos].TYPE != 'SNP' and read.base_quality[pos] >= stats_params["base_quality"]:
                        count_MS(sites[pos], snp, read.bases[pos], snp_read.bases[snp.POS], sample_name)


def count_MS(site, snp, site_base, snp_base, sample_name):
    if snp.POS not in site.samples[sample_name].MSP.keys():
        site.samples[sample_name].MSP[snp.POS] = MSP()

    if site_base == site.REF and snp_base == snp.REF:
            site.samples[sample_name].MSP[snp.POS].RR += 1
    elif site_base == site.REF and snp_base == snp.ALTS['A1']:
            site.samples[sample_name].MSP[snp.POS].RA += 1
    elif site_base == site.ALTS['A1'] and snp_base == snp.REF:
            site.samples[sample_name].MSP[snp.POS].AR += 1
    elif site_base == site.ALTS['A1'] and snp_base == snp.ALTS['A1']:
            site.samples[sample_name].MSP[snp.POS].AA += 1
    ## todo: count A2, A3...


def bulk_stats(site, bulk_bam):
    sample_reads = dict()

    bases = {"A":0, "C":0, "G":0, "T":0}
    for read in bulk_bam.fetch(site.CHROM,  site.POS, site.POS+1):
        if read.mapping_quality >= stats_params["mapping_quality"]:
            r = Read(read.query_name, None, read.query_sequence, read.get_aligned_pairs(),
                read.reference_start, read.reference_end-1, read.query_qualities, read.mapping_quality, False)
            if site.POS in r.base_quality.keys() and r.base_quality[site.POS] >= stats_params["base_quality"]:
                base = r.bases[site.POS]
                if base in bases.keys():
                    bases[base] += 1
    
    
    T = sum(list(bases.values()))
    if T > 0:
        site.BULK_INFO['A'] = bases['A']
        site.BULK_INFO['C'] = bases['C']
        site.BULK_INFO['G'] = bases['G']
        site.BULK_INFO['T'] = bases['T']
        site.BULK_INFO['SUM'] = T           
        
        bulk_bases = 0
        if site.TYPE != 'SNP' and float(site.BULK_INFO[site.REF])/T < stats_params["bulk_ref_limit"]:
            site.TYPE = 'E'
        # for v in site.BULK_INFO.values():
        #     if v > 0:
        #         bulk_bases += 1
        # if bulk_bases > 1 and site.TYPE != 'SNP':
        #     site.TYPE = 'E'
        # elif bulk_bases == 1 and site.TYPE != 'SNP':
        #     max_base = max(list(site.BULK_INFO.items()), key=lambda t: t[1])[0]
        #     if max_base != site.REF:
        #         site.TYPE = 'E'


def stats_to_json(i, snps_chunk_path, bams_path, sample_names, reference_path, output_name, queue):
    bams, bam_bulk = get_bams(bams_path)
    reference_genome_file = pysam.Fastafile(reference_path)

    sites = dict()
    old_end = 0
    json_path = './.conbase/' + output_name + '_chunk_' + str(i) + '.json'
    json_file = Site2JSON(json_path)

    snps_reader = csv.DictReader(open(snps_chunk_path, 'rU'), delimiter='\t')
    for row in snps_reader:
        snp = Site(row['CHROM'], int(row['POS']) - 1, (row['REF'].strip()), {"A1":row['ALT'].strip()}, 'SNP', dict(), sample_names)
        
        reads = get_reads(snp, bams)
        new_start, new_end = snp_limits(snp, reads)
        
        if (new_start == None and new_end == None):            
            for s in sites.values():
                if s.TYPE == '':
                    bulk_stats(s, bam_bulk)
                json_file.write(s)
            sites = dict()
            old_end = 0
            continue
        
        if new_start - old_end > stats_params["fragment_length"]:
            for s in sites.values():
                if s.TYPE == '':
                    bulk_stats(s, bam_bulk)
                json_file.write(s)
            sites = dict()
        
        reference = get_references(snp.CHROM, new_start, new_end, reference_genome_file)
        for pos in range(new_start, new_end+1):
            if pos not in sites.keys() and reference[pos] != None:
                site = init_site(snp, sample_names, reference, pos) 
                allele_counter(reads, site, pos) 
                if not is_indel(site):
                    sites[pos] = site
                    define_altenative(sites[pos])
        mut_snp(snp, sites, reads)
        old_end = new_end
        
        queue.put(1)
        
    for s in sites.values():
        if s.TYPE == '':
            bulk_stats(s, bam_bulk)
        json_file.write(s)
    json_file.close()
        
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def get_bams(bam_paths):
    bam_reader = csv.DictReader(open(bam_paths, 'rU'), delimiter='\t')
    bams = []
    bam_bulk = None
    for row in bam_reader:
        # print(row['NAME'])
        if row['NAME'] == 'BULK':
            bam_bulk = pysam.AlignmentFile(row['PATH'], 'rb')
        else:
            bams.append((row['NAME'], pysam.AlignmentFile(row['PATH'], 'rb')))  
    return bams, bam_bulk

def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

def new_chunk_file(chunk_number, output_name):
    chunk_path = './.conbase/' + output_name + '_snp_chunk_' + str(chunk_number) + '.tsv'
    chunk_file = open(chunk_path, 'w')
    chunk_file.write('CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\n')
    return chunk_file, chunk_path


def snps_to_chunks(snps_path, nodes, output_name):
    print('Loading SNPS ...')
    if nodes == 1:
        return [snps_path]

    f = open(snps_path, 'r')
    snp_count = sum([ bl.count("\n") for bl in blocks(f)]) - 1
    chunk_size = int(snp_count/nodes)

    current_chunk_number = 0
    current_chunk_file, current_chunk_path = new_chunk_file(current_chunk_number, output_name)
    chunks_path = [current_chunk_path]
    
    i = 0
    prev_row_pos = None
    prev_row_chrom = None

    snp_reader = csv.DictReader(open(snps_path, 'rU'), delimiter='\t')
    for row in snp_reader:
        current_row_pos = int(row['POS'])
        current_row_chrom = int(row['CHROM'])

        if i > chunk_size:
            if (prev_row_pos != None and prev_row_chrom != None) and abs(current_row_pos - prev_row_pos) >= stats_params["fragment_length"]*2 or current_row_chrom != prev_row_chrom:
                i = 0
                current_chunk_file.close()

                current_chunk_number += 1
                current_chunk_file, current_chunk_path = new_chunk_file(current_chunk_number, output_name)
                chunks_path.append(current_chunk_path)

        current_chunk_file.write(row['CHROM'] + '\t' + row['POS'] + '\t' + row['REF'] + '\t' + row['ALT'] + '\n')
        prev_row_pos = current_row_pos
        prev_row_chrom = current_row_chrom
        i += 1

    current_chunk_file.close()
    return chunks_path, snp_count

def get_sample_names(bam_paths):
    print('Loading BAMS ...')
    bam_reader = csv.DictReader(open(bam_paths, 'rU'), delimiter='\t')
    sample_names = []
    for row in bam_reader:
        print(row['NAME'])
        if row['NAME'] != 'BULK':
            sample_names.append(row['NAME'])
    return sample_names  

def progress_bar(nr_snps, queue, bar_width=100):
    sys.stdout.write("[{}]".format(" " * bar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (bar_width+1))
    counter = 0
    prev_progress = 0
    while True :
        if queue.get() == 'Done':
            sys.stdout.write("\n")
            break
        else:
            counter += 1
            current_progress = int(counter/nr_snps*bar_width)
            if current_progress != prev_progress:
                sys.stdout.write("#" * (current_progress - prev_progress))
                sys.stdout.flush()
                prev_progress = current_progress

def stats(snps_path, bam_paths, reference_path, nodes, output_name): 
    if not os.path.exists("./.conbase"):
        os.makedirs("./.conbase")

    if not os.path.exists('../results/'):
        os.makedirs('../results/')

    # os.system("rm ./.conbase/chunk_*")	
    # os.system("rm ./.conbase/snp_chunk_*")
    sample_names = get_sample_names(bam_paths)
    snps_chunks_path, nr_snps = snps_to_chunks(snps_path, nodes, output_name)
    nr_chunks = len(snps_chunks_path)

    jobs = []
    queue = mp.Queue()
    for i, snps_chunk_path in enumerate(snps_chunks_path):
        p = mp.Process(target=stats_to_json, args=(i, snps_chunk_path, bam_paths, sample_names, reference_path, output_name, queue))
        jobs.append(p)
        p.start()
    
    p = mp.Process(target=progress_bar, args=(nr_snps, queue))
    jobs.append(p)
    p.start()

    for i, job in enumerate(jobs):
        if i == len(jobs) - 1:
            queue.put('Done')
        job.join()
        
    print('All done')

    f = open( '../results/' + output_name + '.json', 'w')
    f.write('{' + '"samples":' + json.dumps(sample_names) + '}\n')
    f.write('{' + '"stats_params":' + json.dumps(stats_params) + '}\n')
    f.close()

    for i in range(nr_chunks):
        f = './.conbase/' + output_name + '_chunk_' + str(i) + '.json'
        os.system('cat '+f+' >> ../results/' + output_name + '.json')
    os.system("rm ./.conbase/" + output_name + "_chunk_*")	
    os.system("rm ./.conbase/" + output_name + "_snp_chunk_*")