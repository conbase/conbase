import sys, csv, pysam, os
from File_Output import SiteToJSON
import json
import multiprocessing as mp
from Params import stats_params


class JSONSerializable(object):
    def __repr__(self):
        return json.dumps(self.__dict__, default = lambda o: o.__dict__)

class Site(JSONSerializable):

    """
        Class Site is a representation of a genomic site and may hold information relating to that site,
        such as the reference base and alternative bases (from BAM file), type of site (SNP, undetermined, etc)
        and the Sample objects associated with the specific site. 
    """
    def __init__(self, chrom, pos, ref, alts, kind, samples, sample_names):

        """
            Args:
                chrom (str): chromosome
                pos (int): genomic position (zero-indexed)
                ref (str): reference base (e.g. A, C, G, T)
                alts (dict): dict of alternative bases (A1, A2, A3) where the read depth of A1 >= A2 >= A3 (eg. {'A1' : 'C', 'A2' : 'G', 'A3' : 'T'} if ref is 'A')
                kind (str): type of site (e.g. SNP, homozygous, heterozygous, undefined)
                samples (dict): mapping of a sample name to a Sample object
                sample_names (list): list of sample names in the same order as the BAM file
            
            Attributes:
                true_pos (int): one-indexed genomic position
                bulk (dict): dictionary of read depth for bases in bulk (e.g. {'A' : 30, 'C' : 40, 'G' : 0, 'T' : 0, 'SUM' : 70})
                snp_ms_win (dict): ??
        """
        self.chrom = chrom
        self.pos = pos 
        self.ref = ref
        self.alts = alts
        self.kind = kind
        self.samples = samples
        if not samples:
            self.init_samples(sample_names)

        self.true_pos = pos + 1  
        self.bulk = dict()
        self.snp_ms_win = dict()

    def init_samples(self, sample_names):
        for sample_name in sample_names:
            self.samples[sample_name] = Sample(sample_name, {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}, dict(), 'None')

class Sample(JSONSerializable):
    def __init__(self, name, AD, tuples, info):
        self.name       = name
        self.AD         = AD
        self.snp_reads  = 0
        self.tuples     = tuples
        self.info       = info
        self.indels     = 0

    def get_AD(self, site):
        str_AD = ""
        alts = [k for k, v in site.alts.items() if v != None and self.AD[v] > 0]
        alts.sort()
        alts_AD = [k + ":" + str(self.AD[site.alts[k]]) for k in alts]

        if self.AD[site.ref] != 0 and len(alts) != 0:
            str_AD = 'R:{R}/{alts}'.format(R = self.AD[site.ref], alts = ",".join(alts_AD))
        elif self.AD[site.ref] == 0 and len(alts) != 0:
            str_AD = "./{alts}".format(alts = ",".join(alts_AD))
        elif self.AD[site.ref] != 0 and len(alts) == 0:
            str_AD = 'R:{R}/.'.format(R = self.AD[site.ref])
        elif self.AD[site.ref] == 0 and len(alts) == 0:
            str_AD = "./."

        alts_value = [v for k, v in site.alts.items() if v != None and self.AD[v] > 0]
        other = sum([v for k, v in self.AD.items() if k != site.ref and k not in alts_value and v > 0])


        str_AD += ", O: {other}".format(other = other)
        return str_AD

class Read(object):
    def __init__(self, id, mate, sequence, ind_pos, start, end, base_quality, mapping_quality, has_snp):
        self.id = id
        self.mate = mate
        self.start = start
        self.end = end
        self.bases = self.init_bases(ind_pos, sequence)
        self.base_quality = self.init_base_quality(ind_pos, base_quality)
        self.mapping_quality = mapping_quality
        self.has_snp = has_snp

    def __str__(self):
        return "ID: {ID}, MATE: {MATE}, SEQ: {SEQ}, START: {START}, END: {END}".format(ID = self.id, MATE = self.mate.id, SEQ = self.bases, START = self.start, END = self.end)
    
    # D: AA<C>G -> [(0,100), (1,101), (None,102), (2,103)]
    # I: AA<C>G -> [(0,100), (1,101), (2,None), (3,102)]
    def init_bases(self,ind_pos,sequence):
        bases = dict()
        base_quality = dict()
        for i, p in ind_pos:
            if p != None and i != None:
                bases[p] = sequence[i]
            elif p != None and i == None: # D
                bases[p] = 'D'
            # elif p == None and i != None: #I
            #     bases[p] == 'I'

            # if self.start <= 76111775 and self.end >= 76111775 and p==None:
        return bases 

    def init_base_quality(self, ind_pos, base_quality_list):
        base_quality = dict()
        for i, p in ind_pos:
            if p != None and i != None:
                base_quality[p] = base_quality_list[i] 
            elif p != None and i == None:
                base_quality[p] = 0      
        return base_quality 

def get_reads(snp, bams):
    sample_reads = dict()
    for bam_name, bam_file in bams:
        reads = list()
        mates = dict()
        left = max(0, snp.pos-stats_params["fragment_length"])
        right = snp.pos+stats_params["fragment_length"]
        # TODO try except (get_ref..)
        for read in bam_file.fetch(snp.chrom, left, right):
            if read.mapping_quality >= stats_params["mapping_quality"] and read.is_paired and read.is_proper_pair:
                r = Read(read.query_name, None, read.query_sequence, read.get_aligned_pairs(),
                    read.reference_start, read.reference_end-1, read.query_qualities, read.mapping_quality, False)
                if snp.pos in r.bases.keys():
                    r.has_snp = True
                try:
                    m = mates[r.id]
                    if m.mate != None:
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
        if base not in stats_params["acceptable_bases"]:
            base = None
        aligned_pairs[pos] = base
    return aligned_pairs

def init_site(snp, sample_names, reference, pos):
    ref = reference[pos]
    if pos == snp.pos:
        site = snp
    else:
        site = Site(snp.chrom, pos, ref, None, '', dict(), sample_names)
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
                if base in stats_params["acceptable_bases"] and read.bases[pos] != None and read.base_quality[pos] > stats_params["base_quality"]:
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
    if site.alts is None:
        allele_vote = {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}
        for sample in site.samples.values():
            lst = []
            for b in sample.AD.keys(): 
                if b != site.ref:
                    if ratio(sample.AD[b], sample.AD[site.ref]) >= stats_params["alt_ratio_limit"] and \
                     (sample.AD[site.ref]+sample.AD[b]) >= stats_params["dp_limit"]:
                        lst.append((b, sample.AD[b]))
                    elif ratio(sample.AD[b], sample.AD[site.ref]) <= stats_params["alt_ratio_limit"] and sample.AD[b] >= stats_params["dp_limit"]:
                        lst.append((b, sample.AD[b]))

            if lst != []:
                lst.sort(key = lambda t: t[1], reverse=True)
                #R:10/A1:5, A2:5 <-- will not be allowed to vote
                if len(lst) == 1 or (lst[0][1] > lst[1][1]):
                    allele_vote[lst[0][0]] += 1

        max_vote_allele = sorted([b for b in allele_vote.items() if b[0] != site.ref], reverse = True, key = lambda t: t[1])
        sum_votes = sum([b[1] for b in max_vote_allele])
        vote_ratio = 0
        if sum_votes > 0:
            vote_ratio = max_vote_allele[0][1]/sum_votes

        site.alts = {'A1' : None, 'A2': None, 'A3' : None}
        if vote_ratio >= stats_params["vote_ratio_limit"] and sum_votes >= stats_params["sample_vote_limit"]:
            i = 1
            for b in max_vote_allele:
                if b[1] != 0:
                    alt_str = "A" + str(i)
                    site.alts.update({alt_str:b[0]})
                    i += 1
        else:
            site.kind = 'UNDEF'
        alts = [alt for alt in site.alts.values() if alt != None]
        if len(alts) == 0:
            site.kind = 'HOMO'

def count_tuple(site, snp, site_base, snp_base, sample_name):
    if snp.pos not in site.samples[sample_name].tuples.keys():
        site.samples[sample_name].tuples[snp.pos] = dict()

    if site_base == site.ref and snp_base == snp.ref:
            site.samples[sample_name].tuples[snp.pos]['RR'] += 1
    elif site_base == site.ref and snp_base == snp.alts['A1']:
            site.samples[sample_name].tuples[snp.pos]['RA'] += 1
    elif site_base == site.alts['A1'] and snp_base == snp.ref:
            site.samples[sample_name].tuples[snp.pos]['AR'] += 1
    elif site_base == site.alts['A1'] and snp_base == snp.alts['A1']:
            site.samples[sample_name].tuples[snp.pos]['AA'] += 1

def tuple_counter(snp, sites, reads):
    for sample_name, sample_reads in reads.items():
        for read in sample_reads:
            if read.has_snp:
                snp_read = read
                if snp.pos not in snp_read.bases.keys():
                    snp_read = read.mate
                    if snp.pos not in snp_read.bases.keys():
                        print(snp_read.id)
                        print("Base not in snp mate pos", snp.pos)

                for pos in read.bases.keys():
                    if pos in sites.keys() and sites[pos].kind != 'SNP' and read.base_quality[pos] >= stats_params["base_quality"]:
                        count_tuple(sites[pos], snp, read.bases[pos], snp_read.bases[snp.pos], sample_name)

def bulk_stats(site, bulk_bam):
    sample_reads = dict()

    bases = {"A" : 0, "C" : 0, "G" : 0, "T" : 0}
    for read in bulk_bam.fetch(site.chrom,  site.pos, site.pos+1):
        if read.mapping_quality >= stats_params["mapping_quality"]:
            r = Read(read.query_name, None, read.query_sequence, read.get_aligned_pairs(),
                read.reference_start, read.reference_end-1, read.query_qualities, read.mapping_quality, False)
            if site.pos in r.base_quality.keys() and r.base_quality[site.pos] >= stats_params["base_quality"]:
                base = r.bases[site.pos]
                if base in bases.keys():
                    bases[base] += 1

    T = sum(list(bases.values()))
    if T > 0:
        site.bulk['A'] = bases['A']
        site.bulk['C'] = bases['C']
        site.bulk['G'] = bases['G']
        site.bulk['T'] = bases['T']
        site.bulk['SUM'] = T           
        
        bulk_bases = 0
        if site.kind != 'SNP' and float(site.bulk[site.ref])/T < stats_params["bulk_ref_limit"]:
            site.kind = 'ERROR'

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

def stats_to_json(i, snps_chunk_path, bams_path, sample_names, reference_path, output_name, queue):
    bams, bam_bulk = get_bams(bams_path)
    reference_genome_file = pysam.Fastafile(reference_path)

    sites = dict()
    old_end = 0
    json_path = './.conbase/' + output_name + '_chunk_' + str(i) + '.json'
    json_file = SiteToJSON(json_path)

    snps_reader = csv.DictReader(open(snps_chunk_path, 'rU'), delimiter='\t')
    for row in snps_reader:
        snp = Site(row['CHROM'], int(row['POS']) - 1, (row['REF'].strip()), {"A1":row['ALT'].strip()}, 'SNP', dict(), sample_names)
        
        reads = get_reads(snp, bams)
        new_start, new_end = snp_limits(snp, reads)
        
        if (new_start == None and new_end == None):            
            for s in sites.values():
                if s.kind == '':
                    bulk_stats(s, bam_bulk)
                json_file.write(s)
            sites = dict()
            old_end = 0
            continue
        
        if new_start - old_end > stats_params["fragment_length"]:
            for s in sites.values():
                if s.kind == '':
                    bulk_stats(s, bam_bulk)
                json_file.write(s)
            sites = dict()
        
        reference = get_references(snp.chrom, new_start, new_end, reference_genome_file)
        for pos in range(new_start, new_end+1):
            if pos not in sites.keys() and reference[pos] != None:
                site = init_site(snp, sample_names, reference, pos) 
                allele_counter(reads, site, pos) 
                if not is_indel(site):
                    sites[pos] = site
                    define_altenative(sites[pos])
        tuple_counter(snp, sites, reads)
        old_end = new_end
        
        queue.put(1)
        
    for s in sites.values():
        if s.kind == '':
            bulk_stats(s, bam_bulk)
        json_file.write(s)
    json_file.close()
        
# Yield successive n-sized chunks from l.
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

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

    # TODO: print that remove files
    os.system("rm ./.conbase/" + output_name + "_chunk_*")	
    os.system("rm ./.conbase/" + output_name + "_snp_chunk_*")

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
