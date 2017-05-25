import csv
import Stats
import Params as params
import pysam
acceptable_bases = {'A','C','G','T'}

def chrom_alt_sites(chrom, bam_path, reference_path):
    bam_file = pysam.AlignmentFile(bam_path, 'rb')
    reference_genome_file = pysam.Fastafile(reference_path)
    chrom_sites = dict()
    for read in bam_file.fetch(str(chrom),  all):
        if read.mapping_quality >= params.mapping_quality and read.is_paired and read.is_proper_pair:
            r = Stats.Read(read.query_name, None, read.query_sequence, read.get_aligned_pairs(),
                read.reference_start, read.reference_end-1, read.query_qualities, read.mapping_quality, False)
            
            for pos in read.bases.keys():
                if read.base_quality[pos] > params.base_quality:
                    if pos not in chrom_sites.keys():
                        chrom_sites[pos] = {'A':0, 'C':0, 'G':0, 'T':0}
                    chrom_sites[pos][read.bases[pos].upper()] += 1
                    chrom_sites[pos]['REF'] = reference[pos]
    
    reference = Stats.get_references(str(chrom), min(chrom_sites.keys()), max(chrom_sites.keys()), reference_genome_file)
    for pos in chrom_sites.keys():
        ref = reference[pos]
        T = sum(chrom_sites[pos].values())
        if ref not in acceptable_bases or float(chrom_sites[pos][ref])/T >= params.bulk_ref_limit:
            chrom_sites.pop(pos)
    return chrom_sites

def chrom_duplicate_region(chrom, pos_list, alt_nr_limit=10, alt_dist_limit=10000):
    chrom_alt_path = chrom + '.txt'
    site_writer = open(chrom_alt_path, 'w')

    pos_list = sorted(list(pos_list))
    left_pos = None
    right_pos = None
    alt_list = []
    for i, pos in enumerate(pos_list):
        for near_pos in pos_list[i:]:
            diff = near_pos - pos
            if diff <= alt_dist_limit:
                alt_list.append(near_pos)
            else:
                break
        if len(alt_list) >= alt_nr_limit:
            if left_pos == None:
                left_pos = alt_list[0]
            right_pos = alt_list[-1]
        else:
            if right_pos != None:
                site_writer.write(chrom + ':' + str(left_pos) + '-' + str(right_pos) + '\n')
            left_pos = None
            right_pos = None
        alt_list = []
                
    
    if len(alt_list) >= alt_nr_limit:
        if left_pos == None:
            left_pos = alt_list[0]
        right_pos = alt_list[-1]
    if right_pos != None:
        site_writer.write(chrom + ':' + str(left_pos) + '-' + str(right_pos) + '\n')

    site_writer.close()
    return chrom_alt_path

bulk_path = "/media/box2/Experiments/Joanna/Snake_analys/j_frisen_1602/Fibs/Tree2/FibBulk/FibBulk.reAligned.bwa.bam"
reference_path = "/media/box2/reference_assemblies/bundle/2.8/b37/from_pall/human_g1k_v37.fasta"
chrom = '1'
chrom_sites = chrom_alt_sites(chrom, bulk_path, reference_path)
chrom_duplicate_region(chrom,chrom_sites[chrom].keys())

# chrom_duplicate_region("1",[100,200,300,400, 900,1000, 1800,1900,2000,2001, 2299, 2500, 2502, 3000], alt_nr_limit=4, alt_dist_limit=500)