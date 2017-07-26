import csv
import Stats
from Params import misc_params, stats_params
import pysam
import multiprocessing as mp
import argparse
import os

acceptable_bases = {'A','C','G','T'}

def check_duplicate_region(sites):
    
    filtered_sites = list()
    if len(sites) > misc_params["mut_nr_limit"]:
        sites_per_chrome = dict()
        for s in sites:
            if s.CHROM not in sites_per_chrome.keys():
                sites_per_chrome[s.CHROM] = list()
            sites_per_chrome[s.CHROM].append(s)

        for chrom_sites in sites_per_chrome.values():
            if len(chrom_sites) > misc_params["mut_nr_limit"]:
                for index, site in enumerate(chrom_sites):
                    
                    left = index-misc_params["mut_nr_limit"] if index-misc_params["mut_nr_limit"] >= 0 else 0
                    right = index+misc_params["mut_nr_limit"]+1 if index+misc_params["mut_nr_limit"]+1 <= len(chrom_sites) else len(chrom_sites)
                    
                    window_of_sites = chrom_sites[left:right]
                    in_duplicate_region = False

                    for i in range(len(window_of_sites)-misc_params["mut_nr_limit"]):
                        interval = window_of_sites[i:i+misc_params["mut_nr_limit"] + 1]
                        if max(interval, key=lambda s: s.POS).POS - min(interval, key=lambda s: s.POS).POS <= misc_params["mut_dist_limit"]:
                            in_duplicate_region = True
                    if not in_duplicate_region:
                        filtered_sites.append(site)
                    else:
                        print('removed site! position:', site.CHROM, site.real_POS())
            else:
                filtered_sites += chrom_sites
    else:
        filtered_sites = sites

    return filtered_sites

def snp_in_duplicate_region(snp, bam_file, reference_genome_file):
    sites = dict()

    left = snp['POS']-misc_params["snp_dist_limit"] if snp['POS']-misc_params["snp_dist_limit"] > 0 else 0
    right = snp['POS']+misc_params["snp_dist_limit"]

    for read in bam_file.fetch(snp['CHROM'], left, right):
        if read.mapping_quality >= stats_params["mapping_quality"] and read.is_paired and read.is_proper_pair:
            r = Stats.Read(read.query_name, None, read.query_sequence, read.get_aligned_pairs(),
                read.reference_start, read.reference_end-1, read.query_qualities, read.mapping_quality, False)

            for pos in r.bases.keys():
                if pos >= left and pos <= right and r.base_quality[pos] > stats_params["base_quality"]:
                    if pos not in sites.keys():
                        sites[pos] = {'A':0, 'C':0, 'G':0, 'T':0}
                    sites[pos][r.bases[pos].upper()] += 1

    reference = Stats.get_references(snp['CHROM'], left, right, reference_genome_file)
    pos_list = list(sites.keys())
    for pos in pos_list:
        ref = reference[pos]
        T = sum(sites[pos].values())
        if ref not in acceptable_bases or float(sites[pos][ref])/T >= stats_params["bulk_ref_limit"]:
            sites.pop(pos)
    
    pos_list = sorted(list(sites.keys()))
    in_duplicate_region = False
    if len(pos_list) > misc_params["snp_nr_limit"]:
        for i in range(len(pos_list)-misc_params["snp_nr_limit"] + 1):
            interval = pos_list[i:i+misc_params["snp_nr_limit"]]
            if max(interval) - min(interval) <= misc_params["snp_dist_limit"]:
                in_duplicate_region = True
                break
    return in_duplicate_region

def SNP_duplicate_region(snp_path, bam_path, reference_path, queue):
    SNP_reader = open(snp_path, 'r')
    SNP_writer = open(snp_path[:-4] + '_not_duplicate_region_.tsv', 'w')

    snps = []
    SNP_reader.readline()
    for line in SNP_reader:
        CHROM, POS, REF, ALT = line.rstrip('\n').strip().split('\t')
        snps.append({'CHROM':CHROM, 'POS':int(POS), 'REF':REF, 'ALT':ALT})

    bam_file = pysam.AlignmentFile(bam_path, 'rb')
    reference_genome_file = pysam.Fastafile(reference_path)
    for snp in snps:
        if not snp_in_duplicate_region(snp, bam_file, reference_genome_file):
            SNP_writer.write(snp['CHROM'] + '\t' + str(snp['POS']) + '\t' +  snp['REF'] + '\t' + snp['ALT'] + '\n')
    queue.put(snp_path)

def duplicate_regions(snps_path, bam_path, reference_path, nodes=1, output_name="duplicate_regions"):

    if not os.path.exists("./.conbase"):
        os.makedirs("./.conbase")
    if not os.path.exists("../results"):
        os.makedirs("../results")

    # os.system("rm ./.conbase/duplicate_region_*")
    # os.system("rm ./.conbase/" + output_name + "_snp_chunk_*")

    snps_chunks_path, _ = Stats.snps_to_chunks(snps_path, int(nodes), output_name)

    jobs = []
    queue = mp.Queue()
    for snps_chunk_path in snps_chunks_path:
        p = mp.Process(target=SNP_duplicate_region, args=(snps_chunk_path, bam_path, reference_path, queue))
        jobs.append(p)
        p.start()

    for job in jobs:
        job.join()

    while not queue.empty():
        queue.get()
    print('all done')

    f = open( '../results/' + output_name + '.tsv', 'w')
    f.write('CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\n')
    f.close()

    for snps_chunk_path in snps_chunks_path:
        f = snps_chunk_path[:-4] + '_not_duplicate_region_.tsv'
        os.system('cat '+f+' >> ../results/' + output_name + '.tsv')
    os.system("rm ./.conbase/duplicate_region_*")
    os.system("rm ./.conbase/" + output_name + "_snp_chunk_*")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Conbase preprocessing-tool for removing duplicate regions')
    parser.add_argument('--duplicate_regions', nargs=5, metavar=("<snp path>", "<bam path>", "<reference path>", "<number of nodes>", "<output name>"))
    args = parser.parse_args()
    if args.duplicate_regions is not None:
        duplicate_regions(*args.duplicate_regions)


# bulk_path = "/media/box2/Experiments/Joanna/Snake_analys/j_frisen_1602/Fibs/Tree2/FibBulk/FibBulk.reAligned.bwa.bam"
# reference_path = "/media/box2/reference_assemblies/bundle/2.8/b37/from_pall/human_g1k_v37.fasta"
# chrom = '1'
# chrom_sites = chrom_alt_sites(chrom, bulk_path, reference_path)
# chrom_duplicate_region(chrom,chrom_sites.keys())
