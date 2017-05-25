import csv
import Stats
import Params as params
import pysam
import multiprocessing as mp
import argparse

acceptable_bases = {'A','C','G','T'}

def chrom_alt_sites(chrom, bam_path, reference_path):
    bam_file = pysam.AlignmentFile(bam_path, 'rb')
    reference_genome_file = pysam.Fastafile(reference_path)
    chrom_sites = dict()
    for read in bam_file.fetch(str(chrom)):
        if read.mapping_quality >= params.mapping_quality and read.is_paired and read.is_proper_pair:
            r = Stats.Read(read.query_name, None, read.query_sequence, read.get_aligned_pairs(),
                read.reference_start, read.reference_end-1, read.query_qualities, read.mapping_quality, False)
            
            for pos in r.bases.keys():
                if r.base_quality[pos] > params.base_quality:
                    if pos not in chrom_sites.keys():
                        chrom_sites[pos] = {'A':0, 'C':0, 'G':0, 'T':0}
                    chrom_sites[pos][r.bases[pos].upper()] += 1
    
    reference = Stats.get_references(str(chrom), min(chrom_sites.keys()), max(chrom_sites.keys()), reference_genome_file)
    pos_list = list(chrom_sites.keys())
    i = 0
    for pos in pos_list:
        ref = reference[pos]
        T = sum(chrom_sites[pos].values())
        if ref not in acceptable_bases or float(chrom_sites[pos][ref])/T >= params.bulk_ref_limit:
            chrom_sites.pop(pos)
            i += 1
    print("ref", str(i), alt, len(chrom_sites))
    return chrom_sites

def chrom_duplicate_region(chrom, bulk_path, reference_path, queue, alt_nr_limit=10, alt_dist_limit=10000):
    chrom_alt_path = './.conbase/duplicate_region_' + chrom + '.txt'
    site_writer = open(chrom_alt_path, 'w')

    chrom_sites = chrom_alt_sites(chrom, bulk_path, reference_path)
    pos_list = sorted(list(chrom_sites.keys()))
    left_pos = None
    right_pos = None
    alt_list = []
    for i, pos in enumerate(pos_list):
        if right_pos == None or pos == right_pos:
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
    queue.put(chrom_alt_path)


def duplicate_regions(bulk_path, reference_path, output_name="duplicate_regions.txt"):
    jobs = []
    i = 0
    queue = mp.Queue()
    for chrom in range(1,23):
        p = mp.Process(target=stats_to_json, args=(str(chrom), bulk_path, reference_path, queue))
        i += 1
        jobs.append(p)
        p.start()

    for job in jobs:
        job.join()

    chrom_duplicate_regions_list = []
    while not queue.empty():
        chrom_duplicate_regions_list.append(queue.get())
    print('all done')

    for chrom in range(1,23):
        f = './.conbase/duplicate_region_' + chrom + '.txt'
        os.system('cat '+f+' >> ../results/' + output_name + '.json')
    os.system("rm ./.conbase/duplicate_region_*")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run conbase')
    parser.add_argument('--duplicate_regions', nargs=3, metavar=("<bam path>", "<reference path>", "<output name>"))
    parser.add_argument('--analyze', nargs=2, metavar=("<json path>", "<output name>"))
    parser.add_argument('--params', nargs=1, metavar=("<params path>"))
    args = parser.parse_args()

    if args.duplicate_regions is not None:
        duplicate_regions(args.duplicate_regions[0], args.duplicate_regions[1], args.duplicate_regions[2])


# bulk_path = "/media/box2/Experiments/Joanna/Snake_analys/j_frisen_1602/Fibs/Tree2/FibBulk/FibBulk.reAligned.bwa.bam"
# reference_path = "/media/box2/reference_assemblies/bundle/2.8/b37/from_pall/human_g1k_v37.fasta"
# chrom = '1'
# chrom_sites = chrom_alt_sites(chrom, bulk_path, reference_path)
# chrom_duplicate_region(chrom,chrom_sites.keys())