import csv

def remove_snp_duplicate_region(snps_path, snp_nr_limit=10, snp_dist_limit=10000):
    snp_reader = csv.DictReader(open(snps_path, 'rU'), delimiter='\t')
    snps_path_new = snps_path + '.new'
    snp_writer = open(snps_path_new, 'w')
    snp_writer.write('CHROM\tPOS\n')

    first_row = snp_reader.__next__()
    left_pos = int(first_row['POS'])
    prev_chrom = first_row['CHROM']

    snp_list = [(prev_chrom, left_pos)]

    def empty_snp_list(chrom, pos, snp_list):
        if len(snp_list) < snp_nr_limit:
            for s_chrom, s_pos in snp_list:
                snp_writer.write(s_chrom + '\t' + str(s_pos) + '\n')
        left_pos = pos
        snp_list = [(chrom, pos)]
        return left_pos, snp_list

    for row in snp_reader:
        pos = int(row['POS'])
        chrom = row['CHROM']

        if prev_chrom == chrom:
            diff = pos - left_pos
            if diff < snp_dist_limit:
                snp_list.append((chrom, pos))
            else:
                left_pos, snp_list = empty_snp_list(chrom, pos, snp_list)
        else:
            left_pos, snp_list = empty_snp_list(chrom, pos, snp_list)
   
        prev_chrom = chrom
    empty_snp_list(prev_chrom, left_pos, snp_list)
    snp_writer.close()
    return snps_path_new

