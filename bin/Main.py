import Stats
import Analyze
import File_Output
import sys
import csv
import json

if __name__ == '__main__':
    tag = sys.argv[1]

    if tag == 'stats':
        snps_path = sys.argv[2]
        bam_paths = sys.argv[3]
        reference_path = sys.argv[4]
        nodes = int(sys.argv[5])
        output_file_name = sys.argv[6]
        Stats.stats(snps_path, bam_paths, reference_path, nodes, output_file_name)
        print("Done!", + '../results/' + output_file_name + ".json" )

    elif tag == 'analyze':
        json_path = sys.argv[2]
        output_name = '../results/' + sys.argv[3]
        f = open(json_path)
        samples_names = json.loads(f.readline().strip('\n'))['samples']

        html = File_Output.HTML(output_name + ".html", samples_names)
        phylip = File_Output.Phylip_format(output_name + ".txt", samples_names)
        sites = Analyze.json_to_site(f)
        my_sites = []
        for site in sites:
            Analyze.gt_ratio(site)
            for sample in site.samples.values():
                if site.TYPE == '' and (sample.info == 'HET' or sample.info == 'ADO-A1' or sample.info == 'HOMO-A1'):
                    print(site.CHROM + ':' + str(site.real_POS()))
                    my_sites.append(site)
                    break
        
        my_sites.sort(key= lambda o: (int(o.CHROM), int(o.POS)))
        for site in my_sites:
            html.write_site(site)
        phylip.write_sites(my_sites)
        html.close()

