import Stats, Analyze, File_Output
import Params as params
import argparse, json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run conbase')
    parser.add_argument('--stats', nargs=5, metavar=("<snp path>", "<bam path>", "<reference path>", "<number of nodes>", "<output name>"))
    parser.add_argument('--analyze', nargs=2, metavar=("<json path>", "<output name>"))
    parser.add_argument('--params', nargs=1, metavar=("<params path>"))
    args = parser.parse_args()

    if args.params is not None:
        f = open(args.params[0], 'r')
        lines = f.read().splitlines()
        for line in lines:
            attribute, _, value_str =  line.split()
            try:
                value = int(value_str)
            except ValueError:
                value = float(value_str)
            setattr(params, attribute, value)

    if args.stats is not None:
        snps_path = args.stats[0]
        bam_paths = args.stats[1]
        reference_path = args.stats[2]
        nodes = int(args.stats[3])
        output_file_name = args.stats[4]
        Stats.stats(snps_path, bam_paths, reference_path, nodes, output_file_name)
        print("Done!", + '../results/' + output_file_name + ".json" )

    if args.analyze is not None:
        json_path = args.analyze[0]
        # main_directory = os.path.dirname(os.path.realpath(__file__)) 
        output_name =  '../results/' + args.analyze[1]
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
