import Stats, Analyze, File_Output
from Params import analyze_params, misc_params, stats_params
import argparse, json
import Misc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run conbase')
    parser.add_argument('--stats', nargs=5, metavar=("<snp path>", "<bam path>", "<reference path>", "<number of nodes>", "<output name>"))
    parser.add_argument('--analyze', nargs=2, metavar=("<json path>", "<output name>"))
    parser.add_argument('--params_stats', nargs=1, metavar=("<params dict>"))
    parser.add_argument('--params_analyze', nargs=1, metavar=("<params dict>"))

    args = parser.parse_args()

    if args.params_analyze is not None:
        analyze_params.update(eval(args.params_analyze[0]))
        
    if args.params_stats is not None:
        stats_params.update(eval(args.params_stats[0]))

    if args.stats is not None:
        snps_path = args.stats[0]

        bam_paths = args.stats[1]
        reference_path = args.stats[2]
        nodes = int(args.stats[3])
        output_file_name = args.stats[4]
        Stats.stats(snps_path, bam_paths, reference_path, nodes, output_file_name)
        print("Done!" + '../results/' + output_file_name + ".json" )

    if args.analyze is not None:
        print('run...')
        json_path = args.analyze[0]
        # main_directory = os.path.dirname(os.path.realpath(__file__)) 
        output_name =  '../results/' + args.analyze[1]
        f = open(json_path)
        samples_names = json.loads(f.readline().strip('\n'))['samples']
        stats_params_in_json = json.loads(f.readline().strip('\n'))['stats_params']

        html = File_Output.HTML(output_name + ".html", samples_names, stats_params_in_json, analyze_params, misc_params)
        tsv = File_Output.TSV(output_name + ".tsv", samples_names)

        sites = Analyze.json_to_site(f)
        my_sites = []
        for site in sites:
            Analyze.gt_ratio(site)
            if site.kind == '' and len(site.bulk) > 0 and site.bulk['SUM'] >= analyze_params["bulk_dp_interval"][0] and site.bulk['SUM'] <= analyze_params["bulk_dp_interval"][1]:
                bulk_a1_ratio = float(site.bulk[site.ALTS['A1']])/site.bulk['SUM']
                if bulk_a1_ratio <= (1 - stats_params_in_json["bulk_ref_limit"]):
                    nr_conflicting = 0
                    nr_c3_conflicting = 0
                    nr_a1 = 0
                    for sample in site.samples.values():
                        if sample.info == 'CONFLICT':
                            nr_conflicting += 1
                        elif sample.info == 'C3-CONFLICT':
                            nr_c3_conflicting += 1
                        elif (sample.info == 'HET-C1' or sample.info == 'HET-C2' or sample.info == 'HOMO-A1'):
                            nr_a1 += 1

                    if nr_a1 >= analyze_params["a1_lower_limit"] and nr_conflicting <= analyze_params["conflicting_upper_limit"] and nr_c3_conflicting <= analyze_params["c3_conflicting_upper_limit"]:
                        print(site.chrom + ':' + str(site.true_pos))
                        my_sites.append(site)
        
        my_sites.sort(key= lambda o: (int(o.chrom), int(o.pos)))
        
        my_sites = Misc.check_duplicate_region(my_sites)
        my_sites = Misc.trees_stats(my_sites, output_name)

        for site in my_sites:
            html.write_site(site)
            tsv.write_site(site)
        html.close()
        tsv.close()

        table_plot = File_Output.TABLE_PLOT(output_name + ".pdf", samples_names, my_sites, a1_param=True, border_color='none') # border_color='none'
        
        
        print('Done! ' + output_name + '.html')

