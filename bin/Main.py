import Stats, Analyze, File_Output
from Params import analyze_params, misc_params, stats_params
import argparse, json, sys
import Misc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run conbase')
    parser.add_argument('--stats', nargs=5, metavar=("<snp path>", "<bam path>", "<reference path>", "<number of nodes>", "<outprefix>"))
    parser.add_argument('--analyze', nargs=2, metavar=("<json path>", "<outprefix>"))
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
        outprefix = args.stats[4]
        Stats.stats(snps_path, bam_paths, reference_path, nodes, outprefix)
        print("Done! results in " + outprefix + ".json" )

    if args.analyze is not None:
        print('run...')
        json_path = args.analyze[0]
        outprefix =  args.analyze[1]
        f = open(json_path)
        samples_names = json.loads(f.readline().strip('\n'))['samples']
        stats_params_in_json = json.loads(f.readline().strip('\n'))['stats_params']

        html = File_Output.HTML(outprefix + ".html", samples_names, stats_params_in_json, analyze_params, misc_params)
        tsv = File_Output.TSV(outprefix + ".tsv", samples_names)

        sites = Analyze.json_to_site(f)
        my_sites = []
        for site in sites:
            Analyze.analyze(site)
            if site.kind == '' and len(site.bulk) > 0 and site.bulk['SUM'] >= analyze_params["bulk_dp_interval"][0] and site.bulk['SUM'] <= analyze_params["bulk_dp_interval"][1]:
                bulk_a1_ratio = float(site.bulk[site.alts['A1']])/site.bulk['SUM']
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
                        #print(site.chrom + ':' + str(site.true_pos))
                        my_sites.append(site)
        
        my_sites.sort(key= lambda o: (int(o.chrom), int(o.pos)))
        # filtering gives error if no sites are found, create empty files instead.
        if len(my_sites) == 0:
            print("OBS! No sites found! Creating empty output files:")
            print(outprefix)
            open(outprefix + ".pdf", 'w').close()
            open(outprefix + ".html", 'w').close()                                            
            sys.exit(0)
            
        print("Read sites: %d"%len(my_sites))
        my_sites = Misc.check_duplicate_region(my_sites)
        print("After removed duplicates: %d"%len(my_sites))
        my_sites = Misc.filter_by_trees(my_sites)
        print("After filter by trees: %d"%len(my_sites))

	# does not produce pdf if no sites found, need to create empty files.
        if len(my_sites) == 0:
            print("OBS! No sites found after filtering! Creating empty output files:")
            print(outprefix)
            open(outprefix + ".pdf", 'w').close()
            open(outprefix + ".html", 'w').close()            
            sys.exit(0)
        
        for site in my_sites:
            html.write_site(site)
            tsv.write_site(site)
        html.close()
        tsv.close()

        table_plot = File_Output.TABLE_PLOT(outprefix + ".pdf", samples_names, my_sites, a1_param=True, border_color='none') # border_color='none'
        
        print('Done! results in:')
        print( outprefix + '.html')
        print( outprefix + '.pdf')
        print( outprefix + '.tsv')                

