import os
import multiprocessing as mp
import argparse



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Experiment')
    parser.add_argument('--run', nargs=0, metavar=(""))
    parser.add_argument('--plot', nargs=1, metavar=("<txt path>"))
    args = parser.parse_args()

    if args.run is not None:
        analyze_params.update(eval(args.params_analyze[0]))

        tcell_path = '/media/box2/Experiments/Marez/conbase/results/female_final.json'
        fibs_path = '/media/box2/Experiments/Marez/conbase/results/fibs/fibs_29juli_wo4243.json'
        test_path = '/Users/ezeddin/Projects/work/conbase_ki/results/test_marez_params.json'

        def get_key(dp, internal_ratio, external_ratio, dataset):
            return dataset + str(dp) + str(internal_ratio) + str(external_ratio)

        # bash_code = 'source activate snowflake; python Main.py --analyze {link} {output_name} --params_analyze {param_dict}'.format(link=test_path, output_name='fibs_29juli_wo4243_dp_5', param_dict='\'"dp_ms_limit":{dp}, "msp_ratio":{ratio}\''.format{dp=5,ratio=0.8})
        def run_bash_code(json_path, dp, internal_ratio, external_ratio, dataset, queue):
            param_dict = '"dp_ms_limit":{dp}, "msp_internal_ratio":{internal_ratio}, "msp_ratio":{external_ratio}'.format(dp=dp,internal_ratio=internal_ratio, external_ratio=external_ratio)
            param_dict = '\'{' + param_dict + '}\''
            bash_code = 'python Main.py --analyze {path} {output_name} --params_analyze {param_dict}'.format(path=json_path, output_name=get_key(dp, internal_ratio, external_ratio, dataset), param_dict=param_dict)
            os.system(bash_code)
            queue.put(bash_code)

        dp_range = [5, 10, 15, 20, 30, 40, 70]
        internal_ratio_range = [0.1, 0.2, 0.3, 0.4]
        external_ratio_range = [0.3, 0.5, 0.7, 0.8, 0.9]
        dataset_name = 'fibs'

        jobs = []
        queue = mp.Queue()
        for dp in dp_range:
            for internal_ratio in internal_ratio_range:
                for external_ratio in external_ratio_range:
                    p = mp.Process(target=run_bash_code, args=(fibs_path, dp, internal_ratio, external_ratio, dataset_name, queue))
                    jobs.append(p)
                    p.start()

        for job in jobs:
            job.join()

        while not queue.empty():
            queue.get()
        print('all done')

        output_path = '../results/' + dataset_name + '_experiment.txt'
        f = open(output_path, 'w')
        f.close()

        for dp in dp_range:
            for internal_ratio in internal_ratio_range:
                for external_ratio in external_ratio_range:
                    file_path =  '../results/' + get_key(dp, internal_ratio,external_ratio, dataset_name) + '_trees_stats.txt'
                    os.system('cat ' + file_path + ' >> ' + output_path)
                    os.system('rm ' + file_path)
    
    if args.plot is not None:
        pass