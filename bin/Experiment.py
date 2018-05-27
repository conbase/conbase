import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Experiment')
    parser.add_argument('-r', '--run', action='store_true')
    parser.add_argument('-p', '--plot', nargs=1, metavar=("<txt path>"))
    args = parser.parse_args()

    if args.run:
        import os
        import multiprocessing as mp

        tcell_path = '/media/box2/Experiments/Marez/conbase/results/female_final.json'
        fibs_path = '/media/box2/Experiments/Marez/conbase/results/fibs/fibs_29juli_wo4243.json'
        test_path = '/Users/ezeddin/Projects/work/conbase_ki/results/test_marez_params.json'

        def get_key(dp, internal_ratio, external_ratio, dataset):
            return dataset + str(dp) + str(internal_ratio) + str(external_ratio)

        # bash_code = 'source activate snowflake; python Main.py --analyze {link} {output_name} --params_analyze {param_dict}'.format(link=test_path, output_name='fibs_29juli_wo4243_dp_5', param_dict='\'"dp_tuple_limit":{dp}, "tuples_ratio":{ratio}\''.format{dp=5,ratio=0.8})
        def run_bash_code(json_path, dp, internal_ratio, external_ratio, dataset, queue):
            param_dict = '"dp_tuple_limit":{dp}, "tuples_internal_ratio":{internal_ratio}, "tuples_ratio":{external_ratio}'.format(dp=dp,internal_ratio=internal_ratio, external_ratio=external_ratio)
            param_dict = '\'{' + param_dict + '}\''
            bash_code = 'python Main.py --analyze {path} {output_name} --params_analyze {param_dict}'.format(path=json_path, output_name=get_key(dp, internal_ratio, external_ratio, dataset), param_dict=param_dict)
            os.system(bash_code)
            queue.put(bash_code)

        dp_range = [5, 10, 15, 20, 30, 40, 70]
        internal_ratio_range = [0.1, 0.2, 0.3, 0.4]
        external_ratio_range = [0.3, 0.5, 0.7, 0.8, 0.9]
        # dp_range = [5, 10, 15]
        # internal_ratio_range = [0.1, 0.2]
        # external_ratio_range = [0.3, 0.5]
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
        import numpy as np
        from numpy import genfromtxt
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt

        def plot(data, filter, title):
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            line = ax.scatter(data[:,0], data[:,1], data[:,2], c=filter)
            cb = plt.colorbar(line, label=title)
            ax.set_xlabel('dp limit')
            ax.set_ylabel('ms internal')
            ax.set_zlabel('ms external')

            def forceUpdate(event):
                line.changed()

            fig.canvas.mpl_connect('draw_event', forceUpdate)
            plt.show()

        data = genfromtxt(args.plot[0], delimiter='\t')
        filter_1 = data[:,3]
        filter_2 = data[:,4]

        plot(data, filter_2, 'impossible gt dist')
        plot(data, filter_1, 'positions with >= 1 green')