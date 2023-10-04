import os
import json
from Modified_Sod2d_Files.Whitebox_Analysis.Function_Call_Data_Postprocessing import Postproccesing, Clustering

if __name__ == '__main__':
    # Set GPU ordering based on PCI_BUS_ID
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"

    # Load config file
    with open('JSONs/Whitebox_Data.json') as config_file:
        config = json.load(config_file)

    # Load config for executing
    sod2d_path = config['sod2d_path']
    example_path = config['example_path']
    rank_num = config['rank_num']

    # Set the desired GPU(s) based on configuration
    os.environ["CUDA_VISIBLE_DEVICES"] = config['gpu_ids']

    # Load config for storing
    store_path = config['store_path']
    start_num = config['start_num']
    end_num = config['end_num']
    step = config['step']
    cluster_num = config['cluster_num']

    os.environ['store_path'] = store_path
    os.environ['file_num'] = str(start_num)
    os.environ['max_calls'] = str(end_num)
    os.environ['calls_step'] = str(step)

    # Execute Sod2d
    print("Executing Sod2d Application")
    execute_cmd = 'cd '+ example_path + ' &&'
    execute_cmd += ' mpirun --allow-run-as-root --mca coll ^hcoll'
    execute_cmd += ' -np ' + rank_num + ' '
    execute_cmd += sod2d_path 
    os.system(execute_cmd)

    print('Analysing Data Collected and Creating CSVs...')
    function_names = Postproccesing.create_csvs(store_path, start_num, end_num, step)

    print('Finished creating CSVs.')

    print('Clustering calls...')
    Clustering.create_clusters(store_path, function_names, cluster_num)

    print('Found all representative calls.')