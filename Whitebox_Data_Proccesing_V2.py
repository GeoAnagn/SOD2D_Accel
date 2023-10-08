import os
import json
from Modified_Sod2d_Files.Whitebox_Analysis.Function_Call_Data_Postprocessing import Postproccesing_V2, Clustering_V2

if __name__ == '__main__':
    # Set GPU ordering based on PCI_BUS_ID
    os.environ['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'

    # Load configuration from the JSON file
    with open('JSONs/Whitebox_Data.json') as config_file:
        config = json.load(config_file)

    # Load configuration for executing Sod2d
    sod2d_path = config['sod2d_path']  # Path to Sod2d executable
    example_path = config['example_path']  # Path to example directory
    rank_num = config['rank_num']  # Number of MPI ranks to use

    # Set the desired GPU(s) based on the configuration
    os.environ['CUDA_VISIBLE_DEVICES'] = config['gpu_ids']

    # Load configuration for storing data
    store_path = config['store_path']  # Path to store data
    start_num = config['start_num']  # Start number for data storage
    end_num = config['end_num']  # End number for data storage
    step = config['step']  # Step size for data storage
    cluster_num = config['cluster_num']  # Number of clusters to create

    # Set environment variables for data storage
    os.environ['store_path'] = store_path
    os.environ['file_num'] = str(start_num)
    os.environ['max_calls'] = str(end_num)
    os.environ['calls_step'] = str(step)

    # Execute Sod2d
    print('Executing Sod2d Application')
    execute_cmd = f'cd {example_path} && '
    execute_cmd += 'mpirun --allow-run-as-root --mca coll ^hcoll'
    execute_cmd += f' -np {rank_num} '
    execute_cmd += sod2d_path
    os.system(execute_cmd)
    os.system('cd ..')

    # Analyze Data Collected and Create CSVs
    print('Analyzing Data Collected and Creating CSVs...')
    function_names = Postproccesing_V2.create_csvs(store_path, start_num, end_num, step)

    print('Finished creating CSVs.')

    # Cluster calls
    print('Clustering calls...')
    Clustering_V2.create_clusters(store_path, function_names, cluster_num)

    print('Found all representative calls.')
