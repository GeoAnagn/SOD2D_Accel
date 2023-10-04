import os
import json
import time
import pandas as pd
from Functions import file_management
from Functions.Parsers import openacc_timing_data_parser

if __name__ == '__main__':
    # Enable OpenACC timing analysis environment variable
    os.environ['PGI_ACC_TIME'] = '1'

    # Set GPU ordering based on PCI_BUS_ID
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"

    # Read configuration from JSON file
    json_path = "JSONs/Config_Tester_Info.json"
    sod2d_info = None
    with open(json_path, 'r') as f:
        sod2d_info = json.load(f)
    
    # Set the desired GPU(s) based on configuration
    os.environ["CUDA_VISIBLE_DEVICES"] = sod2d_info['gpu_ids']
    
    # Extract configuration values
    number_of_configs = sod2d_info["number_of_configs"]
    example_path = sod2d_info["example_path"]
    rank_num = sod2d_info["rank_num"]
    sod2d_path = sod2d_info["sod2d_path"]
    configurations_paths = sod2d_info["config_result_path"]
    
    # Load configuration DataFrames
    configurations_dfs = []
    if len(configurations_paths) > 1:
        for path in configurations_paths:
            configurations_dfs.append(pd.read_csv(path))

    # Select top configurations
    configurations = []
    for configuration_df in configurations_dfs:
        configurations.append(configuration_df.sort_values(by=['Mean Squared Calc Time']).head(number_of_configs))

    config_num = 0

    # Create a DataFrame to store results
    results_df = pd.DataFrame(columns=['time', 'gang_num_fcijk', 'worker_num_fcijk', 'vector_num_fcijk', 'gang_num_fdijk', 'worker_num_fdijk', 'vector_num_fdijk'])
    # Iterate through combinations of configurations
    for i, convec_config in configurations[0].iterrows():
        gang_num_fcijk = str(int(convec_config['gang_num_fcijk']))
        worker_num_fcijk = str(int(convec_config['worker_num_fcijk']))
        vector_num_fcijk = str(int(convec_config['vector_num_fcijk']))
        os.environ["gang_num_fcijk"] = gang_num_fcijk
        os.environ['worker_num_fcijk'] = worker_num_fcijk
        os.environ['vector_num_fcijk'] = vector_num_fcijk
        for j, diffu_config in configurations[1].iterrows():
            gang_num_fdijk = str(int(diffu_config['gang_num_fdijk']))
            worker_num_fdijk = str(int(diffu_config['worker_num_fdijk']))
            vector_num_fdijk = str(int(diffu_config['vector_num_fdijk']))
            os.environ["gang_num_fdijk"] = gang_num_fdijk
            os.environ['worker_num_fdijk'] = worker_num_fdijk
            os.environ['vector_num_fdijk'] = vector_num_fdijk

            # Execute Sod2d
            print("Executing Sod2d Application")
            execute_cmd = 'cd '+ example_path + ' &&'
            execute_cmd += ' mpirun --allow-run-as-root --mca coll ^hcoll'
            execute_cmd += ' -np ' + rank_num + ' '
            execute_cmd += sod2d_path + '> openacc_timing.txt 2>&1'

            start_time = time.time()
            os.system(execute_cmd)
            end_time = time.time()
            total_time = end_time - start_time
            results_folder = 'Archive/Whitebox_Analysis/Configurations_Check/Config_' + str(config_num)
            
            # Create results folder if it doesn't exist
            if not os.path.exists(results_folder):
                print("Creating folder: " + results_folder)
                os.makedirs(results_folder)
            
            # Store configuration data as JSON
            config_data = {"time": total_time,
                         "gang_num_fcijk": gang_num_fcijk,
                         "worker_num_fcijk": worker_num_fcijk,
                         "vector_num_fcijk": vector_num_fcijk,
                         "gang_num_fdijk": gang_num_fdijk,
                         "worker_num_fdijk": worker_num_fdijk,
                         "vector_num_fdijk": vector_num_fdijk}
            
            json_filename = results_folder + '/data.json'
            
            with open(json_filename, "w") as json_file:
                json.dump(config_data, json_file)
            
            # Append configuration data to results DataFrame and save it
            results_df.loc[config_num] = config_data
            results_df.to_csv("Archive/Whitebox_Analysis/Configurations_Check/results.csv")
            
            print("Moving results at: " + results_folder)
            # Move results to appropriate folder
            file_management.move_results(example_path, results_folder)

            # Remove .h5 files
            file_management.rm_files(results_folder)

            print("Parsing OpenAcc Timing Analysis")
            # Parse OpenAcc Timing Analysis
            openacc_timing_data_parser.parser(results_folder)

            config_num += 1
