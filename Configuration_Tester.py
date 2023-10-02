import os
import json
import time
import pandas as pd
from Functions import file_management
from Functions.Parsers import openacc_timing_data_parser

if __name__ == '__main__':
    # Read json file
    json_path = "JSONs/Config_Tester_Info.json"
    
    sod2d_info = None

    # Read opentuner user parameters.
    with open(json_path, 'r') as f:
        sod2d_info = json.load(f)
    
    # Set appropriate gpu id.
    os.environ["CUDA_VISIBLE_DEVICES"] = sod2d_info['gpu_ids']
    
    number_of_configs = sod2d_info["number_of_configs"]

    example_path = sod2d_info["example_path"]
    rank_num = sod2d_info["rank_num"]
    sod2d_path = sod2d_info["sod2d_path"]
    
    configurations_paths = sod2d_info["config_result_path"]
    configurations_dfs = []
    if len(configurations_paths) > 1:
        for path in configurations_paths:
            configurations_dfs.append(pd.read_csv(path))

    configurations = []
    for configuration_df in configurations_dfs:
        configurations.append(configuration_df.sort_values(by=['Mean Squared Calc Time']).head(number_of_configs))

    config_num = 0

    for convec_config in configurations[0]:
        os.system('export gang_num_fcijk=' + str(convec_config['gang_num_fcijk']))
        os.system('export vector_num_fcijk=' + str(convec_config['vector_num_fcijk']))
        os.system('export worker_num_fcijk=' + str(convec_config['worker_num_fcijk']))
        for diffu_config in configurations[1]:
            os.system('export gang_num_fdijk=' + str(convec_config['gang_num_fdijk']))
            os.system('export vector_num_fdijk=' + str(convec_config['vector_num_fdijk']))
            os.system('export worker_num_fdijk=' + str(convec_config['worker_num_fdijk']))
    
            # Execute Sod2d
            print("Executing Sod2d Application")
            execute_cmd = 'cd '+ example_path + ' &&'
            execute_cmd += 'time mpirun --allow-run-as-root --mca coll ^hcoll'
            execute_cmd += ' -np ' + rank_num + ' '
            execute_cmd += sod2d_path + '> openacc_timing.txt 2>&1'
            
            start_time = time.time()
            os.system(execute_cmd)
            end_time = time.time()
            
            results_folder = 'Archive/Whitebox_Analysis/Configurations_Check/Config_' + str(config_num)
            
            if not os.path.exists(results_folder):
                print("Creating folder: " + results_folder)
                # If it doesn't exist, create the folder and any intermediate directories
                os.makedirs(results_folder)
            
            time_data = {"time": end_time - start_time}
            json_filename = results_folder + '/time.json'
            
            with open(json_filename, "w") as json_file:
                json.dump(time_data, json_file)
            
            print("Moving results at: " + results_folder)
            # Move results to appropriate folder
            file_management.move_results(example_path, results_folder)

            print("Parsing OpenAcc Timing Analysis")
            # Parse OpenAcc Timing Analysis
            openacc_timing_data_parser.parser(results_folder)

            config_num += 1
