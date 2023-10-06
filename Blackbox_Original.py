import os
import json
import time
from Functions import file_management
from Functions.Parsers import openacc_timing_data_parser

if __name__ == '__main__':
    # Enviroment variable for OpenACC timing analysis.
    os.environ['PGI_ACC_TIME'] = '1'

    # Set id orded based on GPU BUS ID. 
    # Run nvidia-smi to see how gpus will be ordered and pick your poison.
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    
    # Read json file
    json_path = "JSONs/Original_Info.json"
    
    sod2d_info = None

    # Read opentuner user parameters.
    with open(json_path, 'r') as f:
        sod2d_info = json.load(f)
    
    # Set appropriate gpu id.
    os.environ["CUDA_VISIBLE_DEVICES"] = sod2d_info['gpu_ids']
    
    example_path = sod2d_info["example_path"]
    rank_num = sod2d_info["rank_num"]
    sod2d_path = sod2d_info["sod2d_path"]
    
    # Execute Sod2d
    print("Executing Sod2d Application")
    execute_cmd = 'cd '+ example_path + ' &&'
    execute_cmd += 'mpirun --allow-run-as-root --mca coll ^hcoll'
    execute_cmd += ' -np ' + rank_num + ' '
    execute_cmd += sod2d_path + '> openacc_timing.txt 2>&1'
    
    start_time = time.time()
    os.system(execute_cmd)
    end_time = time.time()
    
    results_folder = 'Archive/Blackbox_Analysis/Original_Example'
    
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
