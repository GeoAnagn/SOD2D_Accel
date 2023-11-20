import os
import json
import time
from Functions import file_management_V2
from Functions.Parsers import openacc_timing_data_parser
import subprocess
import signal
# import multiprocessing
# from multiprocessing import Manager

# def run_example(thread_id, example_path, rank_num, sod2d_path, results_folder):
#     if thread_id == 0:
#         execute_cmd += ' nvidia-smi dmon -s pucvmet > Utilization.csv'
#         os.system(execute_cmd)
#     elif thread_id == 1:
        
def main():
    # Environment variable for OpenACC timing analysis.
    os.environ['PGI_ACC_TIME'] = '1'

    # Set id ordered based on GPU BUS ID.
    # Run nvidia-smi to see how GPUs will be ordered and pick your poison.
    os.environ['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'

    # Read JSON file
    json_path = 'JSONs/Original_Info.json'

    sod2d_info = None

    # Read opentuner user parameters.
    with open(json_path, 'r') as json_file:
        sod2d_info = json.load(json_file)

    # Set appropriate GPU id.
    os.environ['CUDA_VISIBLE_DEVICES'] = sod2d_info['gpu_ids']

    example_path = sod2d_info['example_path']
    rank_num = sod2d_info['rank_num']
    sod2d_path = sod2d_info['sod2d_path']
    results_folder = 'Archive/Blackbox_Analysis/Original_Example'

    # Execute Sod2d
    print('Executing Sod2d Application')
    execute_cmd = f'cd {example_path} &&'
    execute_monitor = execute_cmd + ' nvidia-smi dmon -s pucvmet > Utilization.csv'
    execute_cmd += ' mpirun --allow-run-as-root --mca coll ^hcoll'
    execute_cmd += f' -np {rank_num} {sod2d_path} > openacc_timing.txt 2>&1'

    pro = subprocess.Popen(execute_monitor, shell=True, preexec_fn=os.setsid)
    start_time = time.time()
    os.system(execute_cmd)
    end_time = time.time()
    os.killpg(os.getpgid(pro.pid), signal.SIGTERM)
    os.system('cd ..')

    if not os.path.exists(results_folder):
        print('Creating folder:', results_folder)
        # If it doesn't exist, create the folder and any intermediate directories
        os.makedirs(results_folder)

    time_data = {'time': end_time - start_time}
    json_filename = os.path.join(results_folder, 'time.json')

    with open(json_filename, 'w') as json_file:
        json.dump(time_data, json_file)

    print('Moving results at:', results_folder)
    # Move results to the appropriate folder
    file_management_V2.move_results(example_path, results_folder)

    print('Parsing OpenAcc Timing Analysis')
    # Parse OpenAcc Timing Analysis
    openacc_timing_data_parser.parser(results_folder)


    # manager = Manager()
    # stop_monitor = manager.list()
    # thread_id = [0, 1]
    # with multiprocessing.Pool(processes=2) as pool:
    #     done = pool.starmap(run_example , [(thread_id[i], example_path, rank_num, sod2d_path, results_folder, stop_monitor) for i in range(len(thread_id))])
        
if __name__ == '__main__':
    main()
