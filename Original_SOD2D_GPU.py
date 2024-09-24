import os
import json
import time
import signal
import subprocess
from Functions import file_management_V2
from Functions.Parsers import nsys_sqlite_parser, utilization_parser
        
def main():
    # Read JSON file
    json_path = 'JSONs/Original_SOD2D_GPU_Info.json'
    sod2d_info = None
    # Read opentuner user parameters.
    with open(json_path, 'r') as json_file:
        sod2d_info = json.load(json_file)

    # Set id ordered based on GPU BUS ID.
    # Run nvidia-smi to see how GPUs will be ordered and pick your poison.
    os.environ['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'
    # Set appropriate GPU id.
    os.environ['CUDA_VISIBLE_DEVICES'] = sod2d_info['gpu_ids']

    # Set sod2d executable path.
    sod2d_path = sod2d_info['sod2d_path']
    # Set example file, folder path.
    example_path = sod2d_info['example_path']
    # Set folder to move results to.
    results_folder = sod2d_info['results_path']
    # Set number of proccesses to launch.
    rank_num = sod2d_info['rank_num']

    print('Executing Sod2d Application')
    execute_cd = f'cd {example_path} &&'
    # Create command for SOD2D execution with nsys profiling.
    execute_sod2d = execute_cd + ' nsys profile -t cuda,openacc,nvtx --cuda-memory-usage=true --stats=true mpirun --allow-run-as-root -x UCX_NET_DEVICES=mlx5_0:1 -x HCOLL_MAIN_IB=mlx5_0:1'
    execute_sod2d += f' -np {rank_num} {sod2d_path}'
    # Create command for SOD2D utilization capturing.
    execute_monitor = execute_cd + ' nvidia-smi dmon -s pucvmet > Utilization.csv'
    
    # Open a process and start capturing utilization.
    pro = subprocess.Popen(execute_monitor, shell=True, preexec_fn=os.setsid)
    start = time.time()
    # Execute SOD2D.
    os.system(execute_sod2d)
    end = time.time()
    print('Execution Time:', end-start)
    # Stop utilization capturing process.
    os.killpg(os.getpgid(pro.pid), signal.SIGTERM)
    os.system('cd ..')

    # Check if results folder exists.
    if not os.path.exists(results_folder):
        print('Creating folder:', results_folder)
        # If it doesn't exist, create the folder and any intermediate directories.
        os.makedirs(results_folder)

    # Move results to the appropriate folder.
    print('Moving results at:', results_folder)
    file_management_V2.move_results(example_path, results_folder)

    # # Parse Nsight Timing Analysis.
    print('Parsing Nsight Timing Analysis')
    nsys_sqlite_parser.parser(results_folder)

    # # Parse Utilization Analysis.
    # print('Parsing Utilization Analysis')
    # utilization_parser.parser(results_folder)

if __name__ == '__main__':
    main()
