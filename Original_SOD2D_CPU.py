import os
import json
import time
from Functions import file_management_V2
        
def main():
    # Read JSON file
    json_path = 'JSONs/Original_SOD2D_CPU_Info.json'
    sod2d_info = None
    # Read opentuner user parameters.
    with open(json_path, 'r') as json_file:
        sod2d_info = json.load(json_file)

    # Set sod2d executable path.
    sod2d_path = sod2d_info['sod2d_path']
    # Set example file, folder path.
    example_path = sod2d_info['example_path']
    # Set folder to move results to.
    results_folder = sod2d_info['results_path']
    # Set number of proccesses to launch.
    rank_num = sod2d_info['rank_num']
    # Profiling metric.
    metric = sod2d_info['metric']

    os.environ['ACC_NUM_CORES'] = str(rank_num)
    os.environ['OMP_NUM_THREADS'] = str(rank_num)

    print('Executing Sod2d Application')
    execute_cd = f'cd {example_path} &&'
    # Create command for SOD2D execution with nsys profiling.
    execute_sod2d = execute_cd + f' vtune -collect {metric} -r ./{metric} mpirun'
    execute_sod2d = execute_sod2d + f' -np 1 --map-by node:PE={rank_num} {sod2d_path}'
    
    # Execute SOD2D.
    start = time.time()
    os.system(execute_sod2d)
    execution_time = time.time()-start
    os.system('cd ..')

    # Check if results folder exists.
    if not os.path.exists(results_folder):
        print('Creating folder:', results_folder)
        # If it doesn't exist, create the folder and any intermediate directories.
        os.makedirs(results_folder)

    # Move results to the appropriate folder.
    print('Moving results at:', results_folder)
    file_management_V2.move_results(example_path, results_folder)
        
    # Save execution time result
    json_filename = os.path.join(results_folder, 'time.json')
    with open(json_filename, 'w') as json_file:
        json.dump({'time': execution_time}, json_file)
    
if __name__ == '__main__':
    main()