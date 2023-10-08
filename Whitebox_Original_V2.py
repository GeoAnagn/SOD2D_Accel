import os
import json
import time
import numpy as np
import pandas as pd
from Functions.Parsers import openacc_timing_data_parser_V2
from Functions import file_management_V2, error_calculation_V2, array_reader_V2

if __name__ == '__main__':
    # Enable OpenACC timing analysis.
    os.environ['PGI_ACC_TIME'] = '1'

    # Set GPU ordering based on PCI BUS ID.
    os.environ['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'
    
    # Read JSON file for function information.
    json_path = 'JSONs/Whitebox_Info.json'
    function_info = None

    with open(json_path, 'r') as f:
        function_info = json.load(f)
    
    # Set the desired GPU.
    os.environ['CUDA_VISIBLE_DEVICES'] = function_info['gpu_ids']

    # Get file paths.
    paths = file_management_V2.path_definitions()
    example_path = paths[0]

    data_path = function_info['data_path']
    function_path = function_info['func_path']
    function_version = function_info['func_ver']
    function_executable = function_info['func_exec']
    function_name = function_info['func_name']
    function_call_indexes = function_info['func_call_idx']
    rank_num = function_info['rank_num']

    # Create a list to store function call information.
    func_call_info_list = []

    for function_call_index in function_call_indexes:
        # Get the folder path to load data from.
        data_folder_path = f'{data_path}/Data_{str(function_call_index)}/{function_name}'
        os.environ['folder_path'] = data_folder_path

        # Execute the function.
        print(f'Executing {function_name} function for call index {function_call_index}')
        execute_cmd = f'cd {example_path} && '
        execute_cmd += 'mpirun --allow-run-as-root --mca coll ^hcoll '
        execute_cmd += f'-np {rank_num} '
        execute_cmd += f'{function_path}{function_version}{function_executable} > openacc_timing.txt 2>&1'
        start_time = time.time()
        os.system(execute_cmd)
        end_time = time.time()
        os.system('cd ..')

        # Parse OpenACC Timing Analysis.
        print('Parsing OpenACC Timing Analysis')
        openacc_timing_data_parser_V2.parser(example_path)
        openacc_df = pd.read_excel(f'{example_path}/openacc_timing.xlsx')
        calculation_time = openacc_df['Compute Time'].iloc[0]

        # Calculate errors.
        print('Calculating Errors')
        shape_variables = array_reader_V2.initialize_shape_variables(data_folder_path)
        Rmass_error = error_calculation_V2.error_calc(data_folder_path, example_path, 'Rmass', shape_variables)
        Rmom_error = error_calculation_V2.error_calc(data_folder_path, example_path, 'Rmom', shape_variables)
        Rener_error = error_calculation_V2.error_calc(data_folder_path, example_path, 'Rener', shape_variables)
        
        # Create a JSON file with function call information.
        call_data = {'Total Time': end_time - start_time,
                     'Calculation Time': calculation_time,
                     'Rmass Error': Rmass_error,
                     'Rmom Error': Rmom_error,
                     'Rener Error': Rener_error}
        json_filename = example_path + '/function_call.json'
        with open(json_filename, 'w') as json_file:
            json.dump(call_data, json_file)

        # Create a folder to store results.
        results_folder = f'Archive/Whitebox/Original_Results/{function_name}/Call_{str(function_call_index)}'
        if not os.path.exists(results_folder):
            print('Creating folder: ' + results_folder)
            os.makedirs(results_folder)

        # Move results to the appropriate folder.
        print('Moving results at: ' + results_folder)
        file_management_V2.move_results(example_path, results_folder)

        func_call_info_list.append(call_data)

    # Create an Excel file with all function call information.
    func_call_info_df = pd.DataFrame(func_call_info_list)
    func_call_info_df.to_excel(f'Archive/Whitebox/Original_Results/{function_name}/function_call_info.xlsx', index=False)

    # Calculate Mean Squared Calculation Time and Errors.
    func_call_info_df['Squared Calc Time'] = np.power((func_call_info_df['Calculation Time']), 2)
    func_call_info_df['Squared Rmass Error'] = np.power((func_call_info_df['Rmass Error']), 2)
    func_call_info_df['Squared Rmom Error'] = np.power((func_call_info_df['Rmom Error']), 2)
    func_call_info_df['Squared Rener Error'] = np.power((func_call_info_df['Rener Error']), 2)

    mean_squared_time = func_call_info_df['Squared Calc Time'].mean()
    mean_squared_Rmass_Error = func_call_info_df['Squared Rmass Error'].mean()
    mean_squared_Rmom_Error = func_call_info_df['Squared Rmom Error'].mean()
    mean_squared_Rener_Error = func_call_info_df['Squared Rener Error'].mean()
    
    result = {
        "Mean Squared Calc Time": mean_squared_time,
        "Mean Squared Rmass Error": mean_squared_Rmass_Error,
        "Mean Squared Rmom Error": mean_squared_Rmom_Error,
        "Mean Squared Rener Error": mean_squared_Rener_Error
    }

    # Save results to a JSON file.
    print(f'Saving results file at Archive/Whitebox/Original_Results/{function_name}/results.json')
    with open(f'Archive/Whitebox/Original_Results/{function_name}/results.json', 'w') as results:
        json.dump(result, results)
