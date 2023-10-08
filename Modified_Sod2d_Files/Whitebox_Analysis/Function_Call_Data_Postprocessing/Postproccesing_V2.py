import struct
import numpy as np
import pandas as pd
from alive_progress import alive_bar
from Array_Reader import data_read, load_array

def initialize_shape_variables(store_path: str, folder_id: int, function_name: str) -> list:
    """
    Initialize shape variables from binary data.

    Args:
        store_path (str): The path where data files are stored.
        folder_id (int): Folder identifier.
        function_name (str): Name of the function.

    Returns:
        list: List of shape variables.
    """
    variables = []
    for variable_name in ['nelem', 'npoin', 'nnode', 'ngaus', 'porder', 'ndime']:
        variables.append(struct.unpack('<iii', data_read(store_path, folder_id, function_name, variable_name))[1])
        
    return variables

def create_csvs(store_path: str, start_num: int, end_num: int, step: int) -> list:
    """
    Create CSV files for different functions and load data into them.

    Args:
        store_path (str): The path where data files are stored.
        start_num (int): Starting folder index.
        end_num (int): Ending folder index.
        step (int): Step size for folder indices.

    Returns:
        list: List of function names processed.
    """
    function_names = []

    # Initialize shape variables
    folder_idx = start_num
    function_name = 'full_convec_ijk'
    function_names.append(function_name)
    shape_variables = initialize_shape_variables(store_path, folder_idx, function_name)

    # Initialize array names
    fcijk_array_names = ['connec', 'Ngp', 'dNgp', 'He', 'xgp', 'dlxigp_ip', 'gpvol', 'atoIJK', 'invAtoIJK', 'gmshAtoI', 'gmshAtoJ', 'gmshAtoK', 'q', 'u', 'rho', 'pr', 'E', 'Rmass', 'Rmom', 'Rener']

    # Initialize dataframe
    full_convec_ijk_df = pd.DataFrame(columns=fcijk_array_names)
    full_convec_ijk_df.to_excel(store_path + '/full_convec_ijk.xlsx')

    print(f"Loading Data for {function_name} function and saving results to {store_path}/full_convec_ijk.xlsx")
    with alive_bar(int((end_num - start_num) / step)) as bar:
        # Load arrays and save to dataframe
        for folder_idx in range(start_num, end_num + 1, step):
            metric = []
            for array_name in fcijk_array_names:
                array = load_array(folder_idx, function_name, array_name, shape_variables)
                metric.append(np.median(array))

            info_dict = {array_name: value for array_name, value in zip(fcijk_array_names, metric)}
            info_df = pd.DataFrame([info_dict])
            full_convec_ijk_df = pd.concat([full_convec_ijk_df, info_df], ignore_index=True)
            full_convec_ijk_df.to_excel(store_path + '/full_convec_ijk.xlsx', index=False)
            bar()

    print(f'Finished loading data for {function_name} function.\n')

    # Initialize shape variables
    folder_idx = start_num
    function_name = 'full_diffu_ijk'
    function_names.append(function_name)
    shape_variables = initialize_shape_variables(store_path, folder_idx, function_name)

    # Initialize array names and variables
    fdijk_array_names = ['connec', 'Ngp', 'dNgp', 'He', 'xgp', 'dlxigp_ip', 'gpvol', 'atoIJK', 'invAtoIJK', 'gmshAtoI', 'gmshAtoJ', 'gmshAtoK', 'rho', 'u', 'Tem', 'mu_e', 'mu_sgs', 'Ml', 'mu_fluid', 'Rmass', 'Rmom', 'Rener']
    fdijk_var_names = ['Cp', 'Pr']

    # Initialize dataframe
    full_diffu_ijk_df = pd.DataFrame(columns=fdijk_array_names)
    full_diffu_ijk_df.to_excel(store_path + '/full_diffu_ijk.xlsx')

    print(f"Loading Data for {function_name} function and saving results to {store_path}/full_diffu_ijk.xlsx")
    with alive_bar(int((end_num - start_num) / step)) as bar:
        # Load arrays and save to dataframe
        for folder_idx in range(start_num, end_num + 1, step):
            metric = []
            for array_name in fdijk_array_names:
                array = load_array(folder_idx, function_name, array_name, shape_variables)
                metric.append(np.median(array))

            for var_name in fdijk_var_names:
                metric.append(struct.unpack('<fff', data_read(store_path, folder_idx, function_name, var_name))[1])

            info_dict = {array_name: value for array_name, value in zip(fdijk_array_names + fdijk_var_names, metric)}
            info_df = pd.DataFrame([info_dict])
            full_diffu_ijk_df = pd.concat([full_diffu_ijk_df, info_df], ignore_index=True)
            full_diffu_ijk_df.to_excel(store_path + '/full_diffu_ijk.xlsx', index=False)
            bar()
            
    print(f'Finished loading data for {function_name} function.')

    return function_names
