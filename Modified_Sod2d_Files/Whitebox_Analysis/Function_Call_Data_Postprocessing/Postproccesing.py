import struct
import numpy as np
import pandas as pd
from alive_progress import alive_bar
from Array_Reader import data_read, load_array

def initialize_shape_variables(store_path: str, folder_id: int, function_name: str) -> list:
    variables = []
    for variable_name in ['nelem', 'npoin', 'nnode', 'ngaus', 'porder', 'ndime']:
        variables.append(struct.unpack('<iii', data_read(store_path, folder_id, function_name, variable_name))[1])
        
    return variables

def create_csvs(store_path: str, start_num: int, end_num: int, step: int):

    function_names = []
# ------------------------ Full Element Convection ------------------------

    # Initialize shape variables
    folder_idx = start_num
    function_name = 'full_convec_ijk'
    function_names.append(function_name)
    shape_variables = initialize_shape_variables(store_path, folder_idx, function_name)

    # Initialize array names
    fcijk_array_names = ['connec', 'Ngp', 'dNgp', 'He', 'xgp', 'dlxigp_ip', 'gpvol', 'atoIJK', 'invAtoIJK', 'gmshAtoI', 'gmshAtoJ', 'gmshAtoK', 'q', 'u', 'rho', 'pr', 'E', 'Rmass', 'Rmom', 'Rener']

    # Initialize dataframe
    full_convec_ijk_df = pd.DataFrame(columns = fcijk_array_names)
    full_convec_ijk_df.to_csv(store_path + '/full_convec_ijk.csv')

    print(f"Loading Data for {function_name} function and saving results to {store_path}/full_convec_ijk.csv")
    with alive_bar(int(end_num/step)) as bar:
        # Load arrays and save to dataframe
        for folder_idx in range(0, end_num+1, step):
            metric = []
            for array_name in fcijk_array_names:
                array = load_array(folder_idx, function_name, array_name, shape_variables)
                metric.append(np.median(array))

            info_dict = {
                'connec': metric[0], 
                'Ngp': metric[1], 
                'dNgp': metric[2], 
                'He': metric[3], 
                'xgp': metric[4], 
                'dlxigp_ip': metric[5], 
                'gpvol': metric[6], 
                'atoIJK': metric[7], 
                'invAtoIJK': metric[8], 
                'gmshAtoI': metric[9], 
                'gmshAtoJ': metric[10], 
                'gmshAtoK': metric[11], 
                'q': metric[12], 
                'u': metric[13], 
                'rho': metric[14], 
                'pr': metric[15],
                'E': metric[16], 
                'Rmass': metric[17], 
                'Rmom': metric[18], 
                'Rener': metric[19]
            }

            info_df = pd.DataFrame([info_dict])
            full_convec_ijk_df = pd.concat([full_convec_ijk_df, info_df], ignore_index=True)
            full_convec_ijk_df.to_csv(store_path + '/full_convec_ijk.csv', index=False)
            bar()

    print(f'Finished loading data for {function_name} function.\n')

    # ------------------------ Full Element Diffusion ------------------------

    # Initialize shape variables
    folder_idx = start_num
    function_name = 'full_diffu_ijk'
    function_names.append(function_name)
    shape_variables = initialize_shape_variables(store_path, folder_idx, function_name)

    # Initialize array names and variables
    fdijk_array_names = ['connec', 'Ngp', 'dNgp', 'He', 'xgp', 'dlxigp_ip', 'gpvol', 'atoIJK', 'invAtoIJK', 'gmshAtoI', 'gmshAtoJ', 'gmshAtoK', 'rho', 'u', 'Tem', 'mu_e', 'mu_sgs', 'Ml', 'mu_fluid', 'Rmass', 'Rmom', 'Rener']
    fdijk_var_names = ['Cp', 'Pr']

    # Initialize dataframe
    full_diffu_ijk_df = pd.DataFrame(columns = fdijk_array_names)
    full_diffu_ijk_df.to_csv(store_path + '/full_diffu_ijk.csv')

    print(f"Loading Data for {function_name} function and saving results to {store_path}/full_diffu_ijk.csv")
    with alive_bar(int(end_num/step)) as bar:
        # Load arrays and save to dataframe
        for folder_idx in range(0, end_num+1, step):
            metric = []
            for array_name in fdijk_array_names:
                array = load_array(folder_idx, function_name, array_name, shape_variables)
                metric.append(np.median(array))
            
            for var_name in fdijk_var_names:
                metric.append(struct.unpack('<fff', data_read(store_path, folder_idx, function_name, var_name))[1])
            
            info_dict = {
                'connec': metric[0], 
                'Ngp': metric[1], 
                'dNgp': metric[2], 
                'He': metric[3], 
                'xgp': metric[4], 
                'dlxigp_ip': metric[5], 
                'gpvol': metric[6], 
                'atoIJK': metric[7], 
                'invAtoIJK': metric[8], 
                'gmshAtoI': metric[9], 
                'gmshAtoJ': metric[10], 
                'gmshAtoK': metric[11], 
                'rho': metric[12], 
                'u': metric[13],
                'Tem': metric[14], 
                'mu_e': metric[15], 
                'mu_sgs': metric[16],
                'Ml': metric[17], 
                'mu_fluid': metric[18], 
                'Rmass': metric[19], 
                'Rmom': metric[20], 
                'Rener': metric[21],
                'Cp': metric[22], 
                'Pr': metric[23]
            }

            info_df = pd.DataFrame([info_dict])
            full_diffu_ijk_df = pd.concat([full_diffu_ijk_df, info_df], ignore_index=True)
            full_diffu_ijk_df.to_csv(store_path + '/full_diffu_ijk.csv', index=False)
            bar()
            
    print(f'Finished loading data for {function_name} function.')

    return function_names