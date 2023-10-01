import struct
import numpy as np

def data_read(data_folder: str, variable_name: str):
    with open(data_folder + '/' + variable_name + '.bin', 'rb') as data_file:
        return data_file.read()

def array_type(arr_shape: tuple, arr_type: str):
    # Float Type Array with Byte Size 4
    if arr_type == 'f':
       array = np.array(np.zeros(shape = arr_shape, dtype = float))
       byte_size = 4
       
    # Integer Type Array with Byte Size 4
    elif arr_type == 'i':
       array = np.array(np.zeros(shape = arr_shape, dtype = int))
       byte_size = 4
    
    return array, byte_size

def read_array(x: int, val_type: str, array_data):
     # Create empty array of shape -> (x) with appropriate type.
    array, byte_size = array_type((x), val_type)
    
    # Set data index
    data_index = 4
    
    # Unpack and store data on created array.
    for i in range(0, x):
        array[i] = struct.unpack('<' + val_type, array_data[data_index : data_index + byte_size])[0]
        data_index += byte_size
                
    return array

def read_2D_array(x: int, y: int, val_type: str, array_data) -> np.ndarray:
    # Create empty array of shape -> (x,y) with appropriate type.
    array, byte_size = array_type((x, y), val_type)
    
    # Set data index
    data_index = 4
    # Unpack and store data on created array.
    for i in range(0, y):
        for j in range(0, x):
            array[j][i] = struct.unpack('<' + val_type, array_data[data_index : data_index + byte_size])[0]
            data_index += byte_size
                
    return array

def initialize_shape_variables(data_folder: str) -> list:
    variables = []
    for variable_name in ['nelem', 'npoin', 'nnode', 'ngaus', 'porder', 'ndime']:
        variables.append(struct.unpack('<iii', data_read(data_folder, variable_name))[1])
        
    return variables 

def array_decider(variable_name: str, shape_variables: list, array_data):
    nelem = shape_variables[0]
    npoin = shape_variables[1]
    nnode = shape_variables[2]
    ngaus = shape_variables[3]
    porder = shape_variables[4]
    ndime = shape_variables[5]

    # Load Rmass(npoin) array of function full_convec_ijk or full_diffusion_ijk.
    if variable_name == 'Rmass':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    
    # Load Rmom(npoin,ndime) array of function full_convec_ijk or full_diffusion_ijk.
    elif variable_name == 'Rmom':
        array = read_2D_array(x=npoin, y=ndime, val_type='f', array_data=array_data)
    
    # Load Rener(npoin) array of function full_convec_ijk or full_diffusion_ijk.
    elif variable_name == 'Rener':
        array = read_array(x=npoin, val_type='f', array_data=array_data)

    return array
