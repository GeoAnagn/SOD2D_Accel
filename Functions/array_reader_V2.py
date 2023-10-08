import struct
import numpy as np

# Function to read binary data from a file
def data_read(data_folder: str, variable_name: str):
    file_path = f'{data_folder}/{variable_name}.bin'
    with open(file_path, 'rb') as data_file:
        return data_file.read()

# Function to read a 1D array from binary data
def read_array(data, shape, dtype, fortran_order=True):
    # Create a NumPy array from binary data
    array = np.frombuffer(data, dtype=dtype, offset=4)
    # Reshape the array using Fortran order (column-major) or C order (row-major)
    return array.reshape(shape, order='F' if fortran_order else 'C')

# Function to initialize shape variables from binary data
def initialize_shape_variables(data_folder: str):
    # List of variable names for shape variables
    shape_variables = ['nelem', 'npoin', 'nnode', 'ngaus', 'porder', 'ndime']
    variables = []
    
    for variable_name in shape_variables:
        # Read binary data and unpack the second integer (Fortran uses 1-based indexing)
        data = data_read(data_folder, variable_name)
        value = struct.unpack('<iii', data)[1]
        variables.append(value)
        
    return variables 

# Function to decide which array to load based on variable name
def array_decider(variable_name: str, shape_variables: list, array_data):
    npoin, ndime = shape_variables[1], shape_variables[5]

    if variable_name == 'Rmass':
        # Load Rmass(npoin) array and reshape it
        return read_array(array_data, (npoin,), 'f')
    
    elif variable_name == 'Rmom':
        # Load Rmom(npoin, ndime) array and reshape it
        return read_array(array_data, (npoin, ndime), 'f')
    
    elif variable_name == 'Rener':
        # Load Rener(npoin) array and reshape it
        return read_array(array_data, (npoin,), 'f')

    return None  # Handle unsupported variable names gracefully
