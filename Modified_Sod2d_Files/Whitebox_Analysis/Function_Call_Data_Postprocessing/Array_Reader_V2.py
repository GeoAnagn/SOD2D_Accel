import struct
import numpy as np

# Define constants for byte sizes
BYTE_SIZE_FLOAT = 4
BYTE_SIZE_INT = 4

# Define a dictionary for array types and byte sizes
ARRAY_TYPES = {
    'f': (np.float32, BYTE_SIZE_FLOAT),
    'i': (np.int32, BYTE_SIZE_INT)
}

def data_read(path: str, folder_id: int, function_name: str, variable_name: str):
    """
    Read binary data from a file and return it as bytes.

    Args:
        path (str): The base path to the data files.
        folder_id (int): Folder identifier.
        function_name (str): Name of the function.
        variable_name (str): Name of the variable.

    Returns:
        bytes: The binary data read from the file.
    """
    file_path = f'{path}/Data_{folder_id}/{function_name}/{variable_name}.bin'
    with open(file_path, 'rb') as data_file:
        return data_file.read()

def array_type(arr_shape: tuple, arr_type: str):
    """
    Create an empty NumPy array with the specified shape and data type.

    Args:
        arr_shape (tuple): The shape of the array.
        arr_type (str): The data type ('f' for float, 'i' for integer).

    Returns:
        numpy.ndarray: An empty NumPy array.
        int: The byte size of each element in the array.
    """
    if arr_type in ARRAY_TYPES:
        dtype, byte_size = ARRAY_TYPES[arr_type]
        array = np.empty(shape=arr_shape, dtype=dtype)
        return array, byte_size
    else:
        raise ValueError(f"Unsupported data type: {arr_type}")

def read_array(x: int, val_type: str, array_data):
    """
    Read a 1D array from binary data.

    Args:
        x (int): Number of elements in the array.
        val_type (str): Data type ('f' for float, 'i' for integer).
        array_data (bytes): The binary data to read from.

    Returns:
        numpy.ndarray: A 1D NumPy array containing the read data.
    """
    array, byte_size = array_type((x,), val_type)
    data_index = 4
    flat_array = np.frombuffer(array_data, dtype=np.dtype('<' + val_type), offset=data_index)
    array = flat_array.reshape((x,))
    return array

def read_2D_array(x: int, y: int, val_type: str, array_data):
    """
    Read a 2D array from binary data.

    Args:
        x (int): Number of rows in the array.
        y (int): Number of columns in the array.
        val_type (str): Data type ('f' for float, 'i' for integer).
        array_data (bytes): The binary data to read from.

    Returns:
        numpy.ndarray: A 2D NumPy array containing the read data.
    """
    array, byte_size = array_type((x, y), val_type)
    data_index = 4
    flat_array = np.frombuffer(array_data, dtype=np.dtype('<' + val_type), offset=data_index)
    array = flat_array.reshape((y, x)).T
    return array

def read_3D_array(x: int, y: int, z: int, val_type: str, array_data):
    """
    Read a 3D array from binary data.

    Args:
        x (int): Number of elements along the X-axis.
        y (int): Number of elements along the Y-axis.
        z (int): Number of elements along the Z-axis.
        val_type (str): Data type ('f' for float, 'i' for integer).
        array_data (bytes): The binary data to read from.

    Returns:
        numpy.ndarray: A 3D NumPy array containing the read data.
    """
    array, byte_size = array_type((x, y, z), val_type)
    data_index = 4
    flat_array = np.frombuffer(array_data, dtype=np.dtype('<' + val_type), offset=data_index)
    array = flat_array.reshape((z, y, x)).transpose(2, 1, 0)
    return array

def read_4D_array(x: int, y: int, z: int, n: int, val_type: str, array_data):
    """
    Read a 4D array from binary data.

    Args:
        x (int): Number of elements along the X-axis.
        y (int): Number of elements along the Y-axis.
        z (int): Number of elements along the Z-axis.
        n (int): Number of elements along the N-axis.
        val_type (str): Data type ('f' for float, 'i' for integer).
        array_data (bytes): The binary data to read from.

    Returns:
        numpy.ndarray: A 4D NumPy array containing the read data.
    """
    array, byte_size = array_type((x, y, z, n), val_type)
    data_index = 4
    flat_array = np.frombuffer(array_data, dtype=np.dtype('<' + val_type), offset=data_index)
    array = flat_array.reshape((n, z, y, x)).transpose(3, 2, 1, 0)
    return array

def load_array(folder_id: int, function_name: str, array_name, shape_variables: list):
    """
    Load an array from binary data based on its name and dimensions.

    Args:
        folder_id (int): Folder identifier.
        function_name (str): Name of the function.
        array_name (str): Name of the array.
        shape_variables (list): List of shape variables.

    Returns:
        numpy.ndarray: The loaded array.
    """
    nelem, npoin, nnode, ngaus, porder, ndime = shape_variables

    array_data = data_read(folder_id, function_name, array_name)

    if array_name == 'connec':
        array = read_2D_array(x=nelem, y=nnode, val_type='f', array_data=array_data)
    elif array_name == 'Ngp':
        array = read_2D_array(x=ngaus, y=nnode, val_type='f', array_data=array_data)
    elif array_name == 'dNgp':
        array = read_3D_array(x=ndime, y=nnode, z=ngaus, val_type='f', array_data=array_data)
    elif array_name == 'He':
        array = read_4D_array(x=ndime, y=ndime, z=ngaus, n=nelem, val_type='f', array_data=array_data)
    elif array_name == 'xgp':
        array = read_2D_array(x=ngaus, y=ndime, val_type='f', array_data=array_data)
    elif array_name == 'dlxigp_ip':
        array = read_3D_array(x=ngaus, y=ndime, z=porder+1, val_type='f', array_data=array_data)
    elif array_name == 'gpvol':
        array = read_3D_array(x=1, y=ngaus, z=nelem, val_type='f', array_data=array_data)
    elif array_name == 'atoIJK':
        array = read_array(x=nnode, val_type='i', array_data=array_data)
    elif array_name == 'invAtoIJK':
        array = read_3D_array(x=porder+1, y=porder+1, z=porder+1, val_type='i', array_data=array_data)
    elif array_name == 'gmshAtoI':
        array = read_array(x=nnode, val_type='i', array_data=array_data)
    elif array_name == 'gmshAtoJ':
        array = read_array(x=nnode, val_type='i', array_data=array_data)
    elif array_name == 'gmshAtoK':
        array = read_array(x=nnode, val_type='i', array_data=array_data)
    elif array_name == 'q':
        array = read_2D_array(x=npoin, y=ndime, val_type='f', array_data=array_data)
    elif array_name == 'u':
        array = read_2D_array(x=npoin, y=ndime, val_type='f', array_data=array_data)
    elif array_name == 'rho':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    elif array_name == 'pr':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    elif array_name == 'E':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    elif array_name == 'Rmass':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    elif array_name == 'Rmom':
        array = read_2D_array(x=npoin, y=ndime, val_type='f', array_data=array_data)
    elif array_name == 'Rener':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    elif array_name == 'Tem':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    elif array_name == 'mu_e':
        array = read_2D_array(x=nelem, y=ngaus, val_type='f', array_data=array_data)
    elif array_name == 'mu_sgs':
        array = read_2D_array(x=nelem, y=ngaus, val_type='f', array_data=array_data)
    elif array_name == 'Ml':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    elif array_name == 'mu_fluid':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    else:
        raise ValueError(f"Unsupported array name: {array_name}")

    return array
