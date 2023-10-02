import os
import math
import struct
import numpy as np
import pandas as pd

def data_read(path: str, folder_id: int, function_name: str, variable_name: str):
    with open(path + '/data_' + str(folder_id) + '/' + function_name + '/' + variable_name + '.bin', 'rb') as data_file:
        return data_file.read()
    
def initialize_shape_variables(folder_id: int, function_name: str) -> list:
    variables = []
    for variable_name in ['nelem', 'npoin', 'nnode', 'ngaus', 'porder', 'ndime']:
        variables.append(struct.unpack('<iii', data_read(folder_id, function_name, variable_name))[1])
        
    return variables 

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

def read_3D_array(x: int, y: int, z: int, val_type: str, array_data):
    # Create empty array of shape -> (x,y) with appropriate type.
    array, byte_size = array_type((x, y, z), val_type)
    
    # Set data index
    data_index = 4
    
    # Unpack and store data on created array.
    for i in range(0, z):
        for j in range(0, y):
            for k in range(0, x):
                array[k][j][i] = struct.unpack('<' + val_type, array_data[data_index : data_index + byte_size])[0]
                data_index += byte_size
    
    return array

def read_4D_array(x: int, y: int, z: int, n: int, val_type: str, array_data):
    # Create empty array of shape -> (x,y) with appropriate type.
    array, byte_size = array_type((x, y, z, n), val_type)
    
    # Set data index
    data_index = 4
    
    # Unpack and store data on created array.
    for i in range(0, n):
        for j in range(0, z):
            for k in range(0, y):
                for l in range(0, x):
                    array[l][k][j][i] = struct.unpack('<' + val_type, array_data[data_index : data_index + byte_size])[0]
                    data_index += byte_size
    
    return array

def load_array(folder_id: int, function_name: str, array_name, shape_variables: list) -> np.ndarray:
    nelem = shape_variables[0]
    npoin = shape_variables[1]
    nnode = shape_variables[2]
    ngaus = shape_variables[3]
    porder = shape_variables[4]
    ndime = shape_variables[5]
    
    array_data = data_read(folder_id, function_name, array_name)
    
    # Load connec array of function full_convec_ijk or full_diffusion_ijk. connec(nelem,nnode)
    if array_name == 'connec':
        array = read_2D_array(x=nelem, y=nnode, val_type='f', array_data=array_data)
                
    # Load Ngp array of function full_convec_ijk or full_diffusion_ijk. Ngp(ngaus,nnode)
    elif array_name == 'Ngp':
        array = read_2D_array(x=ngaus, y=nnode, val_type='f', array_data=array_data)
                
    # Load dNgp array of full_convec_ijk or full_diffusion_ijk. dNgp(ndime,nnode,ngaus)
    elif array_name == 'dNgp':
        array = read_3D_array(x=ndime, y=nnode, z=ngaus, val_type='f', array_data=array_data)

    # Load He array of full_convec_ijk or full_diffusion_ijk. He(ndime,ndime,ngaus,nelem)
    elif array_name == 'He':
        array = read_4D_array(x=ndime, y=ndime, z=ngaus, n=nelem, val_type='f', array_data=array_data)

    # Load xgp array of function full_convec_ijk or full_diffusion_ijk. xgp(ngaus,ndime)
    elif array_name == 'xgp':
        array = read_2D_array(x=ngaus, y=ndime, val_type='f', array_data=array_data)
                
    # Load dlxigp_ip array of full_convec_ijk or full_diffusion_ijk. dlxigp_ip(ngaus,ndime,porder+1)
    elif array_name == 'dlxigp_ip':
        array = read_3D_array(x=ngaus, y=ndime, z=porder+1, val_type='f', array_data=array_data)

    # Load gpvol array of full_convec_ijk or full_diffusion_ijk. gpvol(1,ngaus,nelem)
    elif array_name == 'gpvol':
        array = read_3D_array(x=1, y=ngaus, z=nelem, val_type='f', array_data=array_data)

    # Load atoIJK array of function full_convec_ijk or full_diffusion_ijk. atoIJK(nnode)
    elif array_name == 'atoIJK':
        array = read_array(x=nnode, val_type='i', array_data=array_data)
    
    # Load invAtoIJK array of full_convec_ijk or full_diffusion_ijk. invAtoIJK(porder+1,porder+1,porder+1)
    elif array_name == 'invAtoIJK':
        array = read_3D_array(x=porder+1, y=porder+1, z=porder+1, val_type='i', array_data=array_data)

    # Load gmshAtoI array of function full_convec_ijk or full_diffusion_ijk. gmshAtoI(nnode)
    elif array_name == 'gmshAtoI':
        array = read_array(x=nnode, val_type='i', array_data=array_data)

    # Load gmshAtoJ array of function full_convec_ijk or full_diffusion_ijk. gmshAtoJ(nnode)
    elif array_name == 'gmshAtoJ':
        array = read_array(x=nnode, val_type='i', array_data=array_data)
            
    # Load gmshAtoK array of function full_convec_ijk or full_diffusion_ijk. gmshAtoK(nnode)
    elif array_name == 'gmshAtoK':
        array = read_array(x=nnode, val_type='i', array_data=array_data)
    
    # Load q(npoin,ndime) array of function full_convec_ijk.
    elif array_name == 'q':
        array = read_2D_array(x=npoin, y=ndime, val_type='f', array_data=array_data)
    
    # Load u(npoin,ndime) array of function full_convec_ijk or full_diffusion_ijk.
    elif array_name == 'u':
        array = read_2D_array(x=npoin, y=ndime, val_type='f', array_data=array_data)
    
    # Load rho(npoin) array of function full_convec_ijk or full_diffusion_ijk.
    elif array_name == 'rho':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
            
    # Load pr(npoin) array of function full_convec_ijk.
    elif array_name == 'pr':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
            
    # Load E(npoin) array of function full_convec_ijk.
    elif array_name == 'E':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    
    # Load Rmass(npoin) array of function full_convec_ijk or full_diffusion_ijk.
    elif array_name == 'Rmass':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    
    # Load Rmom(npoin,ndime) array of function full_convec_ijk or full_diffusion_ijk.
    elif array_name == 'Rmom':
        array = read_2D_array(x=npoin, y=ndime, val_type='f', array_data=array_data)
    
    # Load Rener(npoin) array of function full_convec_ijk or full_diffusion_ijk.
    elif array_name == 'Rener':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
            
    # Load Tem(npoin) array of function full_diffusion_ijk.
    elif array_name == 'Tem':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    
    # Load mu_e(nelem, ngaus) array of function full_diffusion_ijk.
    elif array_name == 'mu_e':
        array = read_2D_array(x=nelem, y=ngaus, val_type='f', array_data=array_data)
    
    # Load mu_sgs(nelem, ngaus) array of function full_diffusion_ijk.
    elif array_name == 'mu_sgs':
        array = read_2D_array(x=nelem, y=ngaus, val_type='f', array_data=array_data)
    
    # Load Ml(npoin) array of function full_diffusion_ijk.
    elif array_name == 'Ml':
        array = read_array(x=npoin, val_type='f', array_data=array_data)
    
    # Load mu_fluid(npoin) array of function full_diffusion_ijk.
    elif array_name == 'mu_fluid':
        array = read_array(x=npoin, val_type='f', array_data=array_data)

    return array