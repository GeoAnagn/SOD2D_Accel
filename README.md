# SOD2D Acceleration Framework

A custom acceleration framework for CFD Simulation Framework SOD2D using the Opentuner library for the exploration of different combinations of OpenACC directives and loop optimization techniques based on execution time. 

# Installation
## Quick Installation (Recommended)
If Docker is installed in system you can pull an image having everything setup (30GB required space):  
- docker pull geoaiagia/refmap:february

## Local Installation
### Prerequisites
- NVIDIA GPU/s. (AMD version pending...)
- NVIDIA HPC SDK  
  └ Follow instructions in [here](https://developer.nvidia.com/nvidia-hpc-sdk-233-downloads) for detailed installation.
- Anaconda  
  └ Follow instructions in [here](https://docs.anaconda.com/free/anaconda/install/linux/) for detailed installation.
- Intel Vtune  
  └ Follow instructions in [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/vtune-profiler-download.html) for detailed installation


### Environment Instalation instructions
Supposing you have not used our Dockerfile and you have all the presequisites setup, follow the following  

Instructions to setup system package presequities:
- apt -y update && apt -y upgrade
- apt install -y build-essential git cmake gfortran ninja-build wget pkg-config vim
- apt install -y binutils binutils-dev libiberty-dev
- apt install -y software-properties-common
- apt install -y libssl-dev libc6-dev
- apt install -y libdrm-dev
- apt install -y python3
- apt install -y gmsh
- apt install -y btop
- apt install -y pip

Instructions to create a python environment and install the needed libraries:
- conda create --name {environment_name} python=3.11
- conda activate {environment_name}
- conda install pip (if not already installed)
- pip install numpy
- pip install pandas
- pip install openpyxl
- pip install seaborn
- pip install matplotlib
- pip install scikit-learn
- pip install opentuner
- pip install alive_progress

Instructions to setup rest of dependencies can be found on [SOD2D Gitlab Repository](https://gitlab.com/bsc_sod2d/sod2d_gitlab/-/tree/02b47d1f91ffbf4c4ac342d0773015b9243c4e2c) make sure it is version 02b47d1f91ffbf4c4ac342d0773015b9243c4e2c

# Usage

Firstly download this repository:
```
git clone https://github.com/GeoAnagn/SOD2D_Accel.git
cd {path}/SOD2D_Accel
```
Next move to the SOD2D_Versions folder and download SOD2D:
```
cd SOD2D_Versions
git clone https://gitlab.com/bsc_sod2d/sod2d_gitlab.git
cd sod2d_gitlab
git checkout 02b47d1f91ffbf4c4ac342d0773015b9243c4e2c
cd ../../
```

**IMPORTANT: Set your Taylor Green Vortex parameters before procceding.**  
Parameters in file **sod2d_gitlab/src/lib_mainBaseClass/sources/TGVSolver.f90**   
Enable tgv in file **sod2d_gitlab/src/app_sod2/sources/sod2d.f90**

Now run the container as follows:
```
docker run --gpus all -dit --privileged -v {path_to_repo}:/home/apps/refmap geoaiagia/refmap:february /bin/bash
```
Lastlymove to the **SOD2D_Accel/Scripts** folder and execute the **Initialize_SOD2D_versions.sh** to create the modified SOD2D versions.
```
cd Scripts
./Initialize_SOD2D_versions.sh
cd ..
```

Now everything is setup to start profiling and auto-tuning.

## Framework Folder Hierarchy
```bash
├── Original_SOD2D_CPU.py
├── Original_SOD2D_GPU.py
├── Blackbox_Tuner.py
├── JSONs
│   ├── Original_SOD2D_CPU_Info.json
│   ├── Original_SOD2D_CPU_Info.json
│   └── Blackbox_Versions
│       ├── Blackbox_Coupled.json
│       └── Blackbox_Decoupled.json
├── Scripts
│   └── Initialize_SOD2D_versions.sh
├── Example
│   └── ...
├── SOD2D_Versions
│   └── ...
├── Functions
│   ├── Parsers
│   └── ...
└── Modified_Sod2d_Files
    ├── Blackbox_Analysis
    │   └── ...
    └── Nsight_Modified_Files
        └── ...
 
```
- **Original_SOD2D_CPU.py**  
  └ Python file for the CPU Version SOD2D execution.
- **Original_SOD2D_GPU.py**  
  └ Python file for the GPU Version SOD2D execution.
- **Blackbox_Tuner.py**  
  └ Main python file for the blackbox exploration of the modified SOD2D execution.
- **Scripts**  
  └ Folder containing the scripts for the SOD2D file replacement.
- **Example**  
  └ Folder for the example file that the SOD2D will use.
- **SOD2D_Versions**  
  └ Folder for the SOD2D cloned repository.
- **JSONs**
  - **Original_SOD2D_CPU_Info.py**  
    └ CPU Version SOD2D info file for parameter setup. (Detailed explanation in usage section.)
  - **Original_SOD2D_GPU_Info.py**  
    └ GPU Version SOD2D info file for parameter setup. (Detailed explanation in usage section.)
  - **Blackbox_Info.json**  
    └ Blackbox info file for parameter setup. (Detailed explanation in usage section.)
- **Functions**  
  └ Folder where the custom functions, of the main python files, are stored.
- **Modified_Sod2d_Files**
  - **Blackbox_Analysis**  
    └ Folder containing the modified SOD2D files for the blackbox exploration. 
  - **Nsight_Modified_Files**  
    └ Folder containing the modified SOD2D files for the SOD2D GPU profiling. 



## Profiling

### CPU Profiling

In order to profile the execution of SOD2D on a CPU the following steps must be followed:

First modify the **Original_SOD2D_CPU_info.json** located in the **JSONs** folder
```json
{
    "sod2d_path": "./SOD2D_Versions/sod2d_cpu/build/src/app_sod2d/sod2d",
    "example_path": "{Example Path}",
    "results_path": "{result_path_name}",
    "rank_num": "{number_of_threads_to_use}",
    "metric": "{intel_vtune_profiling_option}"
}
```
- **sod2d_path**  
  └ Path of sod2d executable.
- **example_path**  
  └ Folder of sod2d example.
- **results_path**  
  └ Folder of sod2d result files.
- **rank_num**  
  └ Number of CPU threads to be used.
- **metric**  
  └ Intel Vtune profiling metric (i.e. hotspots, performance_snapshot, threading)

DO NOT FORGET to enable V-Tune variables:
```
source /opt/intel/oneapi/vtune/2024.1/vtune-vars.sh
```

Next execute the following command to start the profiling:
```
python3 Original_SOD2D_CPU.py
```
The result files will be automatically moved to the **results_path** provided

### GPU Profiling

In order to profile the execution of SOD2D on GPU/s the following steps must be followed:

First modify the **Original_SOD2D_GPU_info.json** located in the **JSONs** folder
```json
{
    "sod2d_path": "./SOD2D_Versions/sod2d_profiling/build/src/app_sod2d/sod2d",
    "example_path": "{Example Folder}",
    "results_path": "{Results Folder}",
    "gpu_ids": "{GPU ID/s}",
    "rank_num": "{Number of GPUs to use}"
}
```
- **sod2d_path**  
  └ Path of sod2d executable.
- **example_path**  
  └ Folder of sod2d example.
- **results_path**  
  └ Folder of sod2d result files.
- **gpu_ids**  
  └ Which GPU/s to use. In the case of multiple GPUs seperate the ids by comma (i.e. 1,2,3,4)  
  **NOTE:** If unsure what GPU ID to use run on the terminal:
  ```
  export CUDA_DEVICE_ORDER=PCI_BUS_ID
  nvidia-smi
  ```
- **rank_num**  
  └ Number of GPUs to be used.

Next execute the following command to start the profiling:
```
python3 Original_SOD2D_GPU.py
```
The result files will be automatically moved to the **results_path** provided  

## Blackbox Analysis

### Coupled Version

In order to autotune the execution of SOD2D on GPU/s using the Blackbox Coupled variant the following steps must be followed:

First modify the **Blackbox_Coupled.json** located in the **JSONs/Blackbox_Versions** folder
```json
{   
    "sod2d_path": "./SOD2D_Versions/sod2d_blackbox_coupled/build/src/app_sod2d/sod2d",
    "og_res_path": "{Original SOD2D Results Folder}",
    "example_path": "{Example Folder}",
    "results_path": "{Results Path}",
    "gpu_ids": "{GPU ID/s}",
    "rank_num": "{Number of GPUs to use}",
    "dataframe_columns":[
                            "gang_num", 
                            "worker_num", 
                            "vector_num",
                            "time"
                        ],
    "repetitions": No.Repetitions,
    "configs_to_check": No.Configurations,
    "parameters": [
        {
            "name": "gang_num",
            "min" : Min_Value, // 1 (Recommended)
            "max" : Max_Value, // 250 (Recommended)
            "multiplier": Value_Multiplier, // 2048 (Recommended)
            "type": "integer"
        },
        {
            "name": "worker_num",
            "min" : 1,
            "max" : 1,
            "multiplier": 1,
            "type": "integer"
        },
        {
            "name": "vector_num",
            "min" : Min_Value, // 1 (Recommended)
            "max" : Max_Value, // 20 (Recommended)
            "multiplier": Value_Multiplier, // 32 (Recommended)
            "type": "integer"
        }
    ]    
}
```
- **sod2d_path**  
  └ Path of sod2d executable.
- **og_res_path**  
  └ Original SOD2D execution results. (Original_SOD2D_GPU.py can be used to generate results)
- **example_path**  
  └ Folder of sod2d example.
- **results_path**  
  └ Folder of sod2d result files.
- **gpu_ids**  
  └ Which GPU/s to use. In the case of multiple GPUs seperate the ids by comma (i.e. 1,2,3,4)  
  **NOTE:** If unsure what GPU ID to use run on the terminal:
  ```
  export CUDA_DEVICE_ORDER=PCI_BUS_ID
  nvidia-smi
  ```
- **rank_num**  
  └ Number of GPUs to be used.
- **dataframe_columns**  
  └ Columns for xlsx containing all configurations tested and timing results.
- **repetitions**  
  └ Integer value used to stop the autotuner in case it repeats the same configurations over and over/
- **configs_to_check**  
  └ Integer value used to define number of successful configurations tested.
- **parameters**  
  └ Parameters that the autotuner will explore. Each parameter consists of:  
  - **name:** Name of variable must exist in dataframe columns.
  - **min:**  Minimum possible value of parameter.
  - **max:**  Maximum possible value of parameter.
  - **multiplier:** Multiplier for the parameter value.
  - **type:** Parameter type. (integer is currently supported)



Next execute the following command to start the profiling:
```
python3 Blackbox_Tuner.py coupled
```
The result files will be automatically moved to the **results_path** provided  

### Decoupled Version
In order to autotune the execution of SOD2D on GPU/s using the Blackbox Coupled variant the following steps must be followed:

First modify the **Blackbox_Decoupled.json** located in the **JSONs/Blackbox_Versions** folder
```json
{   
    "sod2d_path": "./SOD2D_Versions/sod2d_blackbox_decoupled/build/src/app_sod2d/sod2d",
    "og_res_path": "{Original SOD2D Results Folder}",
    "example_path": "{Example Folder}",
    "results_path": "{Results Path}",
    "gpu_ids": "{GPU ID/s}",
    "rank_num": "{Number of GPUs to use}",
    "dataframe_columns": [
                            "gang_num_fcijk", 
                            "worker_num_fcijk", 
                            "vector_num_fcijk",
                            "gang_num_fdijk", 
                            "worker_num_fdijk", 
                            "vector_num_fdijk",
                            "gang_num_vaek",
                            "worker_num_vaek",
                            "vector_num_vaek",
                            "gang_num_svs",
                            "worker_num_svs",
                            "vector_num_svs",
                            "gang_num_vdr",
                            "worker_num_vdr",
                            "vector_num_vdr",
                            "time"
                        ],
    "repetitions": No.Repetitions,
    "configs_to_check": No.Configurations,
    "parameters": [
        {
            "name": "gang_num",
            "min" : Min_Value, // 1 (Recommended)
            "max" : Max_Value, // 250 (Recommended)
            "multiplier": Value_Multiplier, // 2048 (Recommended)
            "type": "integer"
        },
        {
            "name": "worker_num",
            "min" : 1,
            "max" : 1,
            "multiplier": 1,
            "type": "integer"
        },
        {
            "name": "vector_num",
            "min" : Min_Value, // 1 (Recommended)
            "max" : Max_Value, // 20 (Recommended)
            "multiplier": Value_Multiplier, // 32 (Recommended)
            "type": "integer"
        }

        // Repeat for all the variables in the dataframe columns EXEPT time
    ]    
}
```

- **sod2d_path**  
  └ Path of sod2d executable.
- **og_res_path**  
  └ Original SOD2D execution results. (Original_SOD2D_GPU.py can be used to generate results)
- **example_path**  
  └ Folder of sod2d example.
- **results_path**  
  └ Folder of sod2d result files.
- **gpu_ids**  
  └ Which GPU/s to use. In the case of multiple GPUs seperate the ids by comma (i.e. 1,2,3,4)  
  **NOTE:** If unsure what GPU ID to use run on the terminal:
  ```
  export CUDA_DEVICE_ORDER=PCI_BUS_ID
  nvidia-smi
  ```
- **rank_num**  
  └ Number of GPUs to be used.
- **dataframe_columns**  
  └ Columns for xlsx containing all configurations tested and timing results.
- **repetitions**  
  └ Integer value used to stop the autotuner in case it repeats the same configurations over and over/
- **configs_to_check**  
  └ Integer value used to define number of successful configurations tested.
- **parameters**  
  └ Parameters that the autotuner will explore. Each parameter consists of:  
  - **name:** Name of variable must exist in dataframe columns.
  - **min:**  Minimum possible value of parameter.
  - **max:**  Maximum possible value of parameter.
  - **multiplier:** Multiplier for the parameter value.
  - **type:** Parameter type. (integer is currently supported)



Next execute the following command to start the profiling:
```
python3 Blackbox_Tuner.py decoupled
```
The result files will be automatically moved to the **results_path** provided  
