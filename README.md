# SOD2D Acceleration Framework

A custom acceleration framework for CFD Simulation Framework SOD2D using the Opentuner library for the exploration of different combinations of OpenACC directives and loop optimization techniques based on execution time. 

## Instalation

### Prerequisites
- NVIDIA GPU/s. (AMD version pending...)
- SOD2D with GPU support.  
  └ Follow instructions in https://gitlab.com/bsc_sod2d/sod2d_gitlab/-/wikis/home for detailed installation instructions.
- Anaconda  
  └ Follow instructions in https://docs.anaconda.com/free/anaconda/install/linux/ for detailed installation.
- A Dockerfile will be uploaded for plug and play capability.

### Environment Instalation instructions
Supposing you have not used our Dockerfile and you have all the presequisites setup, follow the following instructions to create a python environment and install the needed libraries:
- conda create --name {environment_name} python=3.11
- conda activate {environment_name}
- conda install pip (if not already installed)
- pip install numpy
- pip install pandas
- pip install scikit-learn
- pip install opentuner
- pip install alive_progress

## Framework Folder Hierarchy
```bash
├── Blackbox_Original.py
├── Blackbox_Tuner.py
├── Whitebox_Data_Postprocessing.py
├── Whitebox_Original.py
├── Whitebox_Tuner.py
├── Configuration_Tester.py
├── JSONs
│   ├── Blackbox_Info.json
│   └── Whitebox_Info.json
├── Scripts
├── Example
│   └── ...
├── Functions
│   ├── Parsers
│   └── ...
└── Modified_Sod2d_Files
    ├── Blackbox_Analysis
    │   └── ...
    └── Whitebox_Analysis
        ├── Function_Call_Data_Generators
        │   └── ...
        ├── Function_Call_Data_Postprocessing
        │   └── ...
        ├── Modified_Functions
        │   └── ...
        └── Original_Functions
            └── ...
 
```
- **Blackbox_Original.py**  
  └ Python file for collecting timing information on the original SOD2D execution.
- **Blackbox_Tuner.py**  
  └ Main python file for the blackbox exploration of the modified SOD2D execution.
- **Whitebox_Data_Postprocessing.py**  
  └ Python file for postprocessing the function input and output storing.
- **Whitebox_Original.py**  
  └ Python file for executing SOD2D with function input and output storing.
- **Whitebox_Tuner.py**  
  └ Main python file for the whitebox exploration of the modified SOD2D execution.
- **Configuration_Tester.py**  
  └ Python file for testing the configurations of the SOD2D exploration results.
- **Scripts**  
  └ Folder containing the scripts for the SOD2D file replacement.
- **Example**  
  └ Folder for the example file that the SOD2D will use.
- **JSONs**
  - **Original_Info.json**  
    └ Info file for SOD2D execution. (Detailed explanation in usage section.)
  - **Blackbox_Info.json**  
    └ Blackbox info file for parameter setup. (Detailed explanation in usage section.)
  - **Whitebox_Data.json**  
    └ Whitebox data storing file for parameter setup. (Detailed explanation in usage section.)
  - **Whitebox_Info.json**  
    └ Whitebox info file for parameter setup. (Detailed explanation in usage section.)
- **Functions**  
  └ Folder where the custom functions, of the main python files, are stored.
- **Modified_Sod2d_Files**
  - **Blackbox_Analysis**  
    └ Folder containing the modified SOD2D files for the blackbox exploration. (Detailed explanation in the first steps in usage section.)
  - **Whitebox_Analysis**
    - **Function_Call_Data_Generators**  
      └ Folder containing the modified SOD2D files for the function call input parameters and result outputs storing.
    - **Function_Call_Data_Postprocessing**  
      └ Folder containing the python files for the function call input parameters and result outputs postprocessing.
    - **Modified_Functions**  
      └ Folder containing the standalone versions of the SOD2D functions fot whitebox analysis.
    - **Original_Functions**  
      └ Folder containing the original standalone versions of the SOD2D functions.

## Usage

### Blackbox Analysis

#### First Steps
- Move SOD2D example file to the **Example** folder
- Modify ***JSONs/Original_Info.json*** with the appropriate SOD2D path, example path, GPU ids, and rank number.
- Run ***Blackbox_Original.py***, this will create, if not already existent, an ***Archive/Blackbox_Analysis/Original_Example*** path including two files:
  - **time.json**: Contains the total SOD2D execution time
  - **openacc_timing.csv**: Contains the timing info for each function. (Can be used for checking if analysis has found a better solution.)

#### SOD2D File Replacement

In order for the exploration to happen we need to move our modified SOD2D files to the appropriate SOD2D folders and rebuild. All the files that should be replaced are located in the **Modified_Sod2d_Files/Blackbox_Analysis** folder. Bellow is an explanation of where the files should be placed. 

#### Automated File Replacement

Run the ***Scripts/Blackbox_File_Replacement.sh*** script. This will automatically move the files to the appropriate folders.

#### Manual File Replacement

- **elem_convec.f90**  
  └ In this file the function that is modified and needs analysis is ***full_convec_ijk***.  
  Should be moved to *{sod2d_folder}/src/lib_sod2d/sources*
  ```
  cp Modified_Sod2d_Files/Blackbox_Analysis/elem_convec.f90 {sod2d_folder}/src/lib_sod2d/sources
  ```
- **elem_diffu.f90**  
  └ In this file the function that is modified and needs analysis is ***full_diffu_ijk***.  
  Should be moved to *{sod2d_folder}/src/lib_sod2d/sources*
  ```
  cp Modified_Sod2d_Files/Blackbox_Analysis/elem_diffu.f90 {sod2d_folder}/src/lib_sod2d/sources
  ```
- **mod_analysis.f90**  
  └ In this file the function that is modified and needs analysis is ***visc_dissipationRate***.  
  Should be moved to *{sod2d_folder}/src/lib_sod2d/sources*
  ```
  cp Modified_Sod2d_Files/Blackbox_Analysis/mod_analysis.f90 {sod2d_folder}/src/lib_sod2d/sources
  ```
- **mod_entropy_viscosity.f90**  
  └ In this file the function that is modified and needs analysis is ***smart_visc_spectral***.  
  Should be moved to *{sod2d_folder}/src/lib_sod2d/sources*
  ```
  cp Modified_Sod2d_Files/Blackbox_Analysis/mod_entropy_viscosity.f90 {sod2d_folder}/src/lib_sod2d/sources
  ```
- **TGVSolver**  
  └ This file must be modified only if simulation parameters need to be changed. (save, time, max_steps, etc.)  
  Should be moved to *{sod2d_folder}/src/lib_mainBaseClass/sources*
  ```
  cp Modified_Sod2d_Files/Blackbox_Analysis/TGVSolver.f90 {sod2d_folder}/src/lib_mainBaseClass/sources
  ```
- **time_integ.f90**  
  └ In this file the function that is modified and needs analysis is ***rk_4_main***.  
  Should be moved to *{sod2d_folder}/src/lib_sod2d/sources*
  ```
  cp Modified_Sod2d_Files/Blackbox_Analysis/time_integ.f90 {sod2d_folder}/src/lib_sod2d/sources
  ```

Now in order for SOD2D to support the changes run:
```
cd {sod2d_folder}/build & make clean & make
```

#### Last step before analysis (JSON Explanation)
```json
{   
    "sod2d_path" : "path_of_sod2d_executable", // Default is /home/apps/sod2d/build/src/app_sod2d/sod2d
    "gpu_ids": "0, ..., No.GPUs",
    "rank_num": "No.GPUs",
    "dataframe_columns": ["gang_num_fcijk", 
                          "worker_num_fcijk", 
                          "vector_num_fcijk", 
                          "gang_num_fdijk", 
                          "worker_num_fdijk", 
                          "vector_num_fdijk",
                          "gang_num_vdr",
                          "worker_num_vdr",
                          "vector_num_vdr",
                          "gang_num_svs",
                          "worker_num_svs",
                          "vector_num_svs",
                          "time"],
    "repetitions": No.Repetitions,
    "configs_to_check": No.Configurations,
    "parameters": [
        {
            "name": "gang_num_fcijk",
            "min" : Min_Value, // Default 1
            "max" : Max_Value, // Default 1000
            "multiplier": Value_Multiplier, // Default 512
            "type": "integer"
        },
        {
            "name": "worker_num_fcijk",
            "min" : 1,
            "max" : 1,
            "multiplier": 1,
            "type": "integer"
        },
        {
            "name": "vector_num_fcijk",
            "min" : Min_Value, // Default 1
            "max" : Max_Value, // Default 32
            "multiplier": Value_Multiplier, // Default 32
            "type": "integer"
        },

        // Repeat for all the variables in the dataframe columns EXEPT time
    ]    
}
```

- **sod2d_path**    
  └ Path of the SOD2D executable.
- **gpu_ids**  
  └ GPU ids to be used for the execution. (Separated by commas)
- **rank_num**  
  └ Rank number to be used for the execution.
- **dataframe_columns**  
  └ Columns of the dataframe that will be created. (Should include all the parameters and the time column.)
- **repetitions**  
  └ Number of consecutive executions that have configurations already checked. (If reached end exploration.)
- **configs_to_check**  
  └ Number of configurations to be checked in each exploration. (Should be less than the total number of configurations.)
- **parameters**
  - **name**  
    └ Name of the parameter.
  - **min**  
    └ Minimum value of the parameter.
  - **max**  
    └ Maximum value of the parameter.
  - **multiplier**  
    └ Multiplier of the parameter. (Used for integer parameters.)
  - **type**  
    └ Type of the parameter. (integer or float)

#### Analysis
- Modify ***JSONs/Blackbox_Info.json*** with the appropriate parameters.
- Run ***Blackbox_Tuner.py***, this will create, if not already existent, an ***Results/Tuner_Results/*** path including one file and one folder:
  - **results.csv**: Contains the results of the exploration.
  - **Configs**: Folder containing the results of each configuration seperately.

- All the files in the ***Results/Tuner_Results/*** will be moved to ***Archive/Blackbox_Analysis/Modified_Folder*** upon completion of the exploration.

Using the results of the exploration we can now modify the SOD2D files with the best configuration and rebuild SOD2D.

### Whitebox Analysis

#### SOD2D File Replacement
In order for the exploration to happen we need to move our modified SOD2D files to the appropriate SOD2D folders and rebuild. All the files that should be replaced are located in the **Modified_Sod2d_Files/Whitebox_Analysis/Function_Call_Data_Generators** folder. Bellow is an explanation of where the files should be placed.

#### Automated File Replacement

Run the ***Scripts/Whitebox_Data_Generators.sh*** script. This will automatically move the files to the appropriate folders.

#### Manual File Replacement

- **elem_convec.f90**  
  └ In this file the function that is modified and needs analysis is ***full_convec_ijk***.  
  Should be moved to *{sod2d_folder}/src/lib_sod2d/sources*
  ```
  cp Modified_Sod2d_Files/elem_convec.f90 {sod2d_folder}/src/lib_sod2d/sources
  ```
- **elem_diffu.f90**
  └ In this file the function that is modified and needs analysis is ***full_diffu_ijk***.  
  Should be moved to *{sod2d_folder}/src/lib_sod2d/sources*
  ```
  cp Modified_Sod2d_Files/elem_diffu.f90 {sod2d_folder}/src/lib_sod2d/sources
  ```

These version of the files include code that will store the input parameters and the result of the function call in a file until the program ends or user interrupt occurs. The files are stored in the ***Archive/Whitebox_Analysis/Data/Data_{Func_Call_IDX}/{function_name}*** folder.  

Now in order for SOD2D to support the changes run:
```
cd {sod2d_folder}/build & make clean & make
```

**(Note)** More functions will be added in future versions.

#### Data Collection

In order to collect the input and output data of the functions that are needed to be analysed, the ***Whitebox_Data_Postproccesing.py*** is used. As a first step it executes the SOD2D with the modified files and collects the data. Some data collection parameters are defined in the ***JSONs/Whitebox_Data.json*** file and are explained bellow.

```json
{
    "sod2d_path": "path_of_sod2d_executable",
    "example_path": "Example",
    "gpu_ids": "0, ..., No.GPUs",
    "rank_num": "No.GPUs",

    "store_path": "Archive/Whitebox_Analysis/Data",
    "start_num" : 0,
    "max_num" : No.MaxCalls,
    "step" : No.Step,
    "cluster_num" : No.Cluster
}
```

- **sod2d_path**    
  └ Path of the SOD2D executable.
- **example_path**  
  └ Path of the example executable. Default is ***Example***.
- **gpu_ids**  
  └ GPU ids to be used for the execution. (Separated by commas)
- **rank_num**  
  └ Rank number to be used for the execution.
- **store_path**  
  └ Path of the data storage. Default is ***Archive/Whitebox_Analysis/Data***.
- **start_num**  
  └ Number of the first function call to be stored.
- **max_num**  
  └ Maximum number of function calls to be stored.
- **step**  
  └ Number of function calls to be skipped between each stored function call.
- **cluster_num**  
  └ Number of clusters to be created for each function.

#### Data Analysis

Due to high number of function calls that are stored, the data analysis is seperated in two parts. The first part is the data postprocessing and the second part is the data analysis. The data postprocessing is done in order to create a dataframe that will be used for the analysis. The data analysis is done in order to find the most representative distinct function calls.

##### Data Postprocessing

The data postprocessing is also done by the ***Whitebox_Data_Postprocessing.py*** python file. The file as a first step creates **{function_name}.csv** files that contain the data of each functions' function calls. These files are stored at ***Archive/Whitebox_Analysis/Data*** folder. 

##### Data Clustering

Next these files are loaded and by applying a clustering algortihm the function calls are seperated in clusters based on the information they take as input and output. The number of clusters created is defined in the ***Whitebox_Data.json***. 

Lastly the function call closest to the centroid of each cluster is selected and the a list of the most representative function calls is created for each function and saved to *.json* file for later use. The *.json* file is stored at ***Archive/Whitebox_Analysis/Data*** folder.  

#### Per-function Executables

In order to execute the Whitebox flow versions of the functions, to be analysed, must be made seperately of the whole SOD2D program and their respective fortran files. 

Precompiled executables containing only the needed function have been made for the **full_convec_ijk** and **full_diffu_ijk**. 

These files are located in the:
- ***Modified_Sod2d_Files/Whitebox_Analysis/Original_Functions*** and 
- ***Modified_Sod2d_Files/Whitebox_Analysis/Modified_Functions*** folders. 

More specifically the folder Original_Functions contains the original versions of the functions and the Modified_Functions contains the modified versions of the functions. The files in the folders are:

- full_convec_ijk.f90  
  The standalone version of the full_convec_ijk function of the file elem_convec.f90.

- full_diffu_ijk.f90  
  The standalone version of the full_diffu_ijk function of the file elem_diffu.f90.

- mod_constants.f90 & mod_nvtx.f90  
  The files needed for the compilation of the standalone versions.

- build.sh  
  The script that compiles the standalone versions and creates the executables.

All the executables generated are found in these folders. 

#### Execution

In order to execute the Whitebox flow versions of the functions, to be analysed, the ***Whitebox_Original.py*** and ***Whitebox_Tuner.py*** python files are used. The ***Whitebox_Original.py*** file executes the original versions of the functions and the ***Whitebox_Tuner.py*** file executes the modified versions of the functions.

##### JSON Explanation
```json
{   
    "data_path" : "Archive/Whitebox_Analysis/Data",
    "func_path" : "Modified_Sod2d_Files/Whitebox_Analysis",
    "func_ver" :  "/Modified_Functions",
    "func_exec" : "/run_elem_diffu",
    "func_name" : "elem_diffu",
    "gpu_ids": "0, ..., x",
    "rank_num": "x",
    "func_call_idx": [],


    "dataframe_columns": [
        "var1",
        "varN",
        "time"
    ],
    "program_end": 1000,
    "parameters": [
        {
            "name": "var1",
            "min": 1,
            "max": 1000,
            "multiplier": 512,
            "type": "integer"
        },
        {
            "name": "varN",
            "min": 1,
            "max": 32,
            "multiplier": 32,
            "type": "integer"
        }
    ]
}
```

- **data_path**  
  └ Path of the data storage. Default is ***Archive/Whitebox_Analysis/Data***.
- **func_path**  
  └ Path of the folder containing the function executables. Default is ***Modified_Sod2d_Files/Whitebox_Analysis***.
- **func_ver**  
  └ Version of the function to be executed. (***Original*** or ***Modified***)
- **func_exec**    
  └ Name of the function executable.
- **func_name**  
  └ Name of the function.
- **gpu_ids**  
  └ GPU ids to be used for the execution. (Separated by commas)
- **rank_num**  
  └ Rank number to be used for the execution.
- **dataframe_columns**  
  └ Columns of the dataframe that will be created. (Should include all the parameters and the time column.)
- **repetitions**  
  └ Number of consecutive executions that have configurations already checked. (If reached end exploration.)
- **configs_to_check**  
  └ Number of configurations to be checked in each exploration. (Should be less than the total number of configurations.)
- **parameters**  
  - **name**  
    └ Name of the parameter.
  - **min**  
    └ Minimum value of the parameter.
  - **max**  
    └ Maximum value of the parameter.
  - **multiplier**  
    └ Multiplier of the parameter. (Used for integer parameters.)
  - **type**  
    └ Type of the parameter. (integer or float)
- **func_call_idx**  
  └ List of the function call indices to be executed. This must be filled with the appropriate results found in ***Archive/Whitebox_Analysis/Data/representative_calls.json***.

##### Original Execution

More specifically the ***Whitebox_Original.py*** file executes the original versions of the function needed using the representative stored input parameters and the results of the function calls and collects timing data. The timing data are stored in the ***Archive/Whitebox_Analysis/Original_Results/{Function Name}*** folder.

##### Tuner Execution

More specifically the ***Whitebox_Tuner.py*** file executes the modified versions of the function needed using the representative stored input parameters and the results of the function calls and collects timing data. The timing data are stored in the ***Archive/Whitebox_Analysis/Tuner_Results/{Function Name}*** folder together with the configuration info.




