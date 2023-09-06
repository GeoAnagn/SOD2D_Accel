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
- conda env create -n {environment_name}
- conda activate {environment_name}
- conda install pip
- pip install numpy
- pip install pandas
- pip install scikit-learn
- pip install opentuner

## Framework Folder Hierarchy
```bash
├── Blackbox_Info.json
├── Blackbox_Tuner.py
├── Blackbox_Original.py
├── Whitebox_Info.json
├── Whitebox_Tuner.py
├── Whitebox_Original.py
├── Example
│   └── ...
├── Functions
│   └── ...
├── Modified_Sod2d_Files
│   ├── Blackbox_Analysis
│   │   └── ...
│   └── Whitebox_Analysis
│       ├── Function_Call_Data_Generators
│       │   └── ...
│       ├── Modified_Functions
│       │   └── ...
│       └── Original_Functions
│           └── ...
└── Parsers
    └── ...
```
- **Blackbox_Original.py**  
  └ Python file for collecting timing information on the original SOD2D execution.
- **Blackbox_Tuner.py**  
  └ Main python file for the blackbox exploration of the modified SOD2D execution.
- **Blackbox_Info.json**  
  └ Blackbox info file for parameter setup. (Detailed explanation in usage section.)
- **Whitebox_Original.py**  
  └ Python file for executing SOD2D with function input and output storing.
- **Whitebox_Tuner.py**  
  └ Main python file for the whitebox exploration of the modified SOD2D execution.
- **Whitebox_Info.json**  
  └ Whitebox info file for parameter setup. (Detailed explanation in usage section.)
- **Example**  
  └ Folder for the example file that the SOD2D will use.
- **Functions**  
  └ Folder where the custom functions, of the main python files, are stored.
- **Modified_Sod2d_Files**
  - **Blackbox_Analysis**  
    └ Folder containing the modified SOD2D files for the blackbox exploration. (Detailed explanation in the first steps in usage section.)
  - **Whitebox_Analysis**
    - **Function_Call_Data_Generators**  
      └ Folder containing the modified SOD2D files for the function call input parameters and result outputs storing.
    - **Modified_Functions**  
      └ Folder containing the standalone versions of the SOD2D functions fot whitebox analysis.
    - **Original_Functions**  
      └ Folder containing the original standalone versions of the SOD2D functions.
- **Parsers**  
  └ Folder containing files for postprocessing.

## Usage

### Blackbox Analysis

#### First Steps
- Move SOD2D example file to the **Example** folder
- Run Blackbox_Original.py, this will create, if not already existent, an ***Archive/Blackbox_Analysis/Original_Example*** path including two files:
  - **time.json**: Contains the total SOD2D execution time
  - **openacc_timing.csv**: Contains the timing info for each function. (Can be used for checking if analysis has found a better solution.)

#### SOD2D File Replacement
In order for the exploration to happen we need to move our modified SOD2D files to the appropriate SOD2D folders and rebuild. All the files that should be replaced are located in the **Modified_Sod2d_Files/Blackbox_Analysis** folder. Bellow is an explanation of where the files should be placed. 

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
- **mod_analysis.f90**  
  └ In this file the function that is modified and needs analysis is ***visc_dissipationRate***.  
  Should be moved to *{sod2d_folder}/src/lib_sod2d/sources*
  ```
  cp Modified_Sod2d_Files/mod_analysis.f90 {sod2d_folder}/src/lib_sod2d/sources
  ```
- **mod_entropy_viscosity.f90**  
  └ In this file the function that is modified and needs analysis is ***smart_visc_spectral***.  
  Should be moved to *{sod2d_folder}/src/lib_sod2d/sources*
  ```
  cp Modified_Sod2d_Files/mod_entropy_viscosity.f90 {sod2d_folder}/src/lib_sod2d/sources
  ```
- **TGVSolver**  
  └ This file must be modified only if simulation parameters need to be changed. (save, time, max_steps, etc.)  
  Should be moved to *{sod2d_folder}/src/lib_mainBaseClass/sources*
  ```
  cp Modified_Sod2d_Files/TGVSolver.f90 {sod2d_folder}/src/lib_mainBaseClass/sources
  ```
- **time_integ.f90**  
  └ In this file the function that is modified and needs analysis is ***rk_4_main***.  
  Should be moved to *{sod2d_folder}/src/lib_sod2d/sources*
  ```
  cp Modified_Sod2d_Files/time_integ.f90 {sod2d_folder}/src/lib_sod2d/sources
  ```

Now in order for SOD2D to support the changes run:
```
cd {sod2d_folder}/build & make clean & make
```

#### Last step before analysis (JSON Explanation)
```json
{   
    "sod2d_path" : "path_of_sod2d_executable",
    "gpu_ids": "0, ..., x",
    "rank_num": "x",
    "dataframe_columns": ["var1",
                          ...,
                          "varN", 
                          "time"],
    "program_end": 1000,
    "parameters": [
        {
            "name": "var1",
            "min" : 1,
            "max" : 1000,
            "multiplier": 512,
            "type": "integer"
        },
        ...,
        {
            "name": "varN",
            "min" : 1,
            "max" : 32,
            "multiplier": 32,
            "type": "integer"
        }
    ]    
}
```

### Whitebox Analysis
