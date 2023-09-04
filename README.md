# SOD2D Acceleration Framework

A custom acceleration framework for CFD Simulation Framework SOD2D using the Opentuner library for the exploration of different combinations of OpenACC directives and loop optimization techniques based on execution time. 

## Instalation

### Prerequisites
- NVIDIA GPU/s. (AMD version pending...)
- SOD2D with GPU support.  
  Follow instructions in https://gitlab.com/bsc_sod2d/sod2d_gitlab/-/wikis/home for detailed installation instructions.
- Anaconda  
  Follow instructions in https://docs.anaconda.com/free/anaconda/install/linux/ for detailed installation.
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
    └── results_plotter.py
```
- **Blackbox_Original.py**  
  Python file for collecting timing information on the original SOD2D execution.
- **Blackbox_Tuner.py**  
  Main python file for the blackbox exploration of the SOD2D execution.
- **Blackbox_Info.json**  
  Blackbox info file for parameter setup. (Detailed explanation in usage section.)
- **Whitebox_Original.py**  
  Python file for executing SOD2D with function input and output storing.
- **Whitebox_Tuner.py**  
  Main python file for the whitebox exploration of the SOD2D execution.
- **Whitebox_Info.json**  
  Whitebox info file for parameter setup. (Detailed explanation in usage section.)
- **Example**  
  Folder for the example file that the SOD2D will use.
- **Functions**  
  Folder where the custom functions, of the main python files, are stored.
- **Modified_Sod2d_Files**
  - **Blackbox_Analysis**  
    Folder containing the modified SOD2D files for the blackbox exploration. (Detailed explanation in the first steps in usage section.)
  - **Whitebox_Analysis**
    - **Function_Call_Data_Generators**  
      Folder containing the modified SOD2D files for the function call input parameters and result outputs storing.
    - **Modified_Functions**  
      Folder containing the standalone versions of the SOD2D functions fot whitebox analysis.
    - **Original_Functions**  
      Folder containing the original standalone versions of the SOD2D functions.
- **Parsers**  
  Folder containing files for postprocessing.

## Usage

### Blackbox Analysis

#### First Steps
- Move SOD2D example file to the **Example** folder
- Run Blackbox_Original.py, this will create, if not already existent, an ./Archive/Blackbox_Analysis/Original_Example path including two files:
  - **time.json**: Contains the total SOD2D execution time
  - **openacc_timing.csv**: Contains the timing info for each function. (Can be used for checking if analysis has found a better solution.)

#### Last step before analysis (JSON Explanation)
```json
{   
    "sod2d_path" : "path_of_sod2d_executable",
    "gpu_ids": "0, ..., x",
    "rank_num": "x",
    "dataframe_columns": ["var1"
                          ...
                          "varN" 
                          "time"],
    "program_end": Y,
    "parameters": [
        {
            "name": "var1",
            "min" : 1,
            "max" : 1000,
            "multiplier": 512,
            "type": "integer"
        },
        ...
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
