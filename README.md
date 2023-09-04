# SOD2D Acceleration Framework

A custom acceleration framework for CFD Simulation Framework SOD2D using the Opentuner library for the exploration of different combinations of OpenACC directives and loop optimization techniques based on execution time. 

### Prerequisites
- NVIDIA GPU/s. (AMD version pending...)
- SOD2D with GPU support.  
  Follow instructions in https://gitlab.com/bsc_sod2d/sod2d_gitlab/-/wikis/home for detailed installation instructions.
- Anaconda  
  Follow instructions in https://docs.anaconda.com/free/anaconda/install/linux/ for detailed installation instructions.
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
├── Whitebox_Info.json
├── Whitebox_Tuner.py
├── example
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
└── parsers
    └── results_plotter.py
```

- **Blackbox_Tuner.py**  
  Main python file for the blackbox exploration of the SOD2D execution.
- **Blackbox_Info.json**  
  Blackbox info file for parameter setup. (Detailed explanation in usage section.)
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
