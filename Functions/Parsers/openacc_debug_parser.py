"""
A simple OpenACC Debug parser
"""

import pandas as pd

def parser(results_folder):
    """
    Parse .txt OpenACC Debug file.

    results_folder: Folder containing openacc_debug.txt file.
    """

    # Read timing file
    with open(results_folder + '/openacc_timing.txt') as timing:
        lines = timing.readlines()

    # Declare .csv columns and create a Dataframe.
    column_names = [# Generic Info
                    "Log File Index", "File Path", "Function Name", "Line", "Device ID", "Thread ID", "Operation", 
                    # Data Info
                    "Bytes", "Variable", "Device Address", "Queue", 
                    # Compute Info
                    "No. Gangs", "No. Workers", "Vector Length", "Grid", "Block", "Shared Memory",
                    ]
    debug_df = pd.DataFrame(columns=column_names)

    # Create lists for filepath and function name indexing
    filepath_index = []
    index = 0
    for pathline in lines:
        if "device=" in pathline or "threadid=" in pathline:
            filepath_index.append(index)
        index += 1

    for index in range(len(filepath_index)):
        # Initialize debug info dictionary
        debug_dict = {"Log File Index": filepath_index[index]+1,
                      "File Path": '',
                      "Function Name": '',
                      "Line": '',
                      "Device ID": '',
                      "Thread ID": '',
                      "Operation": '',
                      "Bytes": '',
                      "Variable": '',
                      "Device Address": '',
                      "Queue": '',
                      "No. Gangs": '',
                      "No. Workers": '',
                      "Vector Length": '',
                      "Grid": '',
                      "Block": '',
                      "Shared Memory": ''}
        
        # Get debug line
        debug_line = lines[filepath_index[index]]
        # Replace ending newline for better parsing
        debug_line = debug_line.replace('\n', ' ')

        # Collect Enter enter data debug info
        if "Enter enter data construct" in debug_line:
            debug_dict["Operation"] = "Enter enter data construct"
            debug_dict = get_info(debug_line, debug_dict)

        # Collect leave enter data info    
        elif "Leave enter data" in debug_line:
            debug_dict["Operation"] = "Leave enter data"
            debug_dict = get_leave_info(debug_line, debug_dict, "Enter enter data construct", debug_df)
          
        # Collect Enter exit data info
        elif "Enter exit data construct" in debug_line:
            debug_dict["Operation"] = "Enter exit data construct"
            debug_dict = get_info(debug_line, debug_dict)

        # Collect leave exit data info
        elif "Leave exit data" in debug_line:
            debug_dict["Operation"] = "Leave exit data"
            debug_dict = get_leave_info(debug_line, debug_dict, "Enter exit data construct", debug_df)

        # Collect Enter compute region data info
        elif "Enter compute region" in debug_line:
            debug_dict["Operation"] = "Enter compute region"
            debug_dict = get_info(debug_line, debug_dict)

        # Collect Leave compute region data info
        elif "Leave compute" in debug_line:
            debug_dict["Operation"] = "Leave compute"
            debug_dict = get_leave_info(debug_line, debug_dict, "Enter compute region", debug_df)
        
        # Collect Create CUDA data data info
        elif "create CUDA data" in debug_line:
            debug_dict["Operation"] = "Create CUDA data"
            debug_dict = get_info(debug_line, debug_dict)

        # Collect Alloc CUDA data data info
        elif "alloc  CUDA data" in debug_line:
            debug_dict["Operation"] = "alloc CUDA data"
            debug_dict = get_info(debug_line, debug_dict)

        # Collect Upload CUDA data data info
        elif "upload CUDA data" in debug_line:
            debug_dict["Operation"] = "upload CUDA data"
            debug_dict = get_info(debug_line, debug_dict)
        
        # Collect Download CUDA data data info
        elif "download CUDA data" in debug_line:
            debug_dict["Operation"] = "download CUDA data"
            debug_dict = get_info(debug_line, debug_dict)

        # Collect Delete CUDA data data info
        elif "delete CUDA data" in debug_line:
            debug_dict["Operation"] = "delete CUDA data"
            debug_dict = get_info(debug_line, debug_dict)

        # Collect Implicit wait data info
        elif "Implicit wait" in debug_line:
            debug_dict["Operation"] = "Implicit wait"
            debug_dict = get_info(debug_line, debug_dict)

        # Collect Launch CUDA kernel data info
        elif "launch CUDA kernel" in debug_line:
            debug_dict["Operation"] = "launch CUDA kernel"
            debug_dict = get_info(debug_line, debug_dict)

        temp_df = pd.DataFrame([debug_dict])
        debug_df = pd.concat([debug_df, temp_df], ignore_index=True)    

        # Save info found to a .csv file.
        debug_df.to_csv(results_folder + '/openacc_debug.csv', index=False)


def get_info(debug_line: str, debug_dict: dict) -> dict:
    if "file=" in debug_line:
        debug_dict["File Path"] = debug_line[debug_line.find("file=") + len("file=")  : debug_line.find(" ", debug_line.find("file=") + len("file="))]
    if "function=" in debug_line:
        debug_dict["Function Name"] = debug_line[debug_line.find("function=") + len("function=")  : debug_line.find(" ", debug_line.find("function=") + len("function="))]
    if "line=" in debug_line:
        debug_dict["Line"] = debug_line[debug_line.find("line=") + len("line=")  : debug_line.find(" ", debug_line.find("line=") + len("line="))]
    if "device=" in debug_line:
        debug_dict["Device ID"] = debug_line[debug_line.find("device=") + len("device=")  : debug_line.find(" ", debug_line.find("device=") + len("device="))]
    if "threadid=" in debug_line:
        debug_dict["Thread ID"] = debug_line[debug_line.find("threadid=") + len("threadid=")  : debug_line.find(" ", debug_line.find("threadid=") + len("threadid="))]
    if "bytes=" in debug_line:
        debug_dict["Bytes"] = debug_line[debug_line.find("bytes=") + len("bytes=")  : debug_line.find(" ", debug_line.find("bytes=") + len("bytes="))]
    if "variable=" in debug_line:
        debug_dict["Variable"] = debug_line[debug_line.find("variable=") + len("variable=")  : debug_line.find(" ", debug_line.find("variable=") + len("variable="))]
    if "devaddr=" in debug_line:
        debug_dict["Device Address"] = debug_line[debug_line.find("devaddr=") + len("devaddr=")  : debug_line.find(" ", debug_line.find("devaddr=") + len("devaddr="))]
    if "queue=" in debug_line:
        debug_dict["Queue"] = debug_line[debug_line.find("queue=") + len("queue=")  : debug_line.find(" ", debug_line.find("queue=") + len("queue="))]
    if "num_gangs=" in debug_line:
        debug_dict["No. Gangs"] = debug_line[debug_line.find("num_gangs=") + len("num_gangs=")  : debug_line.find(" ", debug_line.find("num_gangs=") + len("num_gangs="))]
    if "num_workers=" in debug_line:
        debug_dict["No. Workers"] = debug_line[debug_line.find("num_workers=") + len("num_workers=")  : debug_line.find(" ", debug_line.find("num_workers=") + len("num_workers="))]
    if "vector_length=" in debug_line:
        debug_dict["Vector Length"] = debug_line[debug_line.find("vector_length=") + len("vector_length=")  : debug_line.find(" ", debug_line.find("vector_length=") + len("vector_length="))]
    if "grid=" in debug_line:
        debug_dict["Grid"] = debug_line[debug_line.find("grid=") + len("grid=")  : debug_line.find(" ", debug_line.find("grid=") + len("grid="))]
    if "block=" in debug_line:
        debug_dict["Block"] = debug_line[debug_line.find("block=") + len("block=")  : debug_line.find(" ", debug_line.find("block=") + len("block="))]
    if "shared memory=" in debug_line:
        debug_dict["Shared Memory"] = debug_line[debug_line.find("shared memory=") + len("shared memory=")  : debug_line.find(" ", debug_line.find("shared memory=") + len("shared memory="))]
    
    return debug_dict

def get_leave_info(debug_line: str, debug_dict: dict, operation: str, debug_df) -> dict:
    for enter_index in range(len(debug_df)-1, -1, -1):
        row = debug_df.iloc[enter_index]
        if (row["Operation"] == operation  and 
            row["Device ID"] == debug_line[debug_line.find("device=") + len("device=")  : debug_line.find(" ", debug_line.find("device=") + len("device="))] and 
            row["Thread ID"] == debug_line[debug_line.find("threadid=") + len("threadid=")  : debug_line.find(" ", debug_line.find("threadid=") + len("threadid="))]):
            
            debug_dict["File Path"] = row["File Path"]
            debug_dict["Function Name"] = row["Function Name"]
            debug_dict["Line"] = row["Line"]
            debug_dict["Device ID"] = row["Device ID"]
            debug_dict["Thread ID"] = row["Thread ID"]
            break

    return debug_dict