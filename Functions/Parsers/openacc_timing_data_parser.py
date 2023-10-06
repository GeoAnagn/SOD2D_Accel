"""
A simple OpenACC Timing parser
"""

import pandas as pd

def parser(results_folder):
    """
    Parse .txt OpenACC Timing file.
    
    results_folder: Folder containing openacc_timing.txt file.
    """

    # Read timing file
    with open(results_folder + '/openacc_timing.txt') as timing:
        lines = timing.readlines()

    # Declare .csv columns and create a Dataframe.
    column_names = ["File Path", "Function Name", "Total Time", "Compute Time", "Data Time"]
    timing_df = pd.DataFrame(columns=column_names)

    # Create lists for filepath and function name indexing
    filepath_index = []
    function_names = []

    # Get file indexes for filepath and function name
    index = 0
    for pathline in lines:
        if "/" in pathline:
            filepath_index.append(index)
            function_names.append(lines[index+1])
        index += 1

    # For loop fix
    filepath_index.append(len(lines))

    # For every filepath found run:
    for i in range(len(filepath_index) - 1):
        # Declare filepath start and end.
        filepath_index_start = filepath_index[i]
        filepath_index_end = filepath_index[i+1] - 1
        filepath = lines[filepath_index_start].replace('\n', '')

        # Initialize total, compute, data times.
        compute_time = 0
        data_time = 0
        total_time = 0

        # Start searching for timing results.
        for j in range(filepath_index_start, filepath_index_end):
            # Get function name.
            if lines[j] in function_names:
                function_name = lines[j].replace(' ', '')
                function_name = function_name[:function_name.find("NVIDIA")]
            # Get all times for function found
            if "total=" in lines[j]:
                total_str = lines[j].replace(' ', '').replace('\n', '').replace(',', '')
                total_str = total_str[total_str.find("total=") + len("total=")  : total_str.find("max=")]
                total_time += int(total_str)
            # Get all compute times for function found
            if "total=" in lines[j] and "elapsed time" in lines[j]:
                compute_str = lines[j].replace(' ', '').replace('\n', '').replace(',', '')
                compute_str = compute_str[compute_str.find("total=") + len("total=")  : compute_str.find("max=")]
                compute_time += int(compute_str)
            # Get all data times for function found.
            if "total=" in lines[j] and "device time" in lines[j]:
                data_str = lines[j].replace(' ', '').replace('\n', '').replace(',', '')
                data_str = data_str[data_str.find("total=") + len("total=")  : data_str.find("max=")]
                data_time += int(data_str)

        # For debug purposes.
        # print("Printing Stats for file", filepath)
        # print("The function found is", function_name)
        # print("Total time:", total_time * 1e-6, "s")
        # print("Total time:", compute_time * 1e-6, "s")
        # print("Total time:", data_time * 1e-6, "s")

        # Crete dictionary with collected info.
        timing_dict = {}
        timing_dict["File Path"] = filepath
        timing_dict["Function Name"] = function_name
        timing_dict["Total Time"] = total_time * 1e-6
        timing_dict["Compute Time"] = compute_time * 1e-6
        timing_dict["Data Time"] = data_time * 1e-6

        # Append info to dataframe
        temp_df = pd.DataFrame([timing_dict])
        timing_df = pd.concat([timing_df, temp_df], ignore_index=True)

    # Save info found to a .csv file.\
    timing_df.to_excel(results_folder + '/openacc_timing.xlsx', index=False)
