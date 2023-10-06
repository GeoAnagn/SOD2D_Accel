import pandas as pd

def parse_openacc_timing(results_folder):
    """
    Parse .txt OpenACC Timing file.

    results_folder: Folder containing openacc_timing.txt file.
    """

    # Read timing file
    with open(results_folder + '/openacc_timing.txt') as timing_file:
        lines = timing_file.readlines()

    # Declare .csv columns and create a DataFrame.
    column_names = ["File Path", "Function Name", "Total Time", "Compute Time", "Data Time"]
    timing_df = pd.DataFrame(columns=column_names)

    # Create lists for filepath and function name indexing
    filepath_index = []
    function_names = []

    # Get file indexes for filepath and function name
    index = 0
    for path_line in lines:
        if "/" in path_line:
            filepath_index.append(index)
            function_names.append(lines[index + 1].strip())
        index += 1

    # Append the last line index to handle the end of the file
    filepath_index.append(len(lines))

    # Loop through each filepath found
    for i in range(len(filepath_index) - 1):
        # Declare filepath start and end.
        filepath_index_start = filepath_index[i]
        filepath_index_end = filepath_index[i + 1] - 1
        filepath = lines[filepath_index_start].strip()

        # Initialize total, compute, data times.
        compute_time = 0
        data_time = 0
        total_time = 0
        function_name = ""

        # Start searching for timing results.
        for j in range(filepath_index_start, filepath_index_end):
            line = lines[j].strip().replace(' ', '').replace(',', '')

            # Get function name.
            if line in function_names:
                function_name = line.split("NVIDIA")[0]

            # Get all times for function found
            if "total=" in line:
                total_str = line.split("total=")[1].split("max=")[0]
                total_time += int(total_str)

            # Get all compute times for function found
            if "total=" in line and "elapsedtime" in line:
                compute_str = line.split("total=")[1].split("max=")[0]
                compute_time += int(compute_str)

            # Get all data times for function found.
            if "total=" in line and "devicetime" in line:
                data_str = line.split("total=")[1].split("max=")[0]
                data_time += int(data_str)

        # Create dictionary with collected info.
        timing_dict = {
            "File Path": filepath,
            "Function Name": function_name,
            "Total Time": total_time * 1e-6,
            "Compute Time": compute_time * 1e-6,
            "Data Time": data_time * 1e-6,
        }

        # Append info to dataframe
        temp_df = pd.DataFrame([timing_dict])
        timing_df = pd.concat([timing_df, temp_df], ignore_index=True)

    # Save info found to an Excel file.
    timing_df.to_excel(results_folder + '/openacc_timing.xlsx', index=False)

