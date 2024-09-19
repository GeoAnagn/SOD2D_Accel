# import json
# import sqlite3
# import pandas as pd
# from pathlib import Path

# def parser(results_folder: str):
#     """Parse Nsight Timing Analysis and generate summary reports."""
#     results_path = Path(results_folder)
#     sqlite_file = results_path / "report1.sqlite"
    
#     # Connect to the sqlite file.
#     try:
#         database = sqlite3.connect(sqlite_file)
#     except sqlite3.Error as e:
#         print(f"Error connecting to database: {e}")
#         return None

#     database.text_factory = bytes  # Fixes utf-8 error

#     # Read NVTX table into a DataFrame
#     nvtx_df = pd.read_sql_query("SELECT * FROM NVTX_EVENTS", database)

#     # Calculate execution time
#     execution_time = (nvtx_df["end"].max() - nvtx_df["start"].min()) / 1_000_000_000
#     json_filename = results_path / 'time.json'
#     with open(json_filename, 'w') as json_file:
#         json.dump({'time': execution_time}, json_file)

#     # Calculate total time for each NVTX range
#     nvtx_df["time"] = (nvtx_df["end"] - nvtx_df["start"]) / 1_000_000_000
#     nvtx_df['function_name'] = nvtx_df["text"].apply(lambda x: x.decode('utf-8').split(':')[0])

#     # Create a summary directory if it does not exist
#     summary_dir = results_path / 'Timing_Summary'
#     summary_dir.mkdir(parents=True, exist_ok=True)

#     # Group by function name and save to separate Excel files
#     function_groups = nvtx_df.groupby('function_name')

#     for function_name, group in function_groups:
#         # Calculate total time and percentage for each NVTX range within the function
#         Timing_Summary = group.groupby('text').agg(
#             Time=('time', 'sum')
#         ).reset_index()
#         Timing_Summary['Time %'] = (Timing_Summary['Time'] / execution_time) * 100
#         Timing_Summary['Nvtx Range'] = Timing_Summary['text'].apply(lambda x: x.decode('utf-8'))
#         Timing_Summary.drop(columns='text', inplace=True)

#         # Save the timing summary to file
#         function_filename = summary_dir / f'{function_name}.xlsx'
#         Timing_Summary.to_excel(function_filename, index=False)

#     # Create an additional summary for specified categories
#     categories = ["CPU-COMPUTE", "CPU-DATA", "GPU-COMPUTE", "GPU-DATA"]
#     category_df = nvtx_df[nvtx_df["text"].apply(lambda x: any(cat.encode('utf-8') in x for cat in categories))]

#     # Group by function name and sum the time for all categories
#     category_summary = category_df.groupby('function_name').agg(
#         Time=('time', 'sum')
#     ).reset_index()
#     category_summary['Time %'] = (category_summary['Time'] / execution_time) * 100

#     # Save the category summary to file
#     category_filename = summary_dir / 'Category_Summary.xlsx'
#     category_summary.to_excel(category_filename, index=False)

#     # Close the database
#     database.close()

#     return execution_time



import os
import json
import sqlite3
import pandas as pd
from pathlib import Path

def parser(results_folder: str):
    """Parse Nsight Timing Analysis and generate summary reports."""
    results_path = Path(results_folder)
    sqlite_file = results_path / "report1.sqlite"
    
    # Connect to the sqlite file.
    try:
        database = sqlite3.connect(sqlite_file)
    except sqlite3.Error as e:
        print(f"Error connecting to database: {e}")
        return None

    database.text_factory = bytes  # Fixes utf-8 error

    # Read NVTX table into a DataFrame
    nvtx_df = pd.read_sql_query("SELECT * FROM NVTX_EVENTS", database)

    # Calculate execution time
    execution_time = (nvtx_df["end"].max() - nvtx_df["start"].min()) / 1_000_000_000
    json_filename = results_path / 'time.json'
    with open(json_filename, 'w') as json_file:
        json.dump({'time': execution_time}, json_file)

    # Calculate total time for each NVTX range
    nvtx_df["time"] = (nvtx_df["end"] - nvtx_df["start"]) / 1_000_000_000
    nvtx_df['function_name'] = nvtx_df["text"].apply(lambda x: x.decode('utf-8').split(':')[0])

    # Create a summary directory if it does not exist
    summary_dir = results_path / 'Timing_Summary'
    summary_dir.mkdir(parents=True, exist_ok=True)

    # Group by function name and save to separate Excel files
    function_groups = nvtx_df.groupby('function_name')

    for function_name, group in function_groups:
        # Calculate total time and percentage for each NVTX range within the function
        Timing_Summary = group.groupby('text').agg(
            Time=('time', 'sum')
        ).reset_index()
        Timing_Summary['Time %'] = (Timing_Summary['Time'] / execution_time) * 100
        Timing_Summary['Nvtx Range'] = Timing_Summary['text'].apply(lambda x: x.decode('utf-8'))
        Timing_Summary.drop(columns='text', inplace=True)

        # Save the timing summary to file
        function_filename = summary_dir / f'{function_name}.xlsx'
        Timing_Summary.to_excel(function_filename, index=False)

    # Create an additional summary for specified categories
    categories = ["CPU-COMPUTE", "CPU-DATA", "GPU-COMPUTE", "GPU-DATA"]
    category_df = nvtx_df[nvtx_df["text"].apply(lambda x: any(cat.encode('utf-8') in x for cat in categories))].copy()

    # Group by function name and sum the time for all categories
    category_df['category'] = category_df['text'].apply(lambda x: x.decode('utf-8').split(':')[1].strip())
    category_summary = category_df.groupby(['function_name', 'category']).agg(
        Time=('time', 'sum')
    ).reset_index()
    category_summary = category_summary.pivot_table(index='function_name', columns='category', values='Time', fill_value=0).reset_index()

    # Ensure all necessary columns are present
    for cat in categories:
        if cat not in category_summary.columns:
            category_summary[cat] = 0.0

    # Calculate the total time and percentages for each category
    category_summary['Total Time'] = category_summary[categories].sum(axis=1)
    category_summary['Bound'] = category_summary.apply(lambda row: 'CPU' if (row['CPU-COMPUTE'] + row['CPU-DATA']) > (row['GPU-COMPUTE'] + row['GPU-DATA']) else 'GPU', axis=1)
    category_summary['CPU-COMPUTE %'] = (category_summary['CPU-COMPUTE'] / category_summary['Total Time']) * 100
    category_summary['CPU-DATA %'] = (category_summary['CPU-DATA'] / category_summary['Total Time']) * 100
    category_summary['GPU-COMPUTE %'] = (category_summary['GPU-COMPUTE'] / category_summary['Total Time']) * 100
    category_summary['GPU-DATA %'] = (category_summary['GPU-DATA'] / category_summary['Total Time']) * 100
    category_summary['Total Time %'] = (category_summary['Total Time'] / execution_time) * 100

    category_summary.drop(columns=categories, inplace=True)
    category_summary.drop(columns='Total Time', inplace=True)

    # Save the category summary to file
    category_filename = summary_dir / 'Category_Summary.xlsx'
    category_summary.to_excel(category_filename, index=False)

    # Close the database
    database.close()

    return execution_time

