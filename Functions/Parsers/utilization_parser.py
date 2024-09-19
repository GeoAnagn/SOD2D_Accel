import pandas as pd

def parser(results_folder: str):
    """Parse Utilization Analysis."""
    # Get the path to the utilization file.
    utilization_file = results_folder + '/Utilization.csv'

    # Read the utilization file.
    columns = ['GPU ID', 
                'PWR (W)', # Power Draw (W) 
                'Core Temp (C)', # GPU Temperature (C)
                'Mem Temp (C)', # Memory Temperature (C)
                'SM Util (%)', # SM Utilization (%)
                'Mem Util (%)', # Memory Utilization (%)
                'Encoder Util (%)', # Encoder Utilization (%)
                'Decoder Util (%)', # Decoder Utilization (%)
                'JPG Util (%)', # JPG Encoder Utilization (%)
                'OFA Util (%)', # Optical Flow Accelerator Utilization (%)
                'MClk (MHz)', # Memory Clock (MHz)
                'PClk (MHz)', # Graphics Clock (MHz)
                'PWR Viol (%)', # Power Violation (%)
                'Therm Viol (Bool)', # Thermal Violation (%)
                'FB (MB)', # Frame Buffer Utilization (MB)
                'Bar1 (MB)', # BAR1 Memory Utilization (MB)
                'SBECC Errs', # Single Bit ECC Errors
                'DBECC Errs', # Double Bit ECC Errors
                'PCIe Errs', # PCIe Replay Errors
                'PCIe Rx (MB/s)', # PCIe Rx (MB/s)
                'PCIe Tx (MB/s)'] # PCIe Tx (MB/s)
    # Skip the first two rows and set the column names.
    utilization_df = pd.read_csv(utilization_file, skiprows=2, names=columns, delim_whitespace=True)
    
    # Replace '-' with 0.
    utilization_df = utilization_df.replace('-', 0)

    # Convert columns to numeric.
    columns = ['GPU ID', 
                'PWR Mean (W)', 'PWR Std (W)', 'Total PWR (J)',
                'Core Temp Mean (C)', 'Core Temp Std (C)', 'Core Temp Max (C)',
                'Mem Temp Mean (C)', 'Mem Temp Std (C)', 'Mem Temp Max (C)',
                'SM Util Mean (%)', 'SM Util Std (%)',
                'Mem Util Mean (%)', 'Mem Util Std (%)',
                'Encoder Util Mean (%)', 'Encoder Util Std (%)',
                'Decoder Util Mean (%)', 'Decoder Util Std (%)',
                'JPG Util Mean (%)', 'JPG Util Std (%)',
                'OFA Util Mean (%)', 'OFA Util Std (%)',
                'MClk Mean (MHz)', 'MClk Std (MHz)',
                'PClk Mean (MHz)', 'PClk Std (MHz)',
                'PWR Viol Mean (%)', 'PWR Viol Std (%)',
                'Therm Viol Sum', 
                'FB Mean (MB)', 'FB Std (MB)', 
                'Bar1 Mean (MB)', 'Bar1 Std (MB)',
                'SBECC Errs Sum',
                'DBECC Errs Sum',
                'PCIe Errs Sum',
                'PCIe Rx Mean (MB/s)', 'PCIe Rx Std (MB/s)',
                'PCIe Tx Mean (MB/s)', 'PCIe Tx Std (MB/s)']
    # Create a dataframe to store the results.
    results_df = pd.DataFrame(columns=columns)

    # For each GPU, calculate the mean and standard deviation of the metrics.
    for gpu in utilization_df['GPU ID'].unique():
        gpu_df = utilization_df.loc[utilization_df['GPU ID'] == gpu]
        dict = {
            'GPU ID': gpu,
            'PWR Mean (W)': gpu_df['PWR (W)'].mean(), 'PWR Std (W)': gpu_df['PWR (W)'].std(), 'Total PWR (J)': gpu_df['PWR (W)'].sum(),
            'Core Temp Mean (C)': gpu_df['Core Temp (C)'].mean(), 'Core Temp Std (C)': gpu_df['Core Temp (C)'].std(),
            'Mem Temp Mean (C)': gpu_df['Mem Temp (C)'].mean(), 'Mem Temp Std (C)': gpu_df['Mem Temp (C)'].std(),
            'SM Util Mean (%)': gpu_df['SM Util (%)'].mean(), 'SM Util Std (%)': gpu_df['SM Util (%)'].std(),
            'Mem Util Mean (%)': gpu_df['Mem Util (%)'].mean(), 'Mem Util Std (%)': gpu_df['Mem Util (%)'].std(),
            'Encoder Util Mean (%)': gpu_df['Encoder Util (%)'].mean(), 'Encoder Util Std (%)': gpu_df['Encoder Util (%)'].std(),
            'Decoder Util Mean (%)': gpu_df['Decoder Util (%)'].mean(), 'Decoder Util Std (%)': gpu_df['Decoder Util (%)'].std(),
            'JPG Util Mean (%)': gpu_df['JPG Util (%)'].mean(), 'JPG Util Std (%)': gpu_df['JPG Util (%)'].std(),
            'OFA Util Mean (%)': gpu_df['OFA Util (%)'].mean(), 'OFA Util Std (%)': gpu_df['OFA Util (%)'].std(),
            'MClk Mean (MHz)': gpu_df['MClk (MHz)'].mean(), 'MClk Std (MHz)': gpu_df['MClk (MHz)'].std(),
            'PClk Mean (MHz)': gpu_df['PClk (MHz)'].mean(), 'PClk Std (MHz)': gpu_df['PClk (MHz)'].std(),
            'PWR Viol Mean (%)': gpu_df['PWR Viol (%)'].mean(), 'PWR Viol Std (%)': gpu_df['PWR Viol (%)'].std(),
            'Therm Viol Sum': gpu_df['Therm Viol (Bool)'].sum(),
            'FB Mean (MB)': gpu_df['FB (MB)'].mean(), 'FB Std (MB)': gpu_df['FB (MB)'].std(),
            'Bar1 Mean (MB)': gpu_df['Bar1 (MB)'].mean(), 'Bar1 Std (MB)': gpu_df['Bar1 (MB)'].std(),
            'SBECC Errs Sum': gpu_df['SBECC Errs'].sum(),
            'DBECC Errs Sum': gpu_df['DBECC Errs'].sum(),
            'PCIe Errs Sum': gpu_df['PCIe Errs'].sum(),
            'PCIe Rx Mean (MB/s)': gpu_df['PCIe Rx (MB/s)'].mean(), 'PCIe Rx Std (MB/s)': gpu_df['PCIe Rx (MB/s)'].std(),
            'PCIe Tx Mean (MB/s)': gpu_df['PCIe Tx (MB/s)'].mean(), 'PCIe Tx Std (MB/s)': gpu_df['PCIe Tx (MB/s)'].std()
        }
        gpu_temp_df = pd.DataFrame(dict, index=[0])
        # Replace 0 with N/A on Temperature columns.
        gpu_temp_df['Mem Temp Mean (C)'] = gpu_temp_df['Mem Temp Mean (C)'].replace(0, 'N/A')
        gpu_temp_df['Mem Temp Std (C)'] = gpu_temp_df['Mem Temp Std (C)'].replace(0, 'N/A')

        # Append the results to the dataframe.
        results_df = pd.concat([results_df, gpu_temp_df], ignore_index=True)

    # Save the results to file.
    results_df.to_excel(results_folder + f'/GPUs_Utilization.xlsx', index=False)