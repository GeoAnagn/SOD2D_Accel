import os
import sys
import shutil

def create_results_folder(res_config_path: str, config_counter: int) -> str:
    """
    Create a folder based on the configuration that is currently being tested.

    res_config_path: Parent folder of config folders.

    config_counter: Counter for naming config folders.

    returns: Path of the current config folder.
    """

    # Define path of config folder to be created.
    new_folder = os.path.join(res_config_path, f'config_{config_counter}')
    # Make directory with the path provided above.
    print(f'Creating results folder at {new_folder}')
    os.mkdir(new_folder)

    return new_folder

def move_results(example_path: str, results_folder: str) -> None:
    """
    Move Sod2D execution results, logs, and debug files.

    example_path: Folder where Sod2D has been executed.

    results_folder: Folder where to move Sod2D execution files.
    """

    print(f'Moving simulation results at {results_folder}')
    # Define the starting substring of the files to be moved.
    substrings = ["results_", "sod2d_", "timer_", "analysis_", "openacc_", "Utilization", "h5diff"]
    # Get all filenames from example_path.
    files = [f for f in os.listdir(example_path) if os.path.isfile(os.path.join(example_path, f))]
    # Check if filename contains any of the substrings defined above,
    # and move them to the results_folder.
    for filename in files:
        if any(substring in filename for substring in substrings):
            os.rename(os.path.join(example_path, filename), os.path.join(results_folder, filename))

def dump_results(results_folder: str) -> None:
    """
    Dump Sod2D generated .h5 files and adding a "modified_" substring at the start.

    results_folder: Folder containing Sod2D generated .h5 files.
    """

    # Get all filenames from the results folder.
    files = [f for f in os.listdir(results_folder) if os.path.isfile(os.path.join(results_folder, f))]

    # Iterate through the filename list
    for filename in files:
        # If "results_" and .h5 exists in the filename we must dump to read it.
        if "results_" in filename and filename.endswith(".h5"):
            print(f"Dumping {filename} to modified_{filename.replace('.h5', '')}")
            # Execute the command to dump the .h5 file to a same-name file,
            # with "modified_" added to the start of the filename in the results folder.
            h5dump_cmd = f'cd {results_folder} && h5dump {filename} > modified_{filename.replace(".h5", "")}'
            os.system(h5dump_cmd)


def compare_h5(results_folder: str, og_res_path: str):
    """
    Compare unmodified Sod2D .h5 files with the modified Sod2D .h5 files.

    results_folder: The folder where the results, from the modified Sod2D, have been stored.

    og_res_path: The folder where the results, from the unmodified Sod2D, have been stored.

    remove_files: Remove files for space preservation (0 or 1).
    """
    # Get all filenames from the results folder.
    modified_files = [f for f in os.listdir(results_folder) if os.path.isfile(os.path.join(results_folder, f))]
    # Get all filenames from the original results folder.
    original_files = [f for f in os.listdir(og_res_path) if os.path.isfile(os.path.join(og_res_path, f))]
    for modified_filename in modified_files:
        for original_filename in original_files:
            if modified_filename == original_filename and ".h5" in modified_filename:
                # Diff the two files and store the result to a .txt file starting with "diff_"
                print(f"Comparing {original_filename} and {modified_filename}")
                diff_cmd = f'h5diff --delta=10e-8 {os.path.join(og_res_path, original_filename)}'
                diff_cmd += f' {os.path.join(results_folder, modified_filename)}'
                diff_cmd += f' >> {os.path.join(results_folder, "h5diff")}.txt 2>&1 '
                os.system(diff_cmd)
                print(f'Output Saved to {os.path.join(results_folder, "h5diff")}.txt')

def compare_results(results_folder: str, og_res_path: str, remove_files: int = 1) -> None:
    """
    Compare unmodified Sod2D results with the modified Sod2D results.

    results_folder: The folder where the results, from the modified Sod2D, have been stored.
    The dumped .h5 files in the results_folder must have the string "modified_" at the start.

    og_res_path: The folder where the results, from the unmodified Sod2D, have been stored.
    The dumped .h5 files in og_res_path must have the string "original_" at the start.

    remove_files: Remove files for space preservation (0 or 1).
    """

    # Get all filenames from the results folder.
    modified_files = [f for f in os.listdir(results_folder) if os.path.isfile(os.path.join(results_folder, f))]
    # Get all filenames from the original results folder.
    original_files = [f for f in os.listdir(og_res_path) if os.path.isfile(os.path.join(og_res_path, f))]

    for modified_filename in modified_files:
        # Check if file with the starting string "modified_" exists.
        if modified_filename.startswith("modified_"):
            for original_filename in original_files:
                # Check if file with the starting string "original_" exists.
                if original_filename.startswith("original_"):
                    # Strip the two filenames from the "modified_" and "original_" substrings,
                    # and check if they are the same.
                    if modified_filename.replace('modified_', '') == original_filename.replace('original_', ''):
                        # Diff the two files and store the result to a .txt file starting with "diff_"
                        print(f"Comparing {original_filename} and {modified_filename}")
                        diff_cmd = f'diff {os.path.join(og_res_path, original_filename)}'
                        diff_cmd += f' {os.path.join(results_folder, modified_filename)}'
                        diff_cmd += f' > {os.path.join(results_folder, modified_filename.replace("modified_", "diff_"))}.txt 2>&1 '
                        os.system(diff_cmd)
                        print(f'Output Saved to {os.path.join(results_folder, modified_filename.replace("modified_", "diff_"))}.txt')

                        # Remove the .h5 files and the dumped versions of them to save space.
                        if remove_files == 1:
                            os.remove(os.path.join(results_folder, modified_filename))
                            os.remove(os.path.join(results_folder, modified_filename.replace('modified_', '')))
                        break

def rm_files(results_folder: str) -> None:
    """
    Remove Sod2D generated .h5 files.

    results_folder: The folder where the results, from the modified Sod2D, have been stored.
    """

    # Get all filenames from the results folder.
    modified_files = [f for f in os.listdir(results_folder) if os.path.isfile(os.path.join(results_folder, f))]
    for modified_filename in modified_files:
        # Check if file with starting string "results_" exists and remove it.
        if modified_filename.startswith("results_"):
            os.remove(os.path.join(results_folder, modified_filename))

def path_definitions() -> list:
    """
    Initialize paths for OpenTuner to execute, save, and modify files.
    """

    # Path where the example file must be stored.
    example_path = 'Example'
    # Path where the unmodified Sod2D result files must be stored.
    og_res_path = 'Results/Original_Results'
    # Path where OpenTuner must store the results generated.
    tuner_res_path = 'Results/Tuner_Results'
    # Path where the different configuration results will be stored.
    res_config_path = 'Results/Tuner_Results/Configs'

    # Different checks for the above paths' existence.
    if not os.path.isdir(example_path):
        os.mkdir(example_path)
        print(f"Please transfer example files to the generated folder at {os.getcwd()}/Example")
        sys.exit()
    elif not os.listdir(example_path):
        print("Example directory is empty, please transfer example files and run original.py")
        sys.exit()

    if not os.path.isdir('./Results'):
        print("Results directory not found, please run original.py")
        sys.exit()
    else:
        if not os.path.isdir(og_res_path):
            print("Original results directory not found, please run original.py")
            sys.exit()
        elif not os.listdir(og_res_path):
            print("Original results directory is empty, please run original.py")
            sys.exit()

        if not os.path.isdir(tuner_res_path):
            os.mkdir(tuner_res_path)

        if not os.path.isdir(res_config_path):
            os.mkdir(res_config_path)

    return [example_path, og_res_path, tuner_res_path, res_config_path]
