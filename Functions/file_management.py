"""
Library for management of Sod2D generated files.
"""

import os
import sys

def create_results_folder(res_config_path: str, config_counter: str) -> str:
    """
    Create a folder based on config that is currently being tested.

    res_config_path: Parent folder of config folders.

    config_counter: Counter for naming config folders.

    returns: Path of current config folder.
    """

    # Define path of config folder to be created.
    new_folder = res_config_path + '/config_' + str(config_counter)
    # Make directory with path provided above.
    print('Creating results folder at ' + new_folder)
    os.mkdir(new_folder)

    return new_folder

def move_results(example_path: str, results_folder: str) -> None:
    """
    Move Sod2D execution results, logs, and debug files.

    example_path: Folder where Sod2D has been executed.

    results_folder: Folder where to move Sod2D execution files.
    """

    print('Moving simulation results at ' + results_folder)
    # Define the starting substring of the files to be moved.
    substrings = ["results_", "sod2d_", "timer_", "analysis_", "openacc_"]
    # Get all filenames from example_path.
    files = [f for f in os.listdir(example_path) if os.path.isfile(os.path.join(example_path, f))]
    # Check if filename contains any of the substrings defined above,
    # and move them to the resutls_folder.
    for filename in files:
        if any(substring in filename for substring in substrings):
            os.rename(example_path + '/' + filename, results_folder  + '/' +  filename)

def dump_results(results_folder: str) -> None:
    """
    Dump Sod2D genereted .h5 files and adding a "modified_" substring at the start.

    results_folder: Folder containing Sod2D generated .h5 files.
    """

    # Get all filenames from the results folder.
    files = [f for f in os.listdir(results_folder) if os.path.isfile(os.path.join(results_folder, f))]

    # Iterate filename list
    for filename in files:
        # If "results_" and .h5 exists in the filename we must dump in order to read it.
        if(("results_" in filename) and (".h5" in filename)):
            print("Dumping " + filename + " to modified_" + str(filename).replace('.h5', ''))
            # Execute the command to dump the .h5 file to a same name file,
            # with "modified_" added to the start of the filename at the results folder.
            h5dump_cmd = 'cd ' + results_folder + '&&'
            h5dump_cmd += ' h5dump ' + filename + ' > modified_' + str(filename).replace('.h5', '')
            os.system(h5dump_cmd)

def compare_results(results_folder: str, og_res_path: str, remove_files: int = 1) -> None:
    """
    Compare unmodified Sod2D results with the modified Sod2D results.

    results_folder: The folder where the results, from the modified Sod2D, have been stored.
    The dumped .h5 files in the results_folder, must have the string "modified_" at the start. (e.g. modified_results_AVG_cube_per64-1_1)
    
    og_res_path: The folder where the results, from the unmodified Sod2D, have been stored.
    The dumped .h5 files in the og_res_path, must have the string "original_" at the start. (e.g. original_results_AVG_cube_per64-1_1)
    
    remove_files: Remove files for space preservation. (0 or 1)
    """

    # Get all filenames from the results folder.
    modified_files = [f for f in os.listdir(results_folder) if os.path.isfile(os.path.join(results_folder, f))]
    # Get all filenames from the original results folder.
    original_files = [f for f in os.listdir(og_res_path) if os.path.isfile(os.path.join(og_res_path, f))]

    for modifided_filename in modified_files:
        # Check if file with starting string "modified_" exist.
        if "modified_" in modifided_filename:
            for original_filename in original_files:
                # Check if file with starting string "original_" exist.
                if "original_" in original_filename:
                    # Strip the two filenames from the "modified_" and "original_" substrings,
                    # and check if they are the same
                    if str(modifided_filename).replace('modified_', '') == str(original_filename).replace('original_', ''):

                        # Diff the two files and store the result to a .txt file starting with "diff_". (e.g diff_results_AVG_cube_per64-1_1)
                        print("Comparing " + original_filename + " and " + modifided_filename)
                        diff_cmd = 'diff'
                        diff_cmd += ' ' + og_res_path + '/' + original_filename
                        diff_cmd += ' ' + results_folder  + '/' + modifided_filename
                        diff_cmd += ' > ' + results_folder + '/' + str(modifided_filename).replace('modified', 'diff') + '.txt 2>&1 '
                        os.system(diff_cmd)
                        print('Output Saved to ' + results_folder + '/' + str(modifided_filename).replace('modified', 'diff') + '.txt')

                        # Remove the .h5 files and the dumped versions of them to save space.
                        if remove_files == 1:
                            os.remove(results_folder + '/' + modifided_filename)
                            os.remove(results_folder + '/' + str(modifided_filename).replace('modified_', '') + '.h5')
                        break

def rm_files(results_folder: str) -> None:
    """
    Remove Sod2D generated .h5 files.
    
    results_folder: The folder where the results, from the modified Sod2D, have been stored.
    """

    # Get all filenames from the results folder.
    modified_files = [f for f in os.listdir(results_folder) if os.path.isfile(os.path.join(results_folder, f))]
    for modifided_filename in modified_files:
        # Check if file with starting string "results_" exist and remove it.
        if "results_" in modifided_filename:
            os.remove(results_folder + '/' + modifided_filename)

def path_definitions() -> list:
    """
    Initialize paths for OpenTuner to execute, save and modify files.
    """

    # Path where the example file must be stored.
    example_path = './Example'
    # Path where the unmodified Sod2D result files must be stored.
    og_res_path = './Results/Original_Results'
    # Path where OpenTuner must store the results generated.
    tuner_res_path = './Results/Tuner_Results'
    # Path where the different configuration results will be stored.
    res_config_path = './Results/Tuner_Results/Configs'

    # Different checks for above paths existence.
    if not os.path.isdir(example_path):
        os.mkdir(example_path)
        print("Please transfer example files to generated folder at " + os.getcwd() + "/Example")
        sys.exit()
    elif not os.listdir(example_path):
        print("Example directory is empty please transfer example files and run original.py")
        sys.exit()

    if not os.path.isdir('./Results'):
        print("Results direcotry not found please run original.py")
        sys.exit()
    else:
        if not os.path.isdir(og_res_path):
            print("Original results directory not found please run original.py")
            sys.exit()
        elif not os.listdir(og_res_path):
            print("Original results directory is empty please run original.py")
            sys.exit()

        if not os.path.isdir(tuner_res_path):
            os.mkdir(tuner_res_path)

        if not os.path.isdir(res_config_path):
            os.mkdir(res_config_path)

    return [example_path, og_res_path, tuner_res_path, res_config_path]
