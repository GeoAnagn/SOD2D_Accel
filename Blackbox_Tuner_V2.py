import os
import sys
import json
import shutil
import logging
import opentuner
import pandas as pd
from Functions import file_management_V2, config_checks_V2
from Functions.Parsers import openacc_timing_data_parser
from opentuner import Result, IntegerParameter, ConfigurationManipulator, MeasurementInterface

# Set up logging for the tuner
def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)  # Set up logger for the tuner

    return logger

def load_configuration():
    config_path = "JSONs/Blackbox_Info.json"  # Path to the configuration JSON
    with open(config_path, 'r') as f:
        opentuner_info = json.load(f)  # Load the configuration data
    
    return opentuner_info

 # Set up paths for various directories
def setup_paths():
    paths = file_management_V2.path_definitions()  # Get directory paths from file_management
    
    return paths[0], paths[1], paths[2], paths[3]

# Load the results dataframe from the tuner results directory
def load_results_df(tuner_res_path, opentuner_info):
    results_path = f'{tuner_res_path}/results.xlsx'
    try:
        results_df = pd.read_excel(results_path)  # Load results dataframe
        config_counter = len(results_df)  # Update the configuration counter
    except Exception as e:
        results_df = pd.DataFrame(columns=opentuner_info['dataframe_columns'])  # Create an empty dataframe if no results exist
        config_counter = 0
    
    return results_df, config_counter

class GccFlagsTuner(MeasurementInterface):
    logger = setup_logging() # Set up logging for the tuner
    opentuner_info = load_configuration() # Load the configuration from JSON
    # Paths to the example, original results, tuner results, configuration results directories
    example_path, og_res_path, tuner_res_path, res_config_path = setup_paths()
    # Results Dataframe and Counter for tested configurations
    results_df, config_counter = load_results_df(tuner_res_path, opentuner_info)
    repetitions_counter = 0 # Counter for repeated configurations   
    
    # Environment variable for OpenACC timing analysis.
    os.environ['PGI_ACC_TIME'] = '1'
    # Set id ordered based on GPU BUS ID. 
    # Run nvidia-smi to see how GPUs will be ordered and pick your poison.
    os.environ['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'
    # Set appropriate GPU id.
    os.environ['CUDA_VISIBLE_DEVICES'] = opentuner_info['gpu_ids']

    # Initialize the ConfigurationManipulator and add parameters based on their type
    def manipulator(self):
        manipulator = ConfigurationManipulator()
        for parameter in self.opentuner_info['parameters']:
            if parameter["type"] == "integer":
                manipulator.add_parameter(IntegerParameter(parameter["name"], parameter["min"], parameter["max"]))  # Add integer parameters to manipulator
            else:
                self.logger.warning(f"Unsupported parameter type: {parameter['type']}")  # Log a warning for unsupported parameter types
        return manipulator
    
    # Set environment variables based on the configuration values
    def set_env_variables(self, cfg):
        self.logger.info("Setting variable values...")
        for param_name, param_value in cfg.items():
            multiplier = next(p['multiplier'] for p in self.opentuner_info["parameters"] if p["name"] == param_name)  # Find the multiplier for the parameter
            os.environ[param_name] = str(param_value * multiplier)  # Set environment variables based on the configuration values

    # Create a dictionary with the results of the execution
    def result_dict(self, execute_result, cfg, fail):
        results_dict = {}
        for param_name, param_value in cfg.items():
            multiplier = next(p['multiplier'] for p in self.opentuner_info["parameters"] if p["name"] == param_name)  # Find the multiplier for the parameter
            results_dict[param_name] = str(param_value * multiplier)  # Add parameter values to the results dictionary
        results_dict['time'] = 9999999 if fail else execute_result['time']  # Add execution time to the results dictionary, set to a large value if execution failed
        return results_dict
    
    # Main function for running the tuner
    def run(self, desired_result, input, limit):
        cfg = desired_result.configuration.data  # Get the configuration to test from OpenTuner

        if self.config_counter == self.opentuner_info["configs_to_check"]:
            self.finish_execution()
            sys.exit()

        if self.repetitions_counter == self.opentuner_info["repetitions"]:
            self.finish_execution()
            sys.exit()

        if not config_checks_V2.check_existing_configs(self.config_counter, cfg, self.res_config_path):
            self.repetitions_counter = 0
            self.set_env_variables(cfg)  # Set environment variables for the current configuration

            self.logger.info("Executing Sod2d Application")
            execute_cmd = f'cd {self.example_path} &&'
            execute_cmd += ' mpirun --allow-run-as-root --mca coll ^hcoll'
            execute_cmd += f' -np {self.opentuner_info["rank_num"]} {self.opentuner_info["sod2d_path"]} > openacc_timing.txt 2>&1'
            
            execute_result = self.call_program(execute_cmd)  # Execute the Sod2d Application
            try:
                assert execute_result['returncode'] == 0  # Check if execution was successful
                self.logger.info("Successful execution of Sod2d Application")
                self.handle_successful_execution(cfg, execute_result)  # Handle successful execution

                return Result(time=execute_result['time'])
            except AssertionError:
                self.logger.error("Execution of Sod2d Application Failed.")
                self.handle_failed_execution(cfg, execute_result)  # Handle failed execution

                return Result(time=execute_result['time'])
        else:
            self.repetitions_counter += 1
            self.logger.info("Configuration already tested, using cached result")
            cached_result = self.get_cached_result(cfg)
            return Result(time=cached_result)
        
    # Handle successful execution by saving results and moving files
    def handle_successful_execution(self, cfg, execute_result):
        os.system('cd ..')
        results_folder = file_management_V2.create_results_folder(self.res_config_path, self.config_counter)
        
        self.logger.info(f'Saving current config file at {results_folder}/config.json')
        self.manipulator().save_to_file(cfg, f'{results_folder}/config.json')

        self.logger.info(f'Saving results file at {results_folder}/results.json')
        with open(f'{results_folder}/results.json', 'w') as results:
            json.dump(execute_result['time'], results)

        file_management_V2.rm_files(self.example_path)
        file_management_V2.move_results(self.example_path, results_folder)
        openacc_timing_data_parser.parser(results_folder)

        results_dict = self.result_dict(execute_result, cfg, False)
        df_dictionary = pd.DataFrame([results_dict])
        self.results_df = pd.concat([self.results_df, df_dictionary], ignore_index=True)
        self.results_df.to_excel(f'{self.tuner_res_path}/results.xlsx', index=False)
        self.config_counter += 1
        self.logger.info("\n-----------------------------------------------------------\n")

    # Handle failed execution by saving error files and moving files
    def handle_failed_execution(self, cfg, execute_result):
        os.system('cd ..')
        results_folder = file_management_V2.create_results_folder(self.res_config_path, self.config_counter)

        self.logger.info(f'Saving current config file at {results_folder}/config.json')
        self.manipulator().save_to_file(cfg, f'{results_folder}/config.json')

        error_file = open(f'{results_folder}/error.txt', 'wb')
        error_file.write(execute_result['stderr'])
        error_file.close()

        file_management_V2.rm_files(self.example_path)
        file_management_V2.move_results(self.example_path, results_folder)
        openacc_timing_data_parser.parser(results_folder)

        results_dict = self.result_dict(execute_result, cfg, True)
        df_dictionary = pd.DataFrame([results_dict])
        self.results_df = pd.concat([self.results_df, df_dictionary], ignore_index=True)
        self.results_df.to_excel(f'{self.tuner_res_path}/results.xlsx', index=False)
        self.config_counter += 1
        self.logger.info("\n-----------------------------------------------------------\n")

    # Get the cached result for a configuration
    def get_cached_result(self, cfg):
        existing_cfg = pd.read_excel(f'{self.tuner_res_path}/results.xlsx', index_col=0)
        for param_name, param_value in cfg.items():
            multiplier = next(p['multiplier'] for p in self.opentuner_info["parameters"] if p["name"] == param_name)  # Find the multiplier for the parameter
            existing_cfg = existing_cfg.loc[existing_cfg[param_name] == param_value * multiplier]  # Locate the cached result for the configuration
        
        results_dict = existing_cfg.iloc[0].to_dict() 
        minimization_value = 1
        for param in results_dict:
            minimization_value = minimization_value * results_dict[param]
            
        return  minimization_value # Return the minimization value for the configuration
    
    # Finish the execution by moving the tuner results directory
    def finish_execution(self):
        destination_path = "Archive/Blackbox_Analysis/Modified_Folder"
        os.makedirs(os.path.dirname(destination_path), exist_ok=True)
        shutil.move(f'{self.tuner_res_path}/', destination_path)  # Move the tuner results directory to the archive
        self.logger.info("Checked all requested configs. Bye!")

if __name__ == '__main__':
    argparser = opentuner.default_argparser()
    GccFlagsTuner.main(argparser.parse_args())  # Start the tuner with command-line arguments
