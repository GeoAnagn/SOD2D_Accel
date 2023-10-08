import os
import sys
import json
import logging
import opentuner
import numpy as np
import pandas as pd
from opentuner import Result
from opentuner import IntegerParameter
from opentuner import MeasurementInterface
from opentuner import ConfigurationManipulator
from Functions import file_management_V2, config_checks_V2, array_reader_V2, error_calculation_V2
from Functions.Parsers import openacc_timing_data_parser_V2

class GccFlagsTuner(MeasurementInterface):
    def __init__(self, *pargs, **kwargs):
        super(GccFlagsTuner, self).__init__(*pargs, **kwargs)
        self.load_configuration()  # Load the configuration from JSON
        self.setup_logging()  # Set up logging for the tuner
        self.setup_paths()  # Set up paths for various directories
        self.config_counter = 0  # Counter for tested configurations
        self.repetitions_counter = 0  # Counter for repeated configurations
        self.load_results_df()  # Load the results dataframe
        
    def load_configuration(self):
        config_path = "JSONs/Blackbox_Info.json"  # Path to the configuration JSON
        with open(config_path, 'r') as f:
            self.opentuner_info = json.load(f)  # Load the configuration data
            
        os.environ['PGI_ACC_TIME'] = '1'
        os.environ['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'
        os.environ["CUDA_VISIBLE_DEVICES"] = self.opentuner_info['gpu_ids']
        
    # Set up logging for the tuner
    def setup_logging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)  # Set up logger for the tuner

    # Set up paths for various directories
    def setup_paths(self):
        paths = file_management_V2.path_definitions()  # Get directory paths from file_management
        self.example_path = paths[0]  # Path to the example directory
        self.og_res_path = paths[1]  # Path to the original results directory
        self.tuner_res_path = paths[2]  # Path to the tuner results directory
        self.res_config_path = paths[3]  # Path to the configuration results directory

    # Load the results dataframe from the tuner results directory
    def load_results_df(self):
        try:
            results_path = self.tuner_res_path
            self.results_df = pd.read_excel(f'{results_path}/results.xlsx', index_col=0)  # Load results dataframe
            self.config_counter += len(self.results_df.index)  # Update the configuration counter
        except Exception as e:
            self.results_df = pd.DataFrame(columns=self.opentuner_info['dataframe_columns'])  # Create an empty dataframe if no results exist

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
        results_dict['Mean Squared Calc Time'] = 9999999 if fail else execute_result['Mean Squared Calc Time']  # Add execution time to the results dictionary, set to a large value if execution failed
        results_dict['Mean Squared Rmass Error'] = 9999999 if fail else execute_result['Mean Squared Rmass Error']  # Add execution time to the results dictionary, set to a large value if execution failed
        results_dict['Mean Squared Rmom Error'] = 9999999 if fail else execute_result['Mean Squared Rmom Error']  # Add execution time to the results dictionary, set to a large value if execution failed
        results_dict['Mean Squared Rener Error'] = 9999999 if fail else execute_result['Mean Squared Rener Error']  # Add execution time to the results dictionary, set to a large value if execution failed
        return results_dict

    def run(self, desired_result, input, limit):
        cfg = desired_result.configuration.data

        if self.config_counter == self.opentuner_info["configs_to_check"]:
            self.finish_execution()
            sys.exit()

        if self.repetitions_counter == self.opentuner_info["repetitions"]:
            self.finish_execution()
            sys.exit()

        if not config_checks_V2.check_existing_configs(self.config_counter, cfg, self.res_config_path):
            self.repetitions_counter = 0
            self.set_env_variables(cfg)  # Set environment variables for the current configuration
            func_call_info_list = []
            
            data_index = self.opentuner_info['data_index']

            for data_folder in data_index:
                data_folder_path = f""
                os.environ['folder_path'] = data_folder_path
                
        if not config_checks.check_existing_configs(self.config_counter, cfg, self.res_config_path):
            
            
            for data_folder in data_index:
                data_folder_path = self.opentuner_info["data_path"] + '/data_' + str(data_folder) + '/' + self.opentuner_info["func_name"]
                os.environ['folder_path'] = data_folder_path
                
                print("Testing config N." + str(self.config_counter) + " at function call N." + str(data_folder))
                
                print("Executing " + self.opentuner_info["func_name"] + " function")
                execute_cmd = 'cd '+ self.example_path + ' &&'
                execute_cmd += ' mpirun --allow-run-as-root --mca coll ^hcoll'
                execute_cmd += ' -np ' + self.opentuner_info["rank_num"] + ' '
                execute_cmd += self.opentuner_info["func_path"] + self.opentuner_info["func_ver"] + self.opentuner_info["func_exec"] + '> openacc_timing.txt 2>&1'
                execute_result = self.call_program(execute_cmd)
                try:
                    assert execute_result['returncode'] == 0

                    openacc_timing_data_parser.parser(self.example_path)
                    openacc_df = pd.read_csv(self.example_path + "/openacc_timing.csv")
                    
                    print("Calc Time:", openacc_df['Compute Time'].iloc[0])
                    
                    shape_variables = array_reader.initialize_shape_variables(data_folder_path)
                    Rmass_error = error_calculation.error_calc(data_folder_path, self.example_path, 'Rmass', shape_variables)
                    print("Rmass Error:", Rmass_error)

                    Rmom_error = error_calculation.error_calc(data_folder_path, self.example_path, 'Rmom', shape_variables)
                    print("Rmom Error:", Rmom_error)
                    
                    Rener_error = error_calculation.error_calc(data_folder_path, self.example_path, 'Rener', shape_variables)
                    print("Rener Error:", Rener_error)
                    
                    func_call_info_list.append({
                        "Calc Time": openacc_df['Compute Time'].iloc[0],
                        "Rmass Error": Rmass_error,
                        "Rmom Error": Rmom_error,
                        "Rener Error": Rener_error
                    })
                    
                    print("\n-----------------------------------------------------------\n")
                    
                except AssertionError:
                    print("Execution of Sod2d Application Failed.")

                    func_call_info_list.append({
                        "Calc Time": 999999,
                        "Rmass Error": 999999,
                        "Rmom Error": 999999,
                        "Rener Error": 999999
                    })

                    print("\n-----------------------------------------------------------\n")

                    break

            results_folder = file_management.create_results_folder(self.res_config_path, self.config_counter)
            print('Saving current config file at ' + results_folder + '/config.json')
            self.manipulator().save_to_file(cfg, results_folder +'/config.json')

            func_call_info_df = pd.DataFrame(func_call_info_list)
            
            func_call_info_df['Squared Calc Time'] = np.power((func_call_info_df["Calc Time"]), 2)
            func_call_info_df['Squared Rmass Error'] = np.power((func_call_info_df["Rmass Error"]), 2)
            func_call_info_df['Squared Rmom Error'] = np.power((func_call_info_df["Rmom Error"]), 2)
            func_call_info_df['Squared Rener Error'] = np.power((func_call_info_df["Rener Error"]), 2)
            func_call_info_df.to_csv(results_folder + '/calls_results.csv')

            mean_squared_time = func_call_info_df['Squared Calc Time'].mean()
            mean_squared_Rmass_Error = func_call_info_df['Squared Rmass Error'].mean()
            mean_squared_Rmom_Error = func_call_info_df['Squared Rmom Error'].mean()
            mean_squared_Rener_Error = func_call_info_df['Squared Rener Error'].mean()
            result = {
                "Mean Squared Calc Time": mean_squared_time,
                "Mean Squared Rmass Error": mean_squared_Rmass_Error,
                "Mean Squared Rmom Error": mean_squared_Rmom_Error,
                "Mean Squared Rener Error": mean_squared_Rener_Error
            }

            print('Saving results file at ' + results_folder + '/results.json')
            with open(results_folder + '/results.json', "w") as results:
                json.dump(mean_squared_time , results)

            results_dict = self.result_dict(execute_result=result, cfg = cfg, fail = 0)

            df_dictionary = pd.DataFrame([results_dict])
            self.results_df = pd.concat([self.results_df, df_dictionary], ignore_index=True)

            self.results_df.to_csv(self.tuner_res_path + '/results.csv')

            self.config_counter += 1

            if self.opentuner_info["func_ver"] == '/original':
                exit(1)
            
            return Result(time = mean_squared_time * configuration[0] * configuration[1])

        else:
            configuration = []
            
            existing_cfg = pd.read_csv(self.tuner_res_path + '/results.csv', index_col=0)
            for i in cfg:
                for parameter in self.opentuner_info["parameters"]:
                    if parameter["name"] == i:
                        if 'gang_num' in parameter['name'] or 'vector_num' in parameter['name']:
                            configuration.append(cfg[i]*parameter["multiplier"])   
                            
                        existing_cfg = existing_cfg.loc[existing_cfg[i] == cfg[i] * parameter["multiplier"]]

            print(configuration)
            self.nothing_to_do_counter += 1

            print("\n-----------------------------------------------------------\n")
            print(existing_cfg.iloc[0]['Mean Squared Calc Time'])
            return Result(time= existing_cfg.iloc[0]['Mean Squared Calc Time'] * configuration[0] * configuration[1])

if __name__ == '__main__':
    argparser = opentuner.default_argparser()
    GccFlagsTuner.main(argparser.parse_args())
