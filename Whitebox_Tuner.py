import os
import sys
import json
import opentuner
import numpy as np
import pandas as pd
from opentuner import Result
from opentuner import IntegerParameter
from opentuner import MeasurementInterface
from opentuner import ConfigurationManipulator
from Functions import file_management, config_checks, array_reader, error_calculation
from Functions.Parsers import openacc_timing_data_parser

class GccFlagsTuner(MeasurementInterface):
    # Enviroment variable for OpenACC timing analysis.
    os.environ['PGI_ACC_TIME'] = '1'

    # Set id orded based on GPU BUS ID. 
    # Run nvidia-smi to see how gpus will be ordered and pick your poison.
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"

    # Set Json File Path.
    json_path = "Whitebox_Info.json"
    opentuner_info = None

    # Read opentuner user parameters.
    with open(json_path, 'r') as f:
        opentuner_info = json.load(f)

    # Set appropriate gpu id.
    os.environ["CUDA_VISIBLE_DEVICES"] = opentuner_info['gpu_ids']

    # Get all paths to be used.
    paths = file_management.path_definitions()
    # Set Example Folder Path.
    example_path = paths[0]
    # Set Original Results Folder Path.
    og_res_path = paths[1]
    # Set Opentuner Results Folder Path.
    tuner_res_path = paths[2]
    # Set result configs master folder
    res_config_path = paths[3]

    # Counter for how many configs have been tested.
    config_counter = 0
    # Change the confing counter to match the number of configs already tested.
    try:
        config_counter += len(next(os.walk(res_config_path))[1])
    except:
        config_counter = 0

    # Counter for stopping execution.
    nothing_to_do_counter = 0

    try:
        # Read Csv to continue where left off
        results_df = pd.read_csv(tuner_res_path + '/results.csv', index_col=0)
    except Exception as e:
        # Create Dataframe to store results.
        results_df = pd.DataFrame(columns=opentuner_info['dataframe_columns'])

    def manipulator(self):
        # Initialize Parameter Generator.
        manipulator = ConfigurationManipulator()
        
        # Add parameters based on their type.
        # TODO: Add more types.
        for parameter in self.opentuner_info['parameters']:
            if parameter["type"] == "integer":
                manipulator.add_parameter(
                    IntegerParameter(parameter["name"], parameter["min"], parameter["max"])
                )
            else:
                print(parameter["type"] + " type is not supported.")

        return manipulator

    def set_env_variables(self, cfg):
        print("Setting variable values...")
        configuration = []
        for i in cfg:
            not_found = 1
            for parameter in self.opentuner_info["parameters"]:
                if parameter["name"] == i:
                    if 'gang_num' in parameter['name'] or 'vector_num' in parameter['name']:
                        configuration.append(cfg[i]*parameter["multiplier"])    
                    
                    os.environ[i] = str(cfg[i]*parameter["multiplier"])
                    not_found = 0
                    break
            
            if not_found:
                print("Parameter", i, "in configurator not found in json file provided")
                
        return configuration

    def result_dict(self, execute_result, cfg, fail):
        results_dict = {}
        for i in cfg:
            for parameter in self.opentuner_info["parameters"]:
                if parameter["name"] == i:
                    results_dict[i] = str(cfg[i]*parameter["multiplier"])
                    
        results_dict["Mean Squared Calc Time"] = execute_result["Mean Squared Calc Time"]
        results_dict["Mean Squared Rmass Error"] = execute_result["Mean Squared Rmass Error"]
        results_dict["Mean Squared Rmom Error"] = execute_result["Mean Squared Rmom Error"]
        results_dict["Mean Squared Rener Error"] = execute_result["Mean Squared Rener Error"]

        return results_dict

    def run(self, desired_result, input, limit):
        cfg = desired_result.configuration.data

        # Check if nothing to do counter reaches the desired amount end the execution
        if self.nothing_to_do_counter == self.opentuner_info["program_end"]:
            print("Nothing more to do. Bye!")
            sys.exit()

        # Check if current config has been tested
        if not config_checks.check_existing_configs(self.config_counter, cfg, self.res_config_path):
            # Reinitialize end counter
            self.nothing_to_do_counter = 0

            # Create Call list
            func_call_info_list = []
            configuration = []
            
            # Setting environment variables
            configuration = self.set_env_variables(cfg)
            print(configuration)
            
            data_index = self.opentuner_info['data_index']
            
            for data_folder in data_index:
                # Set function call data folder.
                data_folder_path = self.opentuner_info["data_path"] + '/data_' + str(data_folder) + '/' + self.opentuner_info["func_name"]
                os.environ['folder_path'] = data_folder_path
                
                # Debug Info.
                print("Testing config N." + str(self.config_counter) + " at function call N." + str(data_folder))
                
                # Execute function.
                print("Executing " + self.opentuner_info["func_name"] + " function")
                execute_cmd = 'cd '+ self.example_path + ' &&'
                execute_cmd += ' mpirun --allow-run-as-root --mca coll ^hcoll'
                execute_cmd += ' -np ' + self.opentuner_info["rank_num"] + ' '
                execute_cmd += self.opentuner_info["func_path"] + self.opentuner_info["func_ver"] + self.opentuner_info["func_exec"] + '> openacc_timing.txt 2>&1'
                execute_result = self.call_program(execute_cmd)
                try:
                    assert execute_result['returncode'] == 0

                    # Parse OpenAcc Timing Analysis.
                    openacc_timing_data_parser.parser(self.example_path)
                    openacc_df = pd.read_csv(self.example_path + "/openacc_timing.csv")
                    
                    # TODO: Add values to a dataframe.
                    print("Calc Time:", openacc_df['Compute Time'].iloc[0])
                    

                    shape_variables = array_reader.initialize_shape_variables(data_folder_path)
                    # Calculate Rmass Error.
                    Rmass_error = error_calculation.error_calc(data_folder_path, self.example_path, 'Rmass', shape_variables)
                    print("Rmass Error:", Rmass_error)

                    # Calculate Rmass Error.
                    Rmom_error = error_calculation.error_calc(data_folder_path, self.example_path, 'Rmom', shape_variables)
                    print("Rmom Error:", Rmom_error)
                    
                    # Load Rener(npoin) array of function full_convec_ijk or full_diffusion_ijk.
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

            # Create new config results folder.
            results_folder = file_management.create_results_folder(self.res_config_path, self.config_counter)

            # Save config file
            print('Saving current config file at ' + results_folder + '/config.json')
            self.manipulator().save_to_file(cfg, results_folder +'/config.json')

            # Create dataframe with all call info
            func_call_info_df = pd.DataFrame(func_call_info_list)
            
            # Calculate Mean Squared Calc Time and Errors
            func_call_info_df['Squared Calc Time'] = np.power((func_call_info_df["Calc Time"]),2)
            func_call_info_df['Squared Rmass Error'] = np.power((func_call_info_df["Rmass Error"]),2)
            func_call_info_df['Squared Rmom Error'] = np.power((func_call_info_df["Rmom Error"]),2)
            func_call_info_df['Squared Rener Error'] = np.power((func_call_info_df["Rener Error"]),2)
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

            # Save results file
            print('Saving results file at ' + results_folder + '/results.json')
            with open(results_folder + '/results.json', "w") as results:
                # TODO: Give Appropriate value
                json.dump(mean_squared_time , results)

            # Merge results dataframe with new data 
            results_dict = self.result_dict(execute_result=result, cfg = cfg, fail = 0)

            # Append it to the results dataframe
            df_dictionary = pd.DataFrame([results_dict])
            self.results_df = pd.concat([self.results_df, df_dictionary], ignore_index=True)

            # Save results to a csv file
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
  