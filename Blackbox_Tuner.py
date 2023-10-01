import os
import sys
import json
import shutil
import opentuner
import pandas as pd
from opentuner import Result
from opentuner import IntegerParameter
from opentuner import MeasurementInterface
from opentuner import ConfigurationManipulator
from Functions import file_management, config_checks
from Functions.Parsers import openacc_timing_data_parser
# TODO: Change prints to logging.info

class GccFlagsTuner(MeasurementInterface):
    # Enviroment variable for OpenACC timing analysis.
    os.environ['PGI_ACC_TIME'] = '1'

    # Set id orded based on GPU BUS ID. 
    # Run nvidia-smi to see how gpus will be ordered and pick your poison.
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"

    # Set Json File Path.
    json_path = "JSONs\Blackbox_Info.json"
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
    repetitions_counter = 0
    
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

    # Export the configuration file parameters to environment variables.
    def set_env_variables(self, cfg):
        print("Setting variable values...")
        for i in cfg:
            not_found = 1
            for parameter in self.opentuner_info["parameters"]:
                if parameter["name"] == i:
                    os.environ[i] = str(cfg[i]*parameter["multiplier"])
                    not_found = 0
                    break
            
            if not_found:
                print("Parameter", i, "in configurator not found in json file provided")

    # Create a dictionary with the results of the execution.
    def result_dict(self, execute_result, cfg, fail):
        results_dict = {}
        for i in cfg:
            for parameter in self.opentuner_info["parameters"]:
                if parameter["name"] == i:
                    results_dict[i] = str(cfg[i]*parameter["multiplier"])

        if fail:
            results_dict['time'] = 9999999
        else:
            results_dict['time'] = execute_result['time']

        return results_dict

    def run(self, desired_result, input, limit):
        cfg = desired_result.configuration.data
        
        # Check if all configs have been tested and end the execution
        if self.config_counter == self.opentuner_info["configs_to_check"]:
            shutil.move(self.tuner_res_path, "/Archive/Blackbox_Analysis/Modified_Folder")
            print("Checked all requested configs. Bye!")
            sys.exit()
        
        # Check if repetitions limit has been reached
        if self.repetitions_counter == self.opentuner_info["repetitions"]:
            shutil.move(self.tuner_res_path, "/Archive/Blackbox_Analysis/Modified_Folder")
            print("Repetitions limit reached. Bye!")
            sys.exit()

        # Check if current config has been tested
        if not config_checks.check_existing_configs(self.config_counter, cfg, self.res_config_path):
            # Reinitialize end counter
            self.repetitions_counter = 0

            # Setting environment variables
            self.set_env_variables(cfg)

            # Execute Sod2d
            print("Executing Sod2d Application")
            execute_cmd = 'cd '+ self.example_path + ' &&'
            execute_cmd += ' mpirun --allow-run-as-root --mca coll ^hcoll'
            execute_cmd += ' -np ' + self.opentuner_info["rank_num"] + ' '
            execute_cmd += self.opentuner_info["sod2d_path"] + '> openacc_timing.txt 2>&1'
            execute_result = self.call_program(execute_cmd)
            try:
                assert execute_result['returncode'] == 0
                print("Succesful execution of Sod2d Application")

                # Create new config results folder
                results_folder = file_management.create_results_folder(self.res_config_path, self.config_counter)

                # Save config file
                print('Saving current config file at ' + results_folder + '/config.json')
                self.manipulator().save_to_file(cfg, results_folder +'/config.json')

                # Save results file
                print('Saving results file at ' + results_folder + '/results.json')
                with open(results_folder + '/results.json', "w") as results:
                    json.dump(execute_result['time'], results)

                # Move results to appropriate folder
                file_management.move_results(self.example_path, results_folder)

                # Parse OpenAcc Timing Analysis
                openacc_timing_data_parser.parser(results_folder)

                # TODO: ENABLE WHEN DECIDED
                # Dump .h5 files
                # file_management.dump_results(results_folder)

                # Compare Resutls
                # file_management.compare_results(results_folder, self.og_res_path, 1)

                # Remove .h5 files
                file_management.rm_files(results_folder)

                # Merge all results to one dictionary
                results_dict = self.result_dict(execute_result, cfg, 0)

                # Append it to the results dataframe
                df_dictionary = pd.DataFrame([results_dict])
                self.results_df = pd.concat([self.results_df, df_dictionary], ignore_index=True)

                # Save results to a csv file
                self.results_df.to_csv(self.tuner_res_path + '/results.csv')

                self.config_counter += 1

                print("\n-----------------------------------------------------------\n")

                return Result(time=execute_result['time'])

            except AssertionError:
                print("Execution of Sod2d Application Failed.")

                 # Create new config results folder
                results_folder = file_management.create_results_folder(self.res_config_path, self.config_counter)

                # Save config file
                print('Saving current config file at ' + results_folder +'/config.json')
                self.manipulator().save_to_file(cfg, results_folder +'/config.json')

                # Save error file
                error_file = open(results_folder + '/error.txt', "wb")
                error_file.write(execute_result['stderr'])
                error_file.close()

                # Move results to appropriate folder
                file_management.move_results(self.example_path, results_folder)

                # Remove .h5 files
                file_management.rm_files(results_folder)

                # Merge all results to one dictionary
                results_dict = self.result_dict(execute_result, cfg, 1)

                # Append it to the results dataframe
                df_dictionary = pd.DataFrame([results_dict])
                self.results_df = pd.concat([self.results_df, df_dictionary], ignore_index=True)

                # Save results to a csv file
                self.results_df.to_csv(self.tuner_res_path + '/results.csv')

                self.config_counter += 1

                print("\n-----------------------------------------------------------\n")

                return Result(time=9999999)
        else:
            # Find the information of configuration already tested and export the time.
            existing_cfg = pd.read_csv(self.tuner_res_path + '/results.csv', index_col=0)
            for i in cfg:
                for parameter in self.opentuner_info["parameters"]:
                    if parameter["name"] == i:
                        existing_cfg = existing_cfg.loc[existing_cfg[i] == cfg[i] * parameter["multiplier"]]

            self.repetitions_counter += 1

            print("\n-----------------------------------------------------------\n")

            return Result(time=existing_cfg.iloc[0]['time'])

if __name__ == '__main__':
    argparser = opentuner.default_argparser()
    GccFlagsTuner.main(argparser.parse_args())
  