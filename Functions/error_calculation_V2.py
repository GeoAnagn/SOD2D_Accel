from Functions import array_reader_V2
from sklearn.metrics import mean_squared_error

def error_calc(og_data_folder: str, example_folder: str, variable_name: str, shape_variables: str) -> float:
    # Load original data
    og_file = array_reader_V2.data_read(og_data_folder, variable_name)
    # Transform to readable array
    og_array = array_reader_V2.array_decider(variable_name, shape_variables, og_file)

    # Load original data
    new_file = array_reader_V2.data_read(example_folder, variable_name)
    # Transform to readable array
    new_array = array_reader_V2.array_decider(variable_name, shape_variables, new_file)

    return mean_squared_error(og_array, new_array)