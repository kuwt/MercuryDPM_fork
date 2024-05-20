"""!
 @file read_mercury_cg.py
 @brief This file contains functions to read and process MercuryCG .stat files.
 @details This file is a Python script used for reading and processing MercuryCG .stat files.
 The main function in this file is `mercury_cg()`, which takes in a MercuryCG .stat file and processes the data into a
 dictionary.
 The dictionary keys are variable names and the values are numpy arrays of the data.
 The `mercury_cg()` function has three parameters: `name`, `print_summary`, and `row_major_index_order`.
 The function returns a dictionary where the keys are variable names and the values are numpy arrays of the data.
 The file also contains several helper functions such as `get_raw_data()`, `get_dimension()`, `get_coordinates()`,
 `get_variables()`, `extend_variables()`, and `read_header1()`.
 These functions are used to assist in the processing of the data from the MercuryCG .stat file.
"""

import numpy as np
import pandas as pd

# check if the required versions of pandas and numpy are installed
# specify the required versions
required_pandas_version = "1.2.0"
required_numpy_version = "1.20.0"

# get the installed versions
installed_pandas_version = pd.__version__
installed_numpy_version = np.__version__

# check if the installed versions meet the requirements
if installed_pandas_version < required_pandas_version:
    raise Exception(f"This script requires pandas version {required_pandas_version}, but you have {installed_pandas_version} installed.")

if installed_numpy_version < required_numpy_version:
    raise Exception(f"This script requires numpy version {required_numpy_version}, but you have {installed_numpy_version} installed.")


def mercury_cg(name='Chain.stat', print_summary=True, row_major_index_order=False):
    """!
     @brief This function reads a MercuryCG .stat file and returns a dictionary of the data.
     @details This function reads a MercuryCG .stat file, processes the data, and returns a dictionary of the data.
     The dictionary keys are the variable names and the values are numpy arrays of the data.

     @param[in] name
     The name of the .stat file created by MercuryCG to read from. Default is 'Chain.stat'.

     @param[in] print_summary
     If True, print a summary of the data key shapes. Default is True.

     @param[in] row_major_index_order
     If True, the index order for spatial coordinates is xyz, otherwise it is zyx. Default is False.

     @return data:
     A dictionary where the keys are the variable names and the values are numpy arrays of the data.
    """

    # ignore division by zero errors. Division by zero values will be set to NaN.
    # This line is just to suppress the error in the terminal
    np.seterr(divide='ignore', invalid='ignore')

    raw, header1, header2 = get_raw_data(name)
    data = read_header1(header1)
    data['name'] = name
    print(name)
    print(header1)
    dim = get_dimension(data, raw, row_major_index_order)
    data = get_coordinates(data, raw, dim, row_major_index_order)
    data = get_variables(data, raw, header2, dim, row_major_index_order)
    data = extend_variables(data)
    if print_summary:
        for key, value in data.items():
            print(f"{key}: {np.array(value).shape}")
    return data


def get_raw_data(name):
    """!
     @brief This function reads raw data from a file.
     @details This function reads raw data from a file and returns the raw data and the headers.

     @param[in] name
     The name of the file to read from.

     @return raw, header1, header2:
     The raw data of the .stat file, first header with the MercuryCG parameters, and second header with field names.
    """

    with open(name, 'r') as file:
        lines = file.readlines()

    # Extract headers
    header1 = lines[0].strip()
    if len(lines) > 1:
        header2 = lines[1].strip().split()
    else:
        header2 = None

    # Read data
    data_lines = lines[2:]
    raw = pd.DataFrame([(float(x) for x in line.strip().split()) for line in data_lines])

    return raw, header1, header2


def get_dimension(data, raw, row_major_index_order):
    """!
     @brief This function determines the dimension of the data.
     @details This function calculates the dimension of the data based on the number of data points and the row major
     index order.

     @param[in] data
     A dictionary where the keys are the variable names and the values are numpy arrays of the data.

     @param[in] raw
     A pandas DataFrame containing the raw data.

     @param[in] row_major_index_order
     A boolean value indicating whether the index order is row major. If True, the index order is xyz, otherwise it is
     zyx.

     @return dim:
     A list representing the dimension of the data for spatial coordinates and time if applicable.
    """

    if row_major_index_order:
        dim = data['n'] + [int(len(raw) / np.prod(data['n']))]
    else:
        dim = list(reversed(data['n'])) + [int(len(raw) / np.prod(data['n']))]
    dim = [d for d in dim if d != 1]
    dim.append(1) if not dim else None
    return dim


def get_coordinates(data, raw, dim, row_major_index_order):
    """!
     @brief This function gets the spatial coordinates and time of the data.
     @details This function extracts the spatial coordinates and time from the raw data and adds them to the data
     dictionary.

     @param[in,out] data
     A dictionary where the keys are the variable names and the values are numpy arrays of the data.

     @param[in] raw
     A pandas DataFrame containing the raw data.

     @param[in] dim
     A list representing the dimension of the data for spatial coordinates and time if applicable.

     @param[in] row_major_index_order
     A boolean value indicating whether the index order is row major. If True, the index order is xyz, otherwise it is
     zyx.

     @return data:
     The data dictionary updated with the spatial coordinates and time if applicable.
    """

    order = "C" if row_major_index_order else "F"
    data['t'] = np.reshape(raw.iloc[:, 0], dim, order=order)
    coord = 'xyz'
    coord = ''.join(c for c, n in zip(coord, data['n']) if n != 1)
    for i in range(len(coord)):
        data[coord[i]] = np.reshape(raw.iloc[:, i + 1], dim, order=order)
    return data


def get_variables(data, raw, names, dim, row_major_index_order):
    """!
     @brief This function gets the variables of the data from the .stat file.
     @details This function extracts the variables from the raw data and adds them to the data dictionary.

     @param[in] data
     A dictionary where the keys are the variable names and the values are numpy arrays of the data.

     @param[in] raw
     A pandas DataFrame containing the raw data.

     @param[in] names
     A list of strings representing the field names of the input from a MercuryCG .stat file.

     @param[in] dim
     A list representing the dimension of the data for spatial coordinates and time if applicable.

     @param[in] row_major_index_order
     A boolean value indicating whether the index order is row major. If True, the index order is xyz, otherwise it is
     zyx.

     @return data:
     The data dictionary updated with the variables as new keys of the dictionary.
    """

    order = "C" if row_major_index_order else "F"
    is_variable = [':' in name for name in names]
    for i in range(len(names)):
        if is_variable[i]:
            colon = names[i].find(':')
            hyphen = names[i].find('-')
            # if no hyphen is found
            if hyphen == -1:
                a = np.NaN
                b = np.NaN
            else:
                a = int(names[i][:hyphen])
                b = int(names[i][hyphen + 1:colon])
            n = names[i][colon + 1:]
            if hyphen == -1:
                a = int(names[i][:colon])
                data[n] = np.reshape(raw.iloc[:, a - 1], dim, order=order)
            elif b - a == 2:
                data[n + 'X'] = np.reshape(raw.iloc[:, a - 1], dim, order=order)
                data[n + 'Y'] = np.reshape(raw.iloc[:, a], dim, order=order)
                data[n + 'Z'] = np.reshape(raw.iloc[:, a + 1], dim, order=order)
            elif b - a == 5 and n != 'ParticleSize':
                data[n + 'XX'] = np.reshape(raw.iloc[:, a - 1], dim, order=order)
                data[n + 'XY'] = np.reshape(raw.iloc[:, a], dim, order=order)
                data[n + 'XZ'] = np.reshape(raw.iloc[:, a + 1], dim, order=order)
                data[n + 'YY'] = np.reshape(raw.iloc[:, a + 2], dim, order=order)
                data[n + 'YZ'] = np.reshape(raw.iloc[:, a + 3], dim, order=order)
                data[n + 'ZZ'] = np.reshape(raw.iloc[:, a + 4], dim, order=order)
            elif b - a == 8:
                data[n + 'XX'] = np.reshape(raw.iloc[:, a - 1], dim, order=order)
                data[n + 'XY'] = np.reshape(raw.iloc[:, a], dim, order=order)
                data[n + 'XZ'] = np.reshape(raw.iloc[:, a + 1], dim, order=order)
                data[n + 'YX'] = np.reshape(raw.iloc[:, a + 2], dim, order=order)
                data[n + 'YY'] = np.reshape(raw.iloc[:, a + 3], dim, order=order)
                data[n + 'YZ'] = np.reshape(raw.iloc[:, a + 4], dim, order=order)
                data[n + 'ZX'] = np.reshape(raw.iloc[:, a + 5], dim, order=order)
                data[n + 'ZY'] = np.reshape(raw.iloc[:, a + 6], dim, order=order)
                data[n + 'ZZ'] = np.reshape(raw.iloc[:, a + 7], dim, order=order)
            else:
                for j in range(a, b + 1):
                    data[n + str(j - a)] = np.reshape(raw.iloc[:, j - 1], dim, order=order)
    return data


def extend_variables(data):
    """!
     @brief This function extends the variables of the data dictionary by computing several fields from the standard
            fields.

     @param[in,out] data
     The data dictionary to extend the variables of.

     @return data:
     A dictionary with the data with possible extended variables computed from the standard variables.
    """
    if 'ParticleSize0' in data:
        data['ParticleNumber'] = data['ParticleSize0']
        mean = data['ParticleSize1'] / data['ParticleNumber']
        mom2 = data['ParticleSize2'] / data['ParticleNumber']
        mom3 = data['ParticleSize3'] / data['ParticleNumber']
        mom4 = data['ParticleSize4'] / data['ParticleNumber']
        mom5 = data['ParticleSize5'] / data['ParticleNumber']
        mom5 = mom5 - 5 * mean * mom4 + 10 * mean ** 2 * mom3 - 10 * mean ** 3 * mom2 + 4 * mean ** 5
        mom4 = mom4 - 4 * mean * mom3 + 6 * mean ** 2 * mom2 - 3 * mean ** 4
        mom3 = mom3 - 3 * mean * mom2 + 2 * mean ** 3
        mom2 = mom2 - mean ** 2
        std = np.sqrt(mom2)
        data['ParticleSize'] = {'Mean': mean, 'Variance': mom2, 'Skewness': mom3 / std ** 3,
                                'Kurtosis': mom4 / std ** 4, 'Hyperskewness': mom5 / std ** 5}

    if 'Density' in data:
        data['VelocityX'] = data['MomentumX'] / data['Density']
        data['VelocityY'] = data['MomentumY'] / data['Density']
        data['VelocityZ'] = data['MomentumZ'] / data['Density']
        data['VelocityX'][data['Density'] == 0] = 0
        data['VelocityY'][data['Density'] == 0] = 0
        data['VelocityZ'][data['Density'] == 0] = 0
        data['KineticStressXX'] = data['MomentumFluxXX'] - data['MomentumX'] * data['VelocityX']
        data['KineticStressXY'] = data['MomentumFluxXY'] - data['MomentumX'] * data['VelocityY']
        data['KineticStressXZ'] = data['MomentumFluxXZ'] - data['MomentumX'] * data['VelocityZ']
        data['KineticStressYY'] = data['MomentumFluxYY'] - data['MomentumY'] * data['VelocityY']
        data['KineticStressYZ'] = data['MomentumFluxYZ'] - data['MomentumY'] * data['VelocityZ']
        data['KineticStressZZ'] = data['MomentumFluxZZ'] - data['MomentumZ'] * data['VelocityZ']
        data['StressXX'] = data['ContactStressXX'] + data['KineticStressXX']
        data['StressXY'] = data['ContactStressXY'] + data['KineticStressXY']
        data['StressXZ'] = data['ContactStressXZ'] + data['KineticStressXZ']
        data['StressYX'] = data['ContactStressYX'] + data['KineticStressXY']
        data['StressYY'] = data['ContactStressYY'] + data['KineticStressYY']
        data['StressYZ'] = data['ContactStressYZ'] + data['KineticStressYZ']
        data['StressZX'] = data['ContactStressZX'] + data['KineticStressXZ']
        data['StressZY'] = data['ContactStressZY'] + data['KineticStressYZ']
        data['StressZZ'] = data['ContactStressZZ'] + data['KineticStressZZ']
        data['ContactPressure'] = (data['ContactStressXX'] + data['ContactStressYY'] + data['ContactStressZZ']) / 3
        data['KineticPressure'] = (data['KineticStressXX'] + data['KineticStressYY'] + data['KineticStressZZ']) / 3
        data['FluctuatingKineticEnergy'] = data['KineticStressXX'] + data['KineticStressYY'] + data['KineticStressZZ']
        data['KineticEnergy'] = data['MomentumFluxXX'] + data['MomentumFluxYY'] + data['MomentumFluxZZ']
        data['Pressure'] = (data['StressXX'] + data['StressYY'] + data['StressZZ']) / 3
        data['Temperature'] = data['KineticPressure'] / data['Density']
    return data


def read_header1(header1):
    """!
     @brief Parse the first header line of the input file.
     @details This function takes the first header line of the input file, splits it into
     separate words, and extracts the values associated with certain keys. The keys
     it looks for are 'n', 'width', 'timeMin', 'min', and 'max'. The values associated
     with these keys are then converted to floats and stored in a dictionary.

     @param[in] header1
     A string of the first header line of the input file.

     @return data:
     A dictionary where the keys are the keys found in the header line and the
     values are lists of floats associated with these keys. If a key is not found
     in the header line, its value in the dictionary is an empty list.
    """

    opt = header1.split()
    keys = ['n', 'width', 'timeMin', 'min', 'max']
    key_to_num_values = {'n': 3, 'width': 1, 'timeMin': 1, 'min': 3, 'max': 3}
    # List comprehension to get all the values that exist in the header
    values = [[float(opt[opt.index(key) + i]) for i in range(1, key_to_num_values[key] + 1)] if
              key in opt else [] for key in keys]

    data = {key: value if key in opt else [] for key, value in zip(keys, values)}
    data['n'] = [int(n) for n in data['n']]
    return data
