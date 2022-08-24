#!/usr/bin/env python
import sys
import os
from operator import itemgetter
from netCDF4 import Dataset
import numpy as np
import getopt
from datetime import datetime, timezone

'''
 NAME: GeoCSV_2_netCDF.py - a Python script to read a GeoCSV Earth model file (2D, 3D or mixed) and convert it 
       to the netCDF format.

 USAGE:

                            input GeoCSV Earth model file           
                            |       debug mode      
                            |       |   display headers only           
                            |       |   |            
       GeoCSV_2_netCDF_3D  -i FILE -d  -H

 HISTORY:
   2022-08-25 IRIS DMC Manoch: v2022.237 R2 release combines 2D and 3D scripts and also adds support for the 
                               models with projected coordinate systems.
   2022-08-18 IRIS DMC Manoch: v2022.230 updated print formats.
   2022-01-28 IRIS DMC Manoch: v2022.028 updated how it handles X and Y to avid issue with tags. Converted formats
                               to fstring.
   2022-01-26 IRIS DMC Manoch: v2022.026 addressed the issue caused by v2021.088 when the value of a variable is 
                               the same as the variable name.
   2021-03-39 IRIS DMC Manoch: v2021.088 added check for presence of the required columns.
   2021-02-01 IRIS DMC Manoch: v2021.032 fixed the bug that was trying to do a replace on a numeric value.
   2020-09-29 IRIS DMC Manoch: v2020.273 now supports location variables other than latitude and longitude.
   2020-06-16 IRIS DMC Manoch: v2020.168 made sure variable types are set based on VAR_DTYPE. Also, minor style updates
   2020-02-28 IRIS DMC Manoch: v2020.059 float global attribute values are outputted as float and not string. 
   2020-01-06 IRIS DMC Manoch: v2020.006 R1 History now includes the source file name
   2020-01-03 IRIS DMC Manoch: v2020.003 preserves the history and avoids mixing variable names with common
                               characters (like Qp and QpQs)
   2019-11-14 IRIS DMC Manoch: v2019.318 now retains order of variables
   2019-01-28 IRIS DMC Manoch: v2019.148 removed the extra '_' character  behind the coordinate variable parameter
                               names.
   2019-01-22 IRIS DMC Manoch: v2019.022 modified to fill in the missing points (if any) with nan, rather than zeros
   2018-10-25 IRIS DMC Manoch: expanded the error message R.0.5.2018.298
   2018-10-22 IRIS DMC Manoch: initial release r0.5.2018.295
   2018-10-17 IRIS DMC Manoch: fill_value and missing_value attributes are now supported
   2018-10-16 IRIS DMC Manoch: created version 2018.289
'''

script = os.path.basename(sys.argv[0])
version = 'v2022.237'
print(f'\n\n[INFO] {script} version {version}', flush=True)

debug = False

# The root _directory.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# The optional directory to store netCDF and GeoCSV files.
DATA_DIR = os.path.join(ROOT_DIR, 'data')

NETCDF_FORMAT = 'NETCDF3_CLASSIC'
geocsv_file_name = None
BASE_NAME = None
view_header = False
VAR_DTYPE = 'f4'  # set variables as 'f4' (32-bit floating point) or 'f8' (64-bit floating point) variables

# A dictionary of parameter tags.
tag_dict = dict()


def usage():
    print(f"\n{script} reads a 2D or 3D GeoCSV Earth model file and convert it to the netCDF format.\n\n")
    print(f" USAGE:\n\n",
          f"{script} -i FILE -d -H\n"
          f"\t-i [required] FILE is the input GeoCSV Earth model file\n",
          f"\t-d [default:{debug}] debug mode\n",
          f"\t-H [default:{view_header}] display headers only\n"
          )


def dot():
    """print a dot on screen"""

    sys.stdout.write('.')
    sys.stdout.flush()


def get_float(this_string):
    """convert a string to float, otherwise leave it as a string

    Keyword arguments:
    this_string: the string to convert

    Return values:
    this_value: converted value or the original string
    """

    try:
        if VAR_DTYPE == 'f4':
            this_value = np.float32(this_string)
        else:
            this_value = np.float64(this_string)

    except ValueError as e:
        this_value = this_string
    return this_value


def header_tags(par):
    cols = par['header']['value']
    tags = list()
    for col in cols:
        found = False
        for key in par:
            if 'column' in par[key]:
                if par[key]['column'] == col:
                    tags.append(key)
                    found = True
        if not found:
            print(f"\n[ERR] parameter with column '{col}' not found!")
            sys.exit(2)
    return tags


def get_attribute(par, var, att):
    """return an attribute's value

    Keyword arguments:
    par: the working netCDF dictionary of the key-value pairs found in the header
    var: variable to check
    att: attribute to find

    Return values:
    this_value: value of the attribute, if found, otherwise None
    """
    val = None
    if att in par[var]:
        val = par[var][att]
    return val


def str2float(sequence):
    """convert content of a sequence to float, otherwise leave them
    as string

    Keyword arguments:
    sequence: a sequence of values to convert

    Return values:
    a generator of the converted sequence
    """

    for item in sequence:
        try:
            yield float(item)
        except ValueError as e:
            yield item


def read_geocsv(file_name):
    """read a given GeoCSV file and output a dictionary of the header keys along with a sorted
    block of data

    Keyword arguments:
    file_name: the GeoCSV file to read

    Return values:
    header_params: a dictionary of the key-value pairs found in the header
    data: a text block of the body
    """
    fp = open(file_name, 'r')
    print(f'[INFO] Input GeoCSV file: {file_name}', flush=True)
    if not debug:
        dot()

    try:
        content = fp.read()
    except Exception as err:
        print(f'\n[ERR] failed to read the input file!\n{err}\n')
        fp.close()
        sys.exit(2)

    lines = content.split('\n')
    fp.close()

    # The header parameters.
    header_params = dict()
    data = list()
    found_header = False

    # Go through lines in the file.
    for line in lines:
        line = line.strip()
        # Ignore blank or comment lines.
        if not line or line == '#':
            continue
        # Metadata and comment lines.
        elif line.startswith('#'):
            # Remove the hash mark.
            line = line[1:]
            # Show the progress.
            if not debug:
                dot()
            # Split on the first ":" to get the key-value pair.
            key, value = line.split(':', 1)
            key = key.strip()
            value = value.strip()

            # If key has a "_" then it has two parts, split on the first "_".
            if '_' in key:
                main_key, aux_key = key.split('_', 1)
            else:
                main_key = key
                aux_key = 'value'
            if main_key not in header_params:
                header_params[main_key] = dict()
            # This is done in case there are spaces around :.
            header_params[main_key][aux_key] = value

        # Data lines.
        else:
            # The first data line that does not start with an # is the header.
            if not found_header:
                header_params['header'] = dict()
                header_params['header']['value'] = list()
                for head in line.split(header_params['delimiter']['value']):
                    header_params['header']['value'].append(head.strip())
                found_header = True
                continue
            # The rest are data lines.
            data.append(list(str2float(line.split(header_params['delimiter']['value']))))

    # Set the max dimension for the dataset.
    header_params['max_dimensions'] = {'value': 2}
    for key in header_params:
        if 'dimensions' in header_params[key]:
            if int(header_params[key]['dimensions']) == 3:
                header_params['max_dimensions']['value'] = 3
                break

    # The x and y columns are required.
    for axis in ('x', 'y'):
        if axis not in header_params:
            usage()
            print(f'\n[ERR] The "{axis}" parameter is required', flush=True)
            sys.exit()

    if 'z' not in header_params and header_params['max_dimensions']['value'] > 2:
        usage()
        print(f'\n[ERR] For the 3D parameters, the "z" parameter is required', flush=True)
        sys.exit()
    if not debug:
        dot()

    header = header_params['header']['value']
    z_column = None
    try:
        y_column = header.index(header_params['y']['column'])
        x_column = header.index(header_params['x']['column'])
        if 'z' in header_params:
            z_column = header.index(header_params['z']['column'])
    except Exception as ex:
        usage()
        print(f'\n\n[ERR] did not find all the required variables (y, x, '
              f'z) columns for a 3D model!\n{ex}')
        sys.exit(1)

    if not view_header:
        if z_column is not None:
            data.sort(key=itemgetter(z_column, y_column, x_column))
        else:
            data.sort(key=itemgetter(y_column, x_column))

    return header_params, data


def get_dimensions(column_data, header_params):
    """extract sorted unique values of the dimension columns

    Keyword arguments:
    column_data: matrix representing columns of a data matrix
    header_params: GeoCSV parameters

    Return values:
    y: sorted unique y values
    x: sorted unique x values
    z: sorted unique z values
    """

    try:
        header = header_params['header']['value']
        column = 'y'
        y_column = header.index(header_params['y']['column'])
        column = 'x'
        x_column = header.index(header_params['x']['column'])
        z_column = None
        if 'z' in header_params:
            column = 'z'
            z_column = header.index(header_params['z']['column'])
    except Exception as ex:
        print(
            f"\n[ERR] The 'column' field is required for all variables. Missing the column field for '{column}'\n{ex}")
        sys.exit(2)

    if not debug:
        dot()
    y = list(sorted(set(column_data[y_column])))
    if not debug:
        dot()
    x = list(sorted(set(column_data[x_column])))
    if not debug:
        dot()
    z = list()
    if z_column is not None:
        z = list(sorted(set(column_data[z_column])))

    if debug:
        print(f'\n\n[INFO] values:\n\n{header[y_column]}:{y},\n\n{header[x_column]}:{x},'
              f'\n\n{header[z_column]}:{z}', flush=True)
    else:
        dot()
    return y, x, z


def set_dimensions(this_dataset, header_params, x_coord_list, y_coord_list, z_values_list):
    """define variables' dimensions

    Keyword arguments:
    this_dataset: the working netCDF dataset
    header_params: GeoCSV parameters
    y_coord_list: list of unique y values from GeoCSV file
    x_coord_list: list of unique x values from GeoCSV file
    z_values_list: list of unique z values from GeoCSV file

    Return values:
    this_dataset: the input dataset with dimensions set
    """

    z_dim = None

    if header_params['max_dimensions']['value'] > 2:
        z_dim = this_dataset.createDimension(header_params['z']['column'], len(z_values_list))

    y_dim = this_dataset.createDimension(header_params['y']['column'], len(y_coord_list))
    x_dim = this_dataset.createDimension(header_params['x']['column'], len(x_coord_list))

    if debug:
        print(f"\n\n[INFO] dimensions:\n\n{header_params['y']['column']}:{y_dim},"
              f"\n{header_params['x']['column']}:{x_dim},\n{header_params['z']['column']}:{z_dim}",
              flush=True)
    else:
        dot()

    return this_dataset


def create_coordinate_variables(ds, par):
    """ create coordinate variables

    Keyword arguments:
    ds: the working netCDF dataset
    par: GeoCSV parameters

    Return values:
    dataset: the input dataset with coordinate variables set
    y: y variable
    x: x variable
    z: z variable
    """

    coordinate = dict()
    x_column = par['x']['column']
    y_column = par['y']['column']
    if par['max_dimensions']['value'] == 3:
        z_column = f"{par['z']['column']}"
    else:
        z_column = 'z'

    if par['max_dimensions']['value'] == 3:
        this_value = get_attribute(par, 'z', '_FillValue')
        if this_value:
            coordinate[z_column] = ds.createVariable(
                z_column, VAR_DTYPE, (z_column,), fill_value=this_value)
        else:
            coordinate[par['z']['column']] = ds.createVariable(
                z_column, VAR_DTYPE, (z_column,))

    this_value = get_attribute(par, 'y', '_FillValue')
    if this_value:
        coordinate[y_column] = ds.createVariable(
            y_column, VAR_DTYPE, (y_column,), fill_value=this_value)
    else:
        coordinate[y_column] = ds.createVariable(
            y_column, VAR_DTYPE, (y_column,))

    this_value = get_attribute(par, 'x', '_FillValue')
    if this_value:
        coordinate[x_column] = ds.createVariable(
            x_column, VAR_DTYPE, (x_column,), fill_value=this_value)
    else:
        coordinate[x_column] = ds.createVariable(
            x_column, VAR_DTYPE, (x_column,))

    if par['max_dimensions']['value'] == 2:
        column_list = ('y', 'x')
        coordinate['z'] = None
    else:
        column_list = ('y', 'x', 'z')
    for var in column_list:
        # If the variable parameter is not defined, the column will be selected.
        try:
            var_name = par[var]['variable']
        except Exception as ex:
            var_name = par[var]['column']

        for key in par[var]:
            if key not in ('column', 'variable', 'dimensions'):
                # let it default to default values
                if key != '_FillValue':
                    if key == 'missing_value':
                        setattr(coordinate[var_name], key, get_float(par[var][key]))
                    elif key == '_range':
                        this_range = list()
                        this_value = par[var][key].replace('[', '').replace(']', '').strip().split(',')
                        for item in this_value:
                            this_range.append(get_float(item))
                        setattr(coordinate[var_name], key, np.array(this_range, dtype=VAR_DTYPE))
                    else:
                        setattr(coordinate[var_name], key, par[var][key])

    if debug:
        print(f"\n\n[INFO] coordinate variables:\n\n{y_column}:"
              f"{coordinate[y_column]},\n{x_column}:"
              f"{coordinate[x_column]},\n{par['z']['column']}:"
              f"{coordinate[z_column]}", flush=True)
    else:
        dot()

    return ds, coordinate[y_column], coordinate[x_column], \
           coordinate[z_column]


def create_variables(ds, par, data, latitude_list, longitude_list, levels_list):
    """ create variables

    Keyword arguments:
    ds: the working netCDF dataset
    this_params" GeoCSV parameters
    data: block of data
    latitude_list: list of unique latitudes from GeoCSV file
    longitude_list: list of unique longitudes from GeoCSV file
    levels_list: list of unique levels from GeoCSV file

    Return values:
    ds: the input dataset with coordinate variables set
    all_vars: a dictionary of 3D variables
    """

    # var_list = list(set(header) - set([par['latitude_column'], par['longitude_column'],
    #                                   par['level_column']]))

    # We want to retain order of variables.
    # var_list = list(set(header) - set([par['latitude_column'], par['longitude_column']]))
    var_list = list()

    # The header is a list of all variables.
    for tag in header_tags(par):
        # Skip coordinate variables.
        if tag in ['x', 'y', 'z']:
            continue
        if 'dimensions' in par[tag]:
            z_column = 'z'
            if int(par[tag]['dimensions']) > 2:
                z_column = par['z']['column']
        else:
            print(f"\n[ERR] Missing dimension definition for '{tag}'")
            sys.exit(2)

        var_list.append(tag)
        """ We may need to explore the variable field to name the variable
        for key in par[var]:
            if 'column' in par[key]:
                if par[key]['column'] == key:
                    var_list.append(key)
        """
    all_vars = dict()
    var_values = dict()

    variable = dict()
    for var in var_list:
        three_d = False
        if 'dimensions' in par[var]:
            if int(par[var]['dimensions']) > 2:
                three_d = True
        else:
            three_d = True
        # Tag is used when user uses a different tag than the variable name. For now not used.
        tag = var
        # tag = tag_dict[var]
        this_value = get_attribute(par, tag, '_FillValue')
        if 'variable' in par[var]:
            variable[var] = par[var]['variable']
        else:
            variable[var] = par[var]['column']

        if this_value:
            if three_d:
                all_vars[variable[var]] = ds.createVariable(
                    variable[var], VAR_DTYPE, (par['z']['column'], par['y']['column'],
                                               par['x']['column']), fill_value=this_value)
            else:
                all_vars[variable[var]] = ds.createVariable(
                    variable[var], VAR_DTYPE, (par['y']['column'],
                                               par['x']['column']), fill_value=this_value)
        else:
            if three_d:
                all_vars[variable[var]] = ds.createVariable(
                    variable[var], VAR_DTYPE, (par['z']['column'], par['y']['column'],
                                               par['x']['column']))
            else:
                all_vars[variable[var]] = ds.createVariable(
                    variable[var], VAR_DTYPE, (par['y']['column'],
                                               par['x']['column']))

        # Create the data holders and fill in the points with nan.
        if int(par[var]['dimensions']) == 2:
            var_values[variable[var]] = np.empty((len(y_list), len(x_list)))
            var_values[f"{variable[var]}_done"] = np.empty((len(y_list), len(x_list)))
            var_values[f"{variable[var]}_done"][:] = np.nan
        else:
            var_values[variable[var]] = np.empty((len(z_list), len(y_list), len(x_list)))
        var_values[variable[var]][:] = np.nan

        # set variable attributes
        for key in par[var]:
            # Tag is used when user uses a different tag than the variable name. For now not used.
            tag = var
            # tag = tag_dict[var]
            if key not in ('column', 'dimensions'):
                attribute = key

                # let it default to default values
                if attribute != '_FillValue':
                    if attribute == 'missing_value':
                        setattr(all_vars[variable[var]], attribute, get_float(par[var][key]))
                    elif '_range' in attribute:
                        this_range = list()
                        this_value = par[var][key].replace('[', '').replace(']', '').strip()
                        if ',' in this_value:
                            this_value = this_value.split(',')
                        else:
                            this_value = this_value.split()

                        for item in this_value:
                            this_range.append(get_float(item))
                        setattr(all_vars[variable[var]], attribute, np.array(this_range, dtype=VAR_DTYPE))
                    else:
                        setattr(all_vars[variable[var]], attribute, par[var][key])

    header = header_tags(par)

    for line in data:
        y_val = line[header.index('y')]
        y_index = y_list.index(y_val)
        x_val = line[header.index('x')]
        x_index = x_list.index(x_val)
        for var in var_list:
            if int(par[var]['dimensions']) == 3:
                z_val = line[header.index('z')]
                z_index = z_list.index(z_val)
            value = line[header.index(var)]
            if int(par[var]['dimensions']) == 2:
                # For the 2D variables, we only keep the first set of values.
                if np.isnan(var_values[f"{variable[var]}_done"][y_index, x_index]):
                    var_values[variable[var]][y_index, x_index] = value
                    var_values[f"{variable[var]}_done"][y_index, x_index] = 1
            else:
                var_values[variable[var]][z_index, y_index, x_index] = value

    for var in var_list:
        all_vars[variable[var]][:] = var_values[variable[var]]

    if debug:
        print(f"\n\n[INFO] {par[var]['dimensions']}D variables:", flush=True)
        for var in var_list:
            print(f'\n\n{variable[var]}:{all_vars[variable[var]]}', flush=True)
    else:
        dot()

    return ds, all_vars, tag_dict


def set_global_attributes(ds, netcdf_file, par):
    """ set global attributes

    Keyword arguments:
    ds: the working netCDF dataset
    netcdf_file: name of the input GeoCSV file
    par" GeoCSV parameters

    Return values:
    ds: the input dataset with global attributes
    """

    if debug:
        print('\n\n[INFO] global attributes:\n', flush=True)
    else:
        dot()

    history = None
    for key in par['global']:
        attr_name = key
        if debug:
            print(f"{attr_name}: {par['global'][key]}\n")
        else:
            dot()
        setattr(ds, attr_name, get_float(par['global'][key]))
        if attr_name == 'history':
            history = par['global'][key]

    if history is not None:
        ds.history = f'{datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")} ' \
                     f'Converted to netCDF by ' \
                     f'{script} {version} from {os.path.basename(netcdf_file)}\n{history}'
    else:
        ds.history = f'{datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")} ' \
                     f'Converted to netCDF by ' \
                     f'{script} {version} from {os.path.basename(netcdf_file)}'

    return ds


def display_headers(model_data, data_params):
    """extract and display netCDF and GeoCSV header information

    Keyword arguments:
    model_data: Dataset instance of the model_file

    Return values:
    None, output of the header information
    """
    # GeoCSV header
    print('\n\nGeoCSV header information:\n\n', flush=True)
    for key in data_params.keys():
        if key == 'header':
            continue
        print(f'# {key}: {data_params[key]}', flush=True)

    # netCDF header
    print('\n\nnetCDF header information:\n\n', flush=True)

    # dimension information.
    nc_dimensions = [dim for dim in model_data.dimensions]  # list of netCDF dimensions
    print('\tdimensions:', flush=True)
    for dim in nc_dimensions:
        print(f'\t\t{model_data.dimensions[dim].name} {model_data.dimensions[dim].size}', flush=True)

    # variable information.
    nc_vars = [var for var in model_data.variables]  # list of nc variables

    print('\n\tvariables:', flush=True)
    for var in nc_vars:
        if var not in nc_dimensions:
            print(f'\t\t{var}:', flush=True)
            for attr, value in vars(model_data.variables[var]).items():
                print(f'\t\t\t{attr} = {value}', flush=True)

    # global attributes
    print('\n\tglobal attributes:', flush=True)
    for attr, value in vars(model_data).items():
        if isinstance(value, str):
            value = value.replace('\n', '; ')
        print(f'\t\t\t{attr} = {value}', flush=True)


def check_geocsv_file():
    """check the input GeoCSV model file and make sure it exist and extract info

    Return values:
    model_file_name: model file name (full path)
    file_name_tag:  file name with no extension
    """
    # check the model file and extract necessary information
    # must be in the argument list
    if geocsv_file_name is None:
        print('\n[ERR] the netCDF model file name is required', flush=True)
        usage()
        sys.exit(1)

    # user may provide full path
    elif os.path.isfile(geocsv_file_name):
        model_file_name = geocsv_file_name
        file_name_tag, ext = os.path.splitext(model_file_name)

    # user may place it under the data directory
    elif os.path.isfile(os.path.join(DATA_DIR, geocsv_file_name)):
        model_file_name = os.path.join(DATA_DIR, geocsv_file_name)
        file_name_tag, ext = os.path.splitext(model_file_name)

    # could not find the file
    else:
        print(f'\n[ERR] could not find the GeoCSV model file {geocsv_file_name}', flush=True)
        usage()
        sys.exit(1)

    return model_file_name, file_name_tag


# The main code.
try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hHi:d', ['help', 'input=', 'debug', 'header'])
    for opt, arg in options:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ('-i', '--input'):
            geocsv_file_name = arg
        elif opt in ('-d', '--debug'):
            debug = True
        elif opt in ('-H', '--header'):
            view_header = True
except getopt.GetoptError:
    usage()
    sys.exit(2)

# The model file name.
model_file, base_file_name = check_geocsv_file()
params, data_matrix = read_geocsv(model_file)

data_columns = list(zip(*data_matrix))
y_list, x_list, z_list = get_dimensions(data_columns, params)

# Create the netCDF file.
if view_header:
    dataset = Dataset('temp.nc', 'w', diskless=True, format=NETCDF_FORMAT)
else:
    dataset = Dataset(f'{base_file_name}.nc', 'w', format=NETCDF_FORMAT)

# Set the dimensions.
dataset = set_dimensions(dataset, params, x_list, y_list, z_list)

# Create the coordinate variables.
dataset, ys, xs, zs = create_coordinate_variables(dataset, params)

# Assign data to the coordinates.
ys[:] = y_list
xs[:] = x_list
if zs is not None:
    zs[:] = z_list
if debug:
    print(f"\n\n[INFO] coordinates:\n\n{params['y']['variable']}:{ys},"
          f"\n{params['x']['variable']}:{xs},\n{params['z']['variable']}:{zs}", flush=True)
else:
    dot()

# Set the global attributes.
dataset = set_global_attributes(dataset, model_file, params)

# Create variables.
dataset, variables, tag_dict = create_variables(dataset, params, data_matrix, y_list, x_list, z_list)

if view_header:
    display_headers(dataset, params)
    sys.exit(0)

if debug:
    print(f'\n\n[INFO] DATASET:\n{dataset}', flush=True)
else:
    dot()
    print(' done')

print(f'\n[INFO] Output netCDF file: {base_file_name}.nc\n', flush=True)

dataset.close()
