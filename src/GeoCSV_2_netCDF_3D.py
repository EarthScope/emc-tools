#!/usr/bin/env python
import sys
import os
from operator import itemgetter
from netCDF4 import Dataset
import numpy as np
import getopt
from datetime import datetime, timezone

'''
 NAME: GeoCSV_2_netCDF_3D.py - a Python script to read a 3D GeoCSV Earth model file and convert it to netCDF format.

 USAGE:
        
                            input GeoCSV Earth model file           
                            |       debug mode      
                            |       |   display headers only           
                            |       |   |            
       GeoCSV_2_netCDF_3D  -i FILE -d  -H

 HISTORY:
   2018-10-25 IRIS DMC Manoch: expanded the error message R.0.5.2018.298
   2018-10-22 IRIS DMC Manoch: initial release R.0.5.2018.295
   2018-10-17 IRIS DMC Manoch: fill_value and missing_value attributes are now supported
   2018-10-16 IRIS DMC Manoch: created VERSION 2018.289
'''

SCRIPT = os.path.basename(sys.argv[0])
VERSION = 'R.0.5.2018.298'
print('\n\n[INFO] {} version {}'.format(SCRIPT, VERSION), flush=True)

DEBUG = False

# root _directory
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# optional directory to store netCDF and GeoCSV files
DATA_DIR = os.path.join(ROOT_DIR, 'data')

NETCDF_FORMAT = 'NETCDF3_CLASSIC'
GEOCSV_FILE_NAME = None
BASE_NAME = None
VIEW_HEADER = False
VAR_DTYPE = 'f4'  # set variables as 'f4' (32-bit floating point) or 'f8' (64-bit floating point) variables


def usage():
    print('\n{} reads a 3D GetCSV Earth model file and convert it to netCDF format.\n\n'.format(SCRIPT))
    print(' USAGE:\n\n',
          '                         input getCSV Earth model file\n',
          '                         |          \n',
          '                         |          \n',
          '                         |          \n',
          '                         |       debug mode\n',
          '                         |       |    display headers only\n',
          '                         |       |    |\n',
          '{}   -i FILE -d   -H'.format(SCRIPT))


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
        this_value = float(this_string)
    except ValueError as e:
        this_value = this_string
    return this_value


def get_attribute(this_params, this_var, this_attribute):
    """return an attribute's value

    Keyword arguments:
    this_params: the working netCDF dictionary of the key-value pairs found in the header
    this_var: variable to check
    this_attribute: attribute to find

    Return values:
    this_value: value of the attribute, if found, otherwise None
    """
    this_value = None
    for this_key in this_params.keys():
        if this_var in this_key and this_attribute in this_key:
            this_value = this_params[this_key]

    return this_value


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
    print('\n[INFO] Input GeoCSV file: {}\n'.format(file_name), flush=True)
    if not DEBUG:
        dot()

    try:
       content = fp.read()
    except Exception as err:
        print('\n[Error] failed to read the input file!')
        print('{0}\n'.format(err))
        fp.close()
        sys.exit(2)

    lines = content.split('\n')
    fp.close()
    header_params = dict()
    data = list()
    found_header = False
    for i in range(len(lines)):
        if len(lines[i].strip()) <= 0:
            continue
        elif lines[i].strip()[0] == '#':
            this_line = lines[i].strip()[1:]
            if len(this_line.strip()) <= 0:
                continue
            if not DEBUG:
                dot()
            values = this_line.split(':')
            key = values[0].strip()
            header_params[key] = this_line.replace(key + ':', '').strip()
        else:
            # the first data line is the header
            this_line = lines[i].strip()
            if not found_header:
                header_params['header'] = this_line.split(header_params['delimiter'])
                found_header = True
                continue
            data.append(list(str2float(this_line.split(header_params['delimiter']))))
    if 'latitude_column' not in header_params.keys():
        header_params['latitude_column'] = 'latitude'
    if 'longitude_column' not in header_params.keys():
        header_params['longitude_column'] = 'longitude'
    if 'depth_column' not in header_params.keys() and 'elevation_column' not in header_params.keys():
        header_params['level_column'] = 'depth'
        header_params['depth_column'] = 'depth'
    elif 'depth_column' in header_params.keys():
        header_params['level_column'] = header_params['depth_column']
    elif 'elevation_column' in header_params.keys():
        header_params['level_column'] = header_params['elevation_column']

    if not DEBUG:
        dot()

    header = header_params['header']
    latitude_column = header.index(header_params['latitude_column'])
    longitude_column = header.index(header_params['longitude_column'])
    level_column = header.index(header_params['depth_column'])

    if not VIEW_HEADER:
        data.sort(key=itemgetter(level_column, latitude_column, longitude_column))

    return header_params, data


def get_dimensions(column_data, header_params):
    """extract sorted unique values of the dimension columns

    Keyword arguments:
    column_data: matrix representing columns of a data matrix
    header_params: GeoCSV parameters

    Return values:
    lat: sorted unique latitude values
    lon: sorted unique longitude values
    level: sorted unique level values
    """

    header = header_params['header']
    latitude_column = header.index(header_params['latitude_column'])
    longitude_column = header.index(header_params['longitude_column'])
    level_column = header.index(header_params['depth_column'])

    if not DEBUG:
        dot()
    lat = list(sorted(set(column_data[latitude_column])))
    if not DEBUG:
        dot()
    lon = list(sorted(set(column_data[longitude_column])))
    if not DEBUG:
        dot()
    level = list(sorted(set(column_data[level_column])))

    if DEBUG:
        print('\n\n[INFO] values:\n\n{}:{},\n\n{}:{},\n\n{}:{}'.format(header[latitude_column], lat,
                                                                       header[longitude_column], lon,
                                                                       header[level_column], level), flush=True)
    else:
        dot()
    return lat, lon, level


def set_dimensions(this_dataset, header_params, latitude_list, longitude_list, levels_list):
    """define variables' dimensions

    Keyword arguments:
    this_dataset: the working netCDF dataset
    header_params: GeoCSV parameters
    latitude_list: list of unique latitudes from GeoCSV file
    longitude_list: list of unique longitudes from GeoCSV file
    levels_list: list of unique levels from GeoCSV file

    Return values:
    this_dataset: the input dataset with dimensions set
    """
    lat_dim = this_dataset.createDimension(header_params['latitude_column'], len(latitude_list))
    lon_dim = this_dataset.createDimension(header_params['longitude_column'], len(longitude_list))
    level_dim = this_dataset.createDimension(header_params['level_column'], len(levels_list))

    if DEBUG:
        print('\n\n[INFO] dimensions:\n\n{}:{},\n{}:{},\n{}:{}'.format(header_params['latitude_column'], lat_dim,
                                                                       header_params['longitude_column'], lon_dim,
                                                                       header_params['level_column'], level_dim),
              flush=True)
    else:
        dot()

    return this_dataset


def create_coordinate_variables(this_dataset, this_params):
    """ create coordinate variables

    Keyword arguments:
    this_dataset: the working netCDF dataset
    this_params: GeoCSV parameters

    Return values:
    dataset: the input dataset with coordinate variables set
    latitudes: latitude variable
    longitudes: longitude variable
    levels: level variable
    """

    coordinate = dict()
    this_value = get_attribute(this_params, this_params['latitude_column'], '_FillValue')
    if this_value:
        coordinate[params['latitude_column']] = this_dataset.createVariable(
            this_params['latitude_column'], VAR_DTYPE, (this_params['latitude_column'],), fill_value=this_value)
    else:
        coordinate[params['latitude_column']] = this_dataset.createVariable(
            this_params['latitude_column'], VAR_DTYPE, (this_params['latitude_column'],))

    this_value = get_attribute(this_params, this_params['longitude_column'], '_FillValue')
    if this_value:
        coordinate[this_params['longitude_column']] = this_dataset.createVariable(
            this_params['longitude_column'], VAR_DTYPE, (this_params['longitude_column'],), fill_value=this_value)
    else:
        coordinate[this_params['longitude_column']] = this_dataset.createVariable(
            this_params['longitude_column'], VAR_DTYPE, (this_params['longitude_column'],))

    this_value = get_attribute(this_params, this_params['level_column'], '_FillValue')
    if this_value:
        coordinate[this_params['level_column']] = this_dataset.createVariable(
            params['level_column'], VAR_DTYPE, (this_params['level_column'],), fill_value=this_value)
    else:
        coordinate[this_params['level_column']] = this_dataset.createVariable(
            params['level_column'], VAR_DTYPE, (this_params['level_column'],))

    for var in (this_params['latitude_column'], this_params['longitude_column'], this_params['level_column']):
        for key in this_params.keys():
            if key.strip().startswith('{}_'.format(var)):
                attribute = key.replace(var, '').strip()

                # let it default to default values
                if attribute != '_FillValue':
                    if attribute == 'missing_value':
                        setattr(coordinate[var], attribute, get_float(this_params[key]))
                    elif '_range' in attribute:
                        this_range = list()
                        this_value = this_params[key].replace('[','').replace(']','').strip().split(',')
                        for item in this_value:
                            this_range.append(get_float(item))
                        setattr(coordinate[var], attribute, this_range)
                    else:
                        setattr(coordinate[var], attribute, this_params[key])

    if DEBUG:
        print('\n\n[INFO] coordinate variables:\n\n{}:{},\n{}:{},\n{}:{}'.format(
            this_params['latitude_column'], coordinate[this_params['latitude_column']], this_params['longitude_column'],
            coordinate[this_params['longitude_column']],
            this_params['level_column'],
            coordinate[this_params['level_column']]), flush=True)
    else:
        dot()

    return this_dataset, coordinate[this_params['latitude_column']], coordinate[this_params['longitude_column']],\
        coordinate[this_params['level_column']]


def create_3d_variables(this_dataset, this_params, data, latitude_list, longitude_list, levels_list):
    """ create 3D variables

    Keyword arguments:
    this_dataset: the working netCDF dataset
    this_params" GeoCSV parameters
    data: block of data
    latitude_list: list of unique latitudes from GeoCSV file
    longitude_list: list of unique longitudes from GeoCSV file
    levels_list: list of unique levels from GeoCSV file

    Return values:
    this_dataset: the input dataset with coordinate variables set
    all_vars: a dictionary of 3D variables
    """

    header = this_params['header']

    var_list = list(set(header) - set([this_params['latitude_column'], this_params['longitude_column'],
                                       this_params['level_column']]))
    all_vars = dict()
    var_values = dict()

    for var in var_list:
        this_value = get_attribute(this_params, var, '_FillValue')
        if this_value:
            all_vars[var] = this_dataset.createVariable(
                var, VAR_DTYPE, (this_params['level_column'], this_params['latitude_column'],
                                  this_params['longitude_column']), fill_value=this_value)
        else:
            all_vars[var] = this_dataset.createVariable(
                var, VAR_DTYPE, (this_params['level_column'], this_params['latitude_column'],
                                  this_params['longitude_column']))

        var_values[var] = np.zeros((len(levels_list), len(latitude_list), len(longitude_list)))

        # set variable attributes
        for key in this_params.keys():
            if '{}_'.format(var) in key and key != '{}_column'.format(var):
                attribute = key.replace('{}_'.format(var), '').strip()

                # let it default to default values
                if attribute != '_FillValue':
                    if attribute == 'missing_value':
                        setattr(all_vars[var], attribute, get_float(this_params[key]))
                    elif '_range' in attribute:
                        this_range = list()
                        this_value = this_params[key].replace('[', '').replace(']', '').strip()
                        if ',' in this_value:
                            this_value = this_value.split(',')
                        else:
                            this_value = this_value.split()

                        for item in this_value:
                            this_range.append(get_float(item))
                        setattr(all_vars[var], attribute, this_range)
                    else:
                        setattr(all_vars[var], attribute, this_params[key])

    for line in data:
        lat = line[header.index(this_params['latitude_column'])]
        lat_index = latitude_list.index(lat)
        lon = line[header.index(this_params['longitude_column'])]
        lon_index = longitude_list.index(lon)
        lev = line[header.index(this_params['level_column'])]
        lev_index = level_list.index(lev)
        for var in var_list:
            value = line[header.index(var)]
            var_values[var][lev_index, lat_index, lon_index] = value

    for var in var_list:
        all_vars[var][:] = var_values[var]

    if DEBUG:
        print('\n\n[INFO] 3D variables:', flush=True)
        for var in var_list:
            print('\n\n{}:{}'.format(var, all_vars[var]), flush=True)
    else:
        dot()

    return this_dataset, all_vars


def set_global_attributes(this_dataset, this_file, header_params):
    """ set global attributes

    Keyword arguments:
    this_dataset: the working netCDF dataset
    this_file: name of the input GeoCSV file
    header_params" GeoCSV parameters

    Return values:
    this_dataset: the input dataset with global attributes
    """

    if DEBUG:
        print('\n\n[INFO] global attributes:\n', flush=True)
    else:
        dot()

    for key in header_params.keys():
        if key.strip().startswith('global_'):
            attr_name = key.replace('global_', '').strip()
            if DEBUG:
                print('{}: {}\n'.format(attr_name, header_params[key]))
            else:
                dot()

            setattr(this_dataset, attr_name, header_params[key])

    this_dataset.source = 'Converted from {}'.format(os.path.basename(this_file))
    this_dataset.history = 'Created by {} ({})'.format(SCRIPT, datetime.now(timezone.utc).isoformat(timespec='seconds'))

    return this_dataset


def display_headers(model_data, params):
    """extract and display netCDF and GeoCSV header information

    Keyword arguments:
    model_data: Dataset instance of the model_file

    Return values:
    None, output of the header information
    """
    # GeoCSV header
    print('\n\nGeoCSV header information:\n\n', flush=True)
    for key in params.keys():
        if key == 'header':
            continue
        print('# {}: {}'.format(key,params[key]), flush=True)

    # netCDF header
    print('\n\nnetCDF header information:\n\n', flush=True)

    # dimension information.
    nc_dims = [dim for dim in model_data.dimensions]  # list of netCDF dimensions
    print ('\tdimensions:', flush=True)
    for dim in nc_dims:
        print('\t\t{} {}'.format(model_data.dimensions[dim].name, model_data.dimensions[dim].size), flush=True)

    # variable information.
    nc_vars = [var for var in model_data.variables]  # list of nc variables

    print('\n\tvariables:', flush=True)
    for var in nc_vars:
        if var not in nc_dims:
            print('\t\t{}:'.format(var), flush=True)
            for attr, value in vars(model_data.variables[var]).items():
                print('\t\t\t{} = {}'.format(attr, value), flush=True)

    # global attributes
    print('\n\tglobal attributes:', flush=True)
    for attr, value in vars(model_data).items():
        value = value.replace('\n', ' ')
        print('\t\t\t{} = {}'.format(attr, value), flush=True)


def check_geocsv_file():
    """check the input GeoCSV model file and make sure it exist and extract info

    Return values:
    model_file_name: model file name (full path)
    file_name_tag:  file name with no extension
    """
    # check the model file and extract necessary information
    # must be in the argument list
    if GEOCSV_FILE_NAME is None:
        print('[ERROR] the netCDF model file name is required', flush=True)
        usage()
        sys.exit(1)

    # user may provide full path
    elif os.path.isfile(GEOCSV_FILE_NAME):
        model_file_name = GEOCSV_FILE_NAME
        file_name_tag, ext = os.path.splitext(model_file_name)

    # user may place it under the data directory
    elif os.path.isfile(os.path.join(DATA_DIR, GEOCSV_FILE_NAME)):
        model_file_name = os.path.join(DATA_DIR, GEOCSV_FILE_NAME)
        file_name_tag, ext = os.path.splitext(model_file_name)

    # could not find the file
    else:
        print('[ERROR] could not find the GeoCSV model file {}'.format(GEOCSV_FILE_NAME), flush=True)
        usage()
        sys.exit(1)

    return model_file_name, file_name_tag


try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hHi:d', ['help', 'input=', 'debug', 'header'])
    for opt, arg in options:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ('-i', '--input'):
            GEOCSV_FILE_NAME = arg
        elif opt in ('-d', '--debug'):
            DEBUG = True
        elif opt in ('-H', '--header'):
            VIEW_HEADER = True

except getopt.GetoptError:
    usage()
    sys.exit(2)

model_file, base_file_name = check_geocsv_file()
params, data_matrix = read_geocsv(model_file)

data_columns = list(zip(*data_matrix))

lat_list, lon_list, level_list = get_dimensions(data_columns, params)

# create the netCDF file
if VIEW_HEADER:
    dataset = Dataset('temp.nc', 'w', diskless=True, format=NETCDF_FORMAT)
else:
    dataset = Dataset('{}.nc'.format(base_file_name), 'w', format=NETCDF_FORMAT)

# set the dimensions
dataset = set_dimensions(dataset, params, lat_list, lon_list, level_list)

# create the coordinate variables
dataset, latitudes, longitudes, levels = create_coordinate_variables(dataset, params)

# assign data to coordinates
latitudes[:] = lat_list
longitudes[:] = lon_list
levels[:] = level_list
if DEBUG:
    print('\n\n[INFO] coordinates:\n\n{}:{},\n{}:{},\n{}:{}'.format(
        params['latitude_column'], latitudes, params['longitude_column'], longitudes, params['level_column'],
        levels), flush=True)
else:
    dot()

# set global attributes
set_global_attributes(dataset, model_file, params)

# create the 3D variables
dataset, variables = create_3d_variables(dataset, params, data_matrix, lat_list, lon_list, level_list)

if VIEW_HEADER:
    display_headers(dataset, params)
    sys.exit(0)

if DEBUG:
    print('\n\n[INFO] DATASET:\n{}'.format(dataset), flush=True)
else:
    dot()
    print(' done')

print('\n[INFO] Output netCDF file: {}.nc\n'.format(base_file_name), flush=True)

dataset.close()

