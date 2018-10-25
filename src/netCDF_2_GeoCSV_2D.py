#!/usr/bin/env python
import sys
import os
import getopt
from datetime import datetime
from netCDF4 import Dataset

'''
 NAME: netCDF_2_GeoCSV_2D.py - a Python script to read a 2D netCDF Earth model file and convert it to GeoCSV format.
               
 USAGE:
                         
                           input netCDF Earth model file
                           |    
                           |       longitude variable (default longitude)                    
                           |       | 
                           |       |       latitude variable (default latitude)
                           |       |       |   
                           |       |       |      debug mode
                           |       |       |      |
                           |       |       |      |  display headers only
                           |       |       |      |  |  
    netCDF_2_GeoCSV_3D.py -i FILE -x long -y lat -d -H

 HISTORY:
   2018-10-25 IRIS DMC Manoch: created R.0.5.2018.298
'''

SCRIPT = os.path.basename(sys.argv[0])
VERSION = 'R.0.5.2018.298'
print('\n\n[INFO] {} version {}'.format(SCRIPT, VERSION), flush=True)

DEBUG = False

GEOCSV_VERSION = 'GeoCSV2.0'

# root _directory
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# optional directory to store netCDF and GeoCSV files
DATA_DIR = os.path.join(ROOT_DIR, 'data')

LON_VARIABLE = 'longitude'  # must match the netCDF file's longitude variable
LAT_VARIABLE = 'latitude'  # must match the netCDF file's longitude variable
VALID_MODES = {'depth': 'km', 'single': ''}
DELIMITER = '|'
NETCDF_FILE_NAME = None
BASE_NAME = None
VIEW_HEADER = False

# output mode (depth | lat | lon) as individual files based on depth, lat, lon
# or (single) as a single file
OUTPUT_MODE = 'single'

def dot():
    """print a dot on screen"""

    sys.stdout.write('.')
    sys.stdout.flush()


def get_variable_attributes(model_data, header, variable, variable_name):
    """add  variable attributes to the header

    Keyword arguments:
    model_data: Dataset instance of the model_file
    header: GeoCSV header variable for the model
    variable : variable to add
    variable_name: variable name used to represent this variable

    Return values:
    the geoCSV header for the model
    """
    header.append('# {}_column: {}\n'.format(variable, variable_name))
    for attr, value in vars(model_data.variables[variable]).items():
        if '_range' in attr:
            header.append('# {}_{}: {},{}\n'.format(variable, attr, value[0], value[1]))
        else:
            header.append('# {}_{}: {}\n'.format(variable, attr, value))
    return header


def get_model_header(model_file, model_data):
    """create GeoCSV header for the model

    Keyword arguments:
    model_file: the netCDF model file name
    model_data: Dataset instance of the model_file

    Return values:
    the geoCSV header for the model
    """
    header = list()
    # GeoCSV header
    header.append('# dataset: {}\n'.format(GEOCSV_VERSION))
    header.append('# created: {} UTC ({})\n'.format(datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), SCRIPT))
    header.append('# netCDF_file: {}\n'.format(os.path.basename(model_file)))
    header.append('# delimiter: {}\n'.format(DELIMITER))

    # global attributes
    for attr, value in vars(model_data).items():
        if isinstance(value, str):
            value = value.replace('\n', ' ')
        header.append('# global_{}: {}\n'.format(attr, value))

    # variable s
    header = get_variable_attributes(model_data, header, 'latitude', LAT_VARIABLE)
    header = get_variable_attributes(model_data, header, 'longitude', LON_VARIABLE)

    return ''.join(header)


def get_var_header(model_data, var):
    """create GeoCSV header for a variable

    Keyword arguments:
    model_data: Dataset instance of the model_file
    var: the netCDF model file variable

    Return values:
    the geoCSV header for the variable
    """
    header = list()
    header = get_variable_attributes(model_data, header, var, var)
    return ''.join(header)


def usage_csv():
    print('\n{} reads a 3D netCDF Earth model file and convert it to GeoCSV format.\n\n'.format(SCRIPT))
    print(' USAGE:\n\n',
          '              \n',
          '                        input netCDF Earth model file\n',
          '                        |\n',
          '                        |       longitude variable (default longitude)\n',
          '                        |       |\n',
          '                        |       |       latitude variable (default latitude)\n',
          '                        |       |       |\n',
          '                        |       |       |      debug mode\n',
          '                        |       |       |      |        \n',
          '                        |       |       |      |        display headers only\n',
          '                        |       |       |      |        |       \n',
          '{}   -i FILE -x long -y lat -d   -H'.format(SCRIPT))


def check_netcdf_file():
    """check the input netCDF model file and make sure it exist and extract info

    Return values:
    this_file:  file name with no extension
    """
    # check the model file and extract necessary information
    # must be in the argument list
    if NETCDF_FILE_NAME is None:
        print('[ERROR] the netCDF model file name is required', flush=True)
        usage_csv()
        sys.exit(1)

    # user may provide full path
    elif os.path.isfile(NETCDF_FILE_NAME):
        model_file_name = NETCDF_FILE_NAME
        base_file_name, ext = os.path.splitext(model_file_name)

    # user may place it under the data directory
    elif os.path.isfile(os.path.join(DATA_DIR, NETCDF_FILE_NAME)):
        model_file_name = os.path.join(DATA_DIR, NETCDF_FILE_NAME)
        base_file_name, ext = os.path.splitext(model_file_name)

    # could not find the file
    else:
        print('[ERROR] could not find the netCDF model file {}'.format(NETCDF_FILE_NAME), flush=True)
        usage_csv()
        sys.exit(1)

    return model_file_name, base_file_name


def display_headers(model_file, model_data):
    """extract and display netCDF and GeoCSV header information

    Keyword arguments:
    model_data: Dataset instance of the model_file

    Return values:
    None, output of the header information
    """
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
        if isinstance(value, str):
            value = value.replace('\n', ' ')
        print('\t\t\t{} = {}'.format(attr, value), flush=True)

    # GeoCSV header
    print('\n\nGeoCSV header information:\n\n{}\n\n'.format(get_model_header(model_file, model_data)), flush=True)


def make_model_geocsv():
    """create GeoCSV file from a netCDF model file

    Keyword arguments:
    model_file: the netCDF model file name
    """

    model_file, base_file_name = check_netcdf_file()
    print('[INFO] Input netCDF File: {}'.format(model_file), flush=True)

    data_header = list()
    model_data = Dataset(model_file)
    try:
        # conversion to string is done to preserve precision
        lat = list()
        lon = list()
        depth = list()
        for this_value in model_data.variables[LAT_VARIABLE][:]:
            lat.append("{}".format(str(this_value)))
        for this_value in model_data.variables[LON_VARIABLE][:]:
            lon.append("{}".format(str(this_value)))
    except Exception:
        print('\n[Error] the expected variables ({}, {}) not in the variable list: {}\n'.format(
            LAT_VARIABLE, LON_VARIABLE, str(list(model_data.variables.keys()))))
        sys.exit(1)

    emcin = {}

    # make sure this is a 2D netCDF file
    var_2d = list()
    for var in model_data.variables.keys():
        if len(model_data.variables[var].shape) == 2:
            var_2d.append(var)
    if len(var_2d) <= 0:
        print('\n[ERROR] not a 2D netCDF file\n\n', flush=True)
        sys.exit(1)

    # the standard order is (Z, Y, X) or (depth, latitude, longitude)
    if DEBUG:
        print('[INFO] Mode: {}'.format(OUTPUT_MODE), flush=True)
        print('[INFO] 2D Variables: {}'.format(var_2d), flush=True)

    if VIEW_HEADER:
        display_headers(model_file, model_data)
        sys.exit(0)

    output_data = list()
    # placeholder for future development
    if True:
        data_header = list()
        output_file = '{}.csv'.format(base_file_name)
        fp = open(output_file, 'w')
        print('[INFO] Output file: {}'.format(output_file), flush=True)
        fp.write(get_model_header(model_file, model_data))
        data_header.append('{}{}{}'.format(LAT_VARIABLE, DELIMITER, LON_VARIABLE))

        dot()

        index = [-1, -1]
        for i, this_lat in enumerate(lat):
            for j, this_lon in enumerate(lon):
                output_data.append('{}{}{}'.format(str(this_lat), DELIMITER, str(this_lon)))

                for var in model_data.variables.keys():
                    lat_index = None
                    lon_index = None
                    if var.encode('ascii', 'ignore').decode("utf-8") not in [LAT_VARIABLE, LON_VARIABLE,]:
                        if not i and not j:
                            fp.write(get_var_header(model_data, var))
                            data_header.append('{}{}'.format(DELIMITER, var))
                        # find the variable ordering
                        if lat_index is None:
                            for l in range(len(model_data.variables[var].dimensions)):
                                if model_data.variables[var].dimensions[l].encode('ascii', 'ignore').decode(
                                        "utf-8") == LON_VARIABLE:
                                    lon_index = l
                                else:
                                    lat_index = l

                            if var not in emcin.keys():
                                try:
                                    emcin[var] = model_data.variables[var][:]
                                except Exception as err:
                                    print('\n[Error] problem reading variable "{}"'.format(var))
                                    print('{0}\n'.format(err))
                                    sys.exit(2)

                        index[lat_index] = i
                        index[lon_index] = j
                        # nan values, otherwise we write string to preserve the precision
                        if str(emcin[var][index[0]][index[1]]) == '--':
                            output_data.append('{}{}'.format(DELIMITER,
                                                             float(emcin[var][index[0]][index[1]])))
                        else:
                            # conversion to string is done to preserve precision
                            output_data.append('{}{}'.format(DELIMITER,
                                                             str(emcin[var][index[0]][index[1]])))
                output_data.append('\n')

    fp.write('{}\n'.format(''.join(data_header)))
    fp.write(''.join(output_data))
    fp.close()


try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hdHi:v:x:y:z:m:', ['help', 'input=', 'latitude=', 'longitude=',
                                                                        'debug', 'header'])
    for opt, arg in options:
        if opt == '-h':
            usage_csv()
            sys.exit()
        elif opt in ('-d', '--debug'):
            DEBUG = True
        elif opt in ('-H', '--header'):
            VIEW_HEADER = True
        elif opt in ('-i', '--input'):
            NETCDF_FILE_NAME = arg
        elif opt in ('-x', '--longitude'):
            LON_VARIABLE = arg
        elif opt in ('-y', '--latitude'):
            LAT_VARIABLE = arg
except getopt.GetoptError:
    usage_csv()
    sys.exit(2)


make_model_geocsv()
print(' done')
