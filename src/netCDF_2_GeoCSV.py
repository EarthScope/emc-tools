#!/usr/bin/env python
import sys
import os
import getopt
import math
from datetime import datetime, timezone, UTC
from netCDF4 import Dataset

"""
 NAME: netCDF_2_GeoCSV.py - a Python script to read a 3D netCDF Earth model file and convert it to GeoCSV format.
               
 USAGE:
                         
                           input netCDF Earth model file
                           |    
                           |       longitude variable (default longitude)                    
                           |       | 
                           |       |       latitude variable (default latitude)
                           |       |       |   
                           |       |       |      depth variable, default depth
                           |       |       |      |         
                           |       |       |      |        output mode (single, depth)
                           |       |       |      |        |        
                           |       |       |      |        |        debug mode
                           |       |       |      |        |        |
                           |       |       |      |        |        |  display headers only
                           |       |       |      |        |        |  | 
    netCDF_2_GeoCSV.py -i FILE -x long -y lat -z depth -m depth -d -H

 HISTORY:
   2022-08-25 IRIS DMC Manoch: v2022.237 R2 release combines 2D and 3D scripts and also adds support for the 
                               models with projected coordinate systems.
   2020-09-29 IRIS DMC Manoch: v2020.273 fix for a bug when x,y,or z variables were changed.
   2020-06-16 IRIS DMC Manoch: v2020.168 minor style updates
   2020-01-06 IRIS DMC Manoch: v2020.006 R1 added history if it does not exist. 
                               The history now includes the source file name
   2018-10-25 IRIS DMC Manoch: expanded the error message R.0.5.2018.298
   2018-10-22 IRIS DMC Manoch: initial release R.0.5.2018.295
   2018-09-15 IRIS DMC Manoch: created VERSION 2018.258 
   2018-10-16 IRIS DMC Manoch: updated to include full netCDF header VERSION 2018.289
"""

SCRIPT = os.path.basename(sys.argv[0])
VERSION = "v2022.237"
print(f"\n\n[INFO] {SCRIPT} version {VERSION}", flush=True)

DEBUG = False

GEOCSV_VERSION = "GeoCSV2.0"

# root _directory
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# optional directory to store netCDF and GeoCSV files
DATA_DIR = os.path.join(ROOT_DIR, "data")

X_VARIABLE = None  # must match the netCDF file's x variable
Y_VARIABLE = None  # must match the netCDF file's y variable
Z_VARIABLE = None  # must match the netCDF file's z variable
ndim = 3
ndim_label = "3D"
VALID_MODES = {"depth": "km", "single": ""}
DELIMITER = "|"
NETCDF_FILE_NAME = None
BASE_NAME = None
VIEW_HEADER = False

# output mode (depth | lat | lon) as individual files based on depth, lat, lon
# or (single) as a single file
OUTPUT_MODE = "single"


def dot():
    """print a dot on screen"""

    sys.stdout.write(".")
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
    header.append(f'# {variable_name.replace("_", "-")}_column: {variable}\n')
    header.append(f'# {variable_name.replace("_", "-")}_variable: {variable}\n')
    header.append(
        f'# {variable_name.replace("_", "-")}_dimensions: {len(model_data.variables[variable].shape)}\n'
    )

    for attr, value in vars(model_data.variables[variable]).items():
        if "_range" in attr:
            header.append(
                f'# {variable_name.replace("_", "-")}_{attr}: {value[0]},{value[1]}\n'
            )
        else:
            header.append(f'# {variable_name.replace("_", "-")}_{attr}: {value}\n')
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
    header.append(f"# dataset: {GEOCSV_VERSION}\n")
    header.append(
        f'# created: {datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S")} UTC ({SCRIPT})\n'
    )
    header.append(f"# netCDF_file: {os.path.basename(model_file)}\n")
    header.append(f"# delimiter: {DELIMITER}\n")

    # global attributes
    history_done = False
    history = (
        f'{datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")} Converted to GeoCSV by {SCRIPT} ,'
        f"{VERSION} "
        f"from {NETCDF_FILE_NAME}"
    )
    for attr, value in vars(model_data).items():
        if isinstance(value, str):
            value = value.replace("\n", "; ")
        if attr.lower() == "history":
            value = f"{history}; {value}"
            history_done = True
        header.append(f"# global_{attr}: {value}\n")

    if not history_done:
        header.append(f"# global_history: {history}\n")

    # variable s
    if VIEW_HEADER and (X_VARIABLE is None or Y_VARIABLE is None):
        header.append(
            "\n\nNOTE: For information on variables, the -x, -y, and -z parameters are required\n\n"
        )
    else:
        header = get_variable_attributes(
            model_data, header, Y_VARIABLE, "y"
        )  # , Y_VARIABLE)
        header = get_variable_attributes(
            model_data, header, X_VARIABLE, "x"
        )  # , X_VARIABLE)
        if ndim == 3:
            header = get_variable_attributes(
                model_data, header, Z_VARIABLE, "z"
            )  # , Z_VARIABLE)

    return "".join(header)


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
    return "".join(header)


def usage():
    print(
        f"\n{SCRIPT} reads a 2 or 3 dimensional netCDF Earth model file and convert it to GeoCSV format.\n\n"
    )
    print(
        f" USAGE:\n\n",
        f"{SCRIPT} -i FILE -x longitude -y latitude -z depth -m mode -d -H\n"
        f"\t-i [required] FILE is the input netCDF Earth model file\n",
        f"\t-x [required] the netCDF variable representing the x-axis (must match the netCDF file's x variable)\n",
        f"\t-y [required] the netCDF variable representing the y-axis (must match the netCDF file's y variable)\n",
        f"\t-z [required] the netCDF variable representing the z-axis (must match the netCDF file's z variable"
        f" OR set to 'none' for 2D models)\n",
        f"\t-m [default:{OUTPUT_MODE}] output mode [single: single GeoCSV file OR depth: one GeoCSV file per depth\n",
        f"\t-d [default:{DEBUG}] debug mode\n",
        f"\t-H [default:{VIEW_HEADER}] display headers only\n",
    )


def check_netcdf_file():
    """check the input netCDF model file and make sure it exist and extract info

    Return values:
    this_file:  file name with no extension
    """
    # check the model file and extract necessary information
    # must be in the argument list
    if NETCDF_FILE_NAME is None:
        print("[ERR] the netCDF model file name is required", flush=True)
        usage()
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
        print(
            f"[ERR] could not find the netCDF model file {NETCDF_FILE_NAME}", flush=True
        )
        usage()
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
    print("\n\nnetCDF header information:\n\n", flush=True)

    # dimension information.
    nc_dims = [dim for dim in model_data.dimensions]  # list of netCDF dimensions
    print("\tdimensions:", flush=True)
    for dim in nc_dims:
        print(
            f"\t\t{model_data.dimensions[dim].name} {model_data.dimensions[dim].size}",
            flush=True,
        )

    # variable information.
    nc_vars = [var for var in model_data.variables]  # list of nc variables

    print("\n\tvariables:", flush=True)
    for var in nc_vars:
        if var not in nc_dims:
            print(f"\t\t{var}:", flush=True)
            for attr, value in vars(model_data.variables[var]).items():
                print(f"\t\t\t{attr} = {value}", flush=True)

    # global attributes
    print('"\n\tglobal attributes:', flush=True)
    for attr, value in vars(model_data).items():
        if isinstance(value, str):
            value = value.replace("\n", "; ")
        print(f"\t\t\t{attr} = {value}", flush=True)

    # GeoCSV header
    print(
        f"\n\nGeoCSV header information:\n\n{get_model_header(model_file, model_data)}\n\n",
        flush=True,
    )


def make_model_geocsv():
    """create GeoCSV file from a netCDF model file

    Keyword arguments:
    model_file: the netCDF model file name
    """

    model_file, base_file_name = check_netcdf_file()
    print(f"[INFO] Input netCDF File: {model_file}", flush=True)

    data_header = list()
    model_data = Dataset(model_file)

    if VIEW_HEADER:
        display_headers(model_file, model_data)
        sys.exit(0)

    try:
        # Conversion to string is done to preserve precision
        y = list()
        x = list()
        z = list()
        for this_value in model_data.variables[Y_VARIABLE][:]:
            y.append(f"{str(this_value)}")
        for this_value in model_data.variables[X_VARIABLE][:]:
            x.append(f"{str(this_value)}")

        if Z_VARIABLE is not None:
            for this_value in model_data.variables[Z_VARIABLE][:]:
                z.append(f"{str(this_value)}")
        else:
            z.append(None)
    except Exception as ex:
        if ndim == 2:
            print(
                f"\n[ERR] The expected variables ({Y_VARIABLE} and {X_VARIABLE}) "
                f"are not in the variable list: {list(model_data.variables.keys())}\n{ex}"
            )
        else:
            print(
                f"\n[ERR] The expected variables ({Y_VARIABLE}, {X_VARIABLE}, {Z_VARIABLE}) "
                f"are not in the variable list: {list(model_data.variables.keys())}\n{ex}"
            )
        sys.exit(1)

    emcin = {}

    # This is a dictionary holding the variable dimensions.
    var_ndim = dict()
    dim_var_count = 0
    for var in model_data.variables.keys():
        # Check the model dimensions.
        if var in (X_VARIABLE, Y_VARIABLE, Z_VARIABLE):
            dim_var_count += 1
            var_ndim[var] = len(model_data.variables[var].shape)
        else:
            # Save variable dimensions.
            var_ndim[var] = len(model_data.variables[var].shape)
    # Last check on the model dimension.
    if dim_var_count != ndim:
        print(f"\n[ERR] not a {ndim_label} netCDF file\n\n", flush=True)
        sys.exit(1)

    # The standard order is (Z, Y, X) or (depth, latitude, longitude)
    if DEBUG:
        print(f"[INFO] Mode: {OUTPUT_MODE}", flush=True)
        print(f"[INFO] {ndim_label} Variables: {var}", flush=True)

    output_data = list()
    depth_index = dict()
    lat_index = dict()
    lon_index = dict()

    # Go through each depth.
    index = [-1, -1, -1]
    for k, this_depth in enumerate(z):
        # Single output file option, if requested or if the model is 2D.
        if (OUTPUT_MODE == "single" and k == 0) or ndim == 2:
            data_header = list()
            output_file = f"{base_file_name}.csv"
            fp = open(output_file, "w")
            print(f"[INFO] Output file: {output_file}", flush=True)
            fp.write(get_model_header(model_file, model_data))
            if ndim == 2:
                data_header.append(f"{Y_VARIABLE}{DELIMITER}{X_VARIABLE}")
            else:
                data_header.append(
                    f"{Y_VARIABLE}{DELIMITER}{X_VARIABLE}{DELIMITER}{Z_VARIABLE}"
                )
        # Depth-based output files.
        elif OUTPUT_MODE == "depth":
            output_data = list()
            data_header = list()
            output_file = os.path.join(
                f"{base_file_name}_{this_depth}_{VALID_MODES[OUTPUT_MODE]}.csv"
            )
            fp = open(output_file, "w")
            print(f"[INFO] Output file: {output_file}", flush=True)
            fp.write(get_model_header(model_file, model_data))
            data_header.append(f"# depth: {this_depth}\n")
            data_header.append(f"{Y_VARIABLE}{DELIMITER}{X_VARIABLE}")

        if DEBUG and this_depth is not None:
            print(f"[INFO] Processing depth: {this_depth}", flush=True)
        elif this_depth is not None:
            # Show the progress.
            if k == len(z) - 1:
                print(f"{this_depth}", flush=True)
            elif k == 0:
                zero_depth = this_depth
            elif k == 1:
                if ndim > 2:
                    print(f"[INFO] Depth range: {z[0]} to {z[-1]}", flush=True)
                print(f"{zero_depth}, {this_depth},", end=" ", flush=True)
            else:
                print(f"{this_depth},", end=" ", flush=True)

        # Go through each latitude and convert to string to preserve precision
        for i, this_lat in enumerate(y):
            _lat = str(this_lat)
            # Go through each longitude.
            for j, this_lon in enumerate(x):
                _lon = str(this_lon)
                # Data header line.
                if OUTPUT_MODE == "single" and ndim == 3:
                    _depth = str(this_depth)
                    output_data.append(
                        f"{_lat:s}{DELIMITER}{_lon:s}{DELIMITER}{_depth:s}"
                    )
                else:
                    output_data.append(f"{_lat:s}{DELIMITER}{_lon:s}")

                # Go through each model variable.
                for var_index, var_value in enumerate(list(model_data.variables)):
                    var = var_value.encode("ascii", "ignore").decode("utf-8")

                    # Model variables only.
                    if var.encode("ascii", "ignore").decode("utf-8") not in [
                        X_VARIABLE,
                        Y_VARIABLE,
                        Z_VARIABLE,
                    ]:
                        # Do this for the first point, when all indices are zero.
                        if (
                            OUTPUT_MODE == "single" and (not i and not j and not k)
                        ) or (OUTPUT_MODE == "depth" and (not i and not j)):
                            fp.write(get_var_header(model_data, var))
                            data_header.append(f"{DELIMITER}{var}")

                        # find the variable dimension ordering
                        if var not in lat_index:
                            print(
                                f"[INFO] Processing variable: {var_value}", flush=True
                            )
                            indices_list = ""
                            for dim_index, dim in enumerate(
                                model_data.variables[var].dimensions
                            ):
                                if (
                                    dim.encode("ascii", "ignore").decode("utf-8")
                                    == Z_VARIABLE
                                ):
                                    depth_index[var] = dim_index
                                    indices_list = f"{indices_list}{len(z)} {Z_VARIABLE}s index:{dim_index}  "
                                elif (
                                    dim.encode("ascii", "ignore").decode("utf-8")
                                    == X_VARIABLE
                                ):
                                    lon_index[var] = dim_index
                                    indices_list = f"{indices_list}{len(x)} {X_VARIABLE}s index:{dim_index}  "
                                elif (
                                    dim.encode("ascii", "ignore").decode("utf-8")
                                    == Y_VARIABLE
                                ):
                                    lat_index[var] = dim_index
                                    indices_list = f"{indices_list}{len(y)} {Y_VARIABLE}s index:{dim_index}  "
                                else:
                                    print(
                                        f'\n[ERR] Invalid dimensions "{dim}" in variable {var}'
                                    )
                                    sys.exit(2)
                            print(f"[INFO] {indices_list}", flush=True)
                            # A few checks.
                            if var not in lat_index and var not in lon_index:
                                print(
                                    f'\n[ERR] problem reading x and y variables "{X_VARIABLE}, {Y_VARIABLE}"'
                                )
                                sys.exit(2)

                            if var not in depth_index is None and var_ndim[var] > 2:
                                print(
                                    f'\n[ERR] problem reading the depth variables "{Z_VARIABLE}" for variable {var}'
                                )
                                sys.exit(2)

                            # Assign the variable's data values.
                            if var not in emcin:
                                try:
                                    emcin[var] = model_data.variables[var][:]
                                except Exception as err:
                                    print(f'\n[ERR] problem reading variable "{var}"')
                                    print("{0}\n".format(err))
                                    sys.exit(2)
                        if var in depth_index:
                            index[depth_index[var]] = k
                        index[lat_index[var]] = i
                        index[lon_index[var]] = j
                        # nan values, otherwise we write string to preserve the precision
                        if var_ndim[var] == 2:
                            if math.isnan(emcin[var][index[0]][index[1]]):
                                output_data.append(f"{DELIMITER}{math.nan}")
                            else:
                                # Conversion to string is done to preserve precision
                                output_data.append(
                                    f"{DELIMITER}{str(emcin[var][index[0]][index[1]]):s}"
                                )
                        else:
                            if math.isnan(emcin[var][index[0]][index[1]][index[2]]):
                                output_data.append(f"{DELIMITER}{math.nan}")
                            else:
                                # Conversion to string is done to preserve precision
                                output_data.append(
                                    f"{DELIMITER}{str(emcin[var][index[0]][index[1]][index[2]]):s}"
                                )

                output_data.append("\n")
        if OUTPUT_MODE == "depth":
            fp.write(f'{"".join(data_header)}\n')
            fp.write("".join(output_data))
            fp.close()

    if OUTPUT_MODE == "single":
        fp.write(f'{"".join(data_header)}\n')
        fp.write("".join(output_data))
        fp.close()


try:
    options, remainder = getopt.getopt(
        sys.argv[1:],
        "hdHi:v:x:y:z:m:",
        [
            "help",
            "input=",
            "variables=",
            "latitude=",
            "longitude=",
            "depth=",
            "mode=",
            "debug",
            "header",
        ],
    )
    for opt, arg in options:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt in ("-d", "--debug"):
            DEBUG = True
        elif opt in ("-H", "--header"):
            VIEW_HEADER = True
        elif opt in ("-i", "--input"):
            NETCDF_FILE_NAME = arg
        elif opt in ("-x", "--xvar"):
            X_VARIABLE = arg
        elif opt in ("-y", "--yvar"):
            Y_VARIABLE = arg
        elif opt in ("-z", "--zvar"):
            Z_VARIABLE = arg
        elif opt in ("-m", "--mode"):
            OUTPUT_MODE = arg
            if OUTPUT_MODE not in VALID_MODES:
                print(
                    f"[ERR] could not find the netCDF model file {NETCDF_FILE_NAME}",
                    flush=True,
                )
                usage()
                sys.exit(2)
except getopt.GetoptError:
    usage()
    sys.exit(2)

if not VIEW_HEADER or (VIEW_HEADER and (X_VARIABLE or Y_VARIABLE or Z_VARIABLE)):
    if X_VARIABLE is None or Y_VARIABLE is None:
        usage()
        print(f"[ERR] -x, -y, and -z values are required")
        sys.exit(2)

    if Z_VARIABLE is None:
        ndim = 2
        ndim_label = "2D"
    elif Z_VARIABLE.lower() == "none":
        Z_VARIABLE = None
        ndim = 2
        ndim_label = "2D"
    else:
        ndim = 3
        ndim_label = "3D"

print(f"[INFO] working on the {ndim_label} netCDF file: {NETCDF_FILE_NAME}\n")
make_model_geocsv()
print("...done")
