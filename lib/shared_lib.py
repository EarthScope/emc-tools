import os
import sys
import io
import logging
import time
from datetime import datetime, timezone
from pprint import pprint

from pyproj import Proj
import xarray as xr
import pandas as pd
from tqdm import tqdm

import numpy as np

# Get the directory paths.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop

import os
import stat
import time


def supports_unicode():
    # Check if the output supports Unicode \
    # and ensuring no exception occurs
    encoding = sys.stdout.encoding or ""
    return "UTF-8" in encoding.upper() or "UTF8" in encoding.upper()


def list_files_in_directory_from_file(file_path):
    """
    Lists the files in the directory containing the specified file and returns the message as a string.

    Parameters:
        file_path (str): The path to the file.

    Returns:
        str: A message with the list of files in the directory or an appropriate error message.
    """
    try:
        # Get the directory from the file path
        directory_path = os.path.dirname(file_path)
        if not directory_path.strip():
            directory_path = "."

        # Check if the directory exists
        if not os.path.isdir(directory_path):
            return f"The directory '{directory_path}' does not exist."

        # Get a list of files in the directory
        files = os.listdir(directory_path)

        # Filter out directories from the list
        files = [
            file for file in files if os.path.isfile(os.path.join(directory_path, file))
        ]

        # Create the message string
        if files:
            message = f"Files in '{directory_path}':\n"
            message += "\n".join(f" - {file}" for file in files)
        else:
            message = f"No files found in '{directory_path}'"

    except PermissionError:
        message = f"Permission denied to access the directory '{directory_path}'."
    except Exception as e:
        message = f"An error occurred: {str(e)}"

    return message


def file_info(file_path):
    """
    Shows detailed information about a file on the disk.

    Parameters:
        file_path (str): The path to the file.

    Returns:
        str: A string containing detailed information about the file.
    """
    try:
        # Get the file stats
        file_stats = os.stat(file_path)

        # Get the file mode (permissions)
        file_mode = stat.filemode(file_stats.st_mode)

        # Get the file size
        file_size = file_stats.st_size

        # Get the last access time
        last_access_time = time.ctime(file_stats.st_atime)

        # Get the last modification time
        last_modification_time = time.ctime(file_stats.st_mtime)

        # Get the last status change time
        last_status_change_time = time.ctime(file_stats.st_ctime)

        # Get the file owner's user ID
        owner_user_id = file_stats.st_uid

        # Get the file owner's group ID
        owner_group_id = file_stats.st_gid

        # Format the detailed information
        info = (
            f"File: {file_path}\n"
            f"Permissions: {file_mode}\n"
            f"Size: {file_size} bytes\n"
            f"Last accessed: {last_access_time}\n"
            f"Last modified: {last_modification_time}\n"
            f"Last status change: {last_status_change_time}\n"
            f"Owner user ID: {owner_user_id}\n"
            f"Owner group ID: {owner_group_id}"
        )

        return info

    except FileNotFoundError:
        return f"The file '{file_path}' does not exist."
    except PermissionError:
        return f"Permission denied to access the file '{file_path}'."
    except Exception as e:
        return f"An error occurred: {str(e)}"


def get_key_value(line):
    """split a line of model metadata to its line_type, key, value elements.

    Keyword arguments:
    line -- [required] a single input line.
    """

    line = line.strip()
    line_type = None
    key = None
    value = None
    if not line or line.startswith("#"):
        return line_type, key, value

    # This is really checking for either > or >>.
    if "=" not in line and line[0] != ">":
        logging.error(f"[ERR] did not find '=' in: {line}")
        sys.exit(2)

    if line[0] not in (">>", ">", "-"):
        line_type = ">>>"
    elif line.startswith(">>"):
        line_type = line[0:2]
        line = line[2:].strip()
    else:
        line_type = line[0]
        line = line[1:].strip()

    if line_type not in ("-", ">", ">>") and "=" not in line:
        logging.error(f"[ERR] did not find '=' in: {line}")
        sys.exit(2)

    if line_type in (">", ">>"):
        key = line.strip()
        value = None
    else:
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()

    return line_type, key, value


def read_model_metadata(model_file, params=dict()):
    """Get model parameters from a parameter file and convert it to a dictionary.

    Keyword arguments:
    model_file -- [required] the model metadata file
    """
    with open(model_file, "r") as fp:
        data = fp.read()
        lines = data.split("\n")

        group_key = None
        subgroup_key = None
        for line in lines:
            line_type, key, value = get_key_value(line)
            if line_type is None:
                continue
            elif line_type == ">":
                group_key = key
                subgroup_key = None
                params[group_key] = dict()
            if line_type == ">>" and group_key is None:
                logging.error(f"[ERR] found an orphaned subgroup: {line}")
                sys.exit(2)
            elif line_type == ">>":
                subgroup_key = key
                params[group_key][subgroup_key] = dict()
            elif line_type == ">>>" and group_key is None and subgroup_key is None:
                print(f"[ERR] item: {line} not in a group.subgroup")
                sys.exit(2)
            elif line_type == ">>>" and subgroup_key is None:
                params[group_key][key] = value
            elif line_type == ">>>":
                params[group_key][subgroup_key][key] = value
            elif line.startswith("-"):
                group_key = None
                subgroup_key = None

                params[key] = value

    return params


def cf_coordinate_names(x, y, aux=False):
    """Produce CF standard coordinate variable names.

    Keyword arguments:
    x -- [required] the current x-coordinate variable name
    y -- [required] the current y-coordinate variable name
    aux -- [default=False] flag to indicate if these variables are for the auxiliary coordinates.
    """
    if x.lower() == "longitude" and y.lower() == "latitude":
        return {"x": "longitude", "y": "latitude"}
    elif aux:
        return {"x": x, "y": y}
    elif x.lower() != "longitude" and y.lower() != "latitude":
        return {"x": "x", "y": "y"}
    else:
        logging.error(
            f"[cf_coordinate_names ERR] invalid coordinate name combinations ({x},{y})"
        )
        raise


def get_filename_without_extension(file_path):
    """
    Get the filename without the extension or path.

    Args:
        file_path (str): The full path of the file.

    Returns:
        str: The filename without its extension. If the filename is only an extension (e.g., ".csv"), an empty string is returned.

    Example:
        >>> get_filename_without_extension("/path/to/your/file.txt")
        'file'
        >>> get_filename_without_extension("/path/to/your/.csv")
        ''
    """
    base_name = os.path.basename(
        file_path
    )  # Get the base name (filename with extension)

    # Handle case where the file name is only an extension (e.g., ".csv")
    if base_name.startswith(".") and base_name.count(".") == 1:
        return ""

    file_name, _ = os.path.splitext(
        base_name
    )  # Split the base name into name and extension
    return file_name


def get_geocsv_metadata_from_ds(ds):
    """Compose GeoCSV style metadata for a given Dataset.

    Keyword arguments:
    ds -- [required] the xarray dataset
    """
    geocsv_metadata = list()
    for row in ds.attrs:
        # For some reason history lines are split (investigate).
        if row == "history":
            ds.attrs[row] = ds.attrs[row].replace("\n", ";")
        geocsv_metadata.append(f"# global_{row}: {ds.attrs[row]}")
    for var in ds.variables:
        if "variable" not in ds[var].attrs:
            geocsv_metadata.append(f"# {var}_variable: {var}")
            geocsv_metadata.append(f"# {var}_dimensions: {len(ds[var].dims)}")

        for att in ds[var].attrs:
            geocsv_metadata.append(f"# {var}_{att}: {ds[var].attrs[att]}")
            if att == "missing_value":
                geocsv_metadata.append(f"# {var}_missing_value: {prop.na_rep}")
            if att == "variable":
                geocsv_metadata.append(f"# {var}_dimensions: {len(ds[var].dims)}")
                geocsv_metadata.append(f"# {var}_column: {var}")
    metadata = "\n".join(geocsv_metadata)
    metadata = f"{metadata}\n"
    return metadata


def view_list(lst, indent=14, chunk=20):
    """Print a chunk from start and end of a list to view.

    Keyword arguments:
    lst -- [required] list of the numbers
    chunk -- [default: 20] print chunk elements from start and end of the list.
    """
    file = io.StringIO()
    # Not needed, but move pointer to start of file.
    file.seek(0)
    if len(lst) <= 60:
        pprint(
            list(lst),
            file,
            width=80,
            compact=True,
            indent=indent,
        )
    else:
        pprint(
            list(lst[0:chunk]),
            file,
            width=80,
            compact=True,
            indent=indent,
        )
        print(f"{indent*' '}...", file=file)
        pprint(
            list(lst[-1 * chunk :]),
            file,
            width=80,
            compact=True,
            indent=indent,
        )
    print(file.getvalue().replace("]", " ").replace("[", " "))


def standard_units(unit):
    """Check an input unit and return the corresponding standard unit."""
    unit = unit.strip().lower()
    if unit in ["m", "meter"]:
        return "m"
    elif unit in ["degrees", "degrees_east", "degrees_north"]:
        return "degrees"
    elif unit in ["km", "kilometer"]:
        return "km"
    elif unit in ["g/cc", "g/cm3", "g.cm-3", "grams.centimeter-3"]:
        return "g/cc"
    elif unit in ["kg/m3", "kh.m-3"]:
        return "kg/m3"
    elif unit in ["km/s", "kilometer/second", "km.s-1", "kilometer/s", "km/s"]:
        return "km/s"
    elif unit in ["m/s", "meter/second", "m.s-1", "meter/s", "m/s"]:
        return "m/s"
    elif unit.strip().lower in ["", "none"]:
        return ""


def is_numeric(in_string):
    """Checks an input string to make sure it only contains numeric characters"""
    return in_string.replace(".", "").replace("-", "").strip().isnumeric()


def is_in_range(value, value_type):
    """Checks a  value to make sure it is in range for the given value_type."""
    value_ranges = {"latitude": [-90.0, 90.0], "longitude": [-180.0, 180.0]}

    # Values without defined ranges are only checked for being a valid number.
    if value_type not in value_ranges:
        return is_numeric(f"{value}")
    else:
        # Make sure the value is a number.
        if not is_numeric(f"{value}"):
            return False
        else:
            range = value_ranges[value_type]

            return range[0] <= float(value) <= range[1]


def closest(lst, value):
    """Find the closest number in a list to the given value

    Keyword arguments:
    lst -- [required] list of the numbers
    value -- [required] value to find the closest list member for.
    """

    # Convert the list to a NumPy array
    lst_array = np.array(lst)

    # Compute the absolute differences
    diff = np.abs(lst_array - value)

    # Find the index of the minimum difference
    min_index = np.unravel_index(np.argmin(diff, axis=None), diff.shape)

    # Get the closest value
    closest_value = lst_array[min_index]
    # lst[min(range(len(lst)), key=lambda i: abs(lst[i] - value))]

    return closest_value


def utc_now():
    """Return the current UTC time."""
    try:
        _utc = datetime.now(tz=timezone.utc)
        utc = {
            "date_time": _utc.strftime("%Y-%m-%dT%H:%M:%S"),
            "datetime": _utc,
            "epoch": _utc.timestamp(),
        }
        return utc
    except:
        logging.error(f"[UTC ERR] Failed to get the current UTC time")
        raise


def check_file_type(file_path):
    """Check a given file's format

    Keyword arguments:
    file_path -- [default None] File to check."""
    file_type = dict()
    # Check for a valid GeoCSV file.
    try:
        with open(file_path, newline="") as csvfile:
            start = csvfile.read(300).strip()
            if start.startswith(f"# dataset:"):
                file_type["engine"] = "geocsv"
                version = start.split("\n")[0].split(":")[1].strip()
                if version == prop.dataset_version:
                    file_type["valid"] = True
                else:
                    file_type["valid"] = (
                        f"[ERR] Invalid {file_type['engine']} file type. For 'dataset', "
                        f"expected {prop.dataset_version}, but received {version}"
                    )
            else:
                file_type["engine"] = None
                file_type["valid"] = (
                    f"[ERR] For a GeoCSV file, the '# dataset: value' should always be present "
                    f"and should be the first line of a dataset."
                )
    except UnicodeDecodeError:
        # Check for a valid netCDF file.
        try:
            nc_dataset = xr.open_dataset(file_path, engine="netcdf4")
            file_type["engine"] = "netcdf"
            file_type["valid"] = True
        except Exception as ex:
            file_type["engine"] = None
            file_type["valid"] = (
                f"[ERR] dataset File {file_path} type not recognized.\n {ex}"
            )
    except Exception as ex:
        file_type["engine"] = None
        file_type["valid"] = f"[ERR] File {file_path} type not recognized.\n {ex}"

    return file_type


def get_logger(log_file=None, log_stream=None, config=True, log_level=0):
    """Set up the logging.

    Keyword arguments:
    log_file -- [default None] name of the log file to write to
    log_stream -- [default None] the log stream to log to
    config -- [default True] should configuring the logging module
    log_level -- [default 0] minimum priority level of messages to log (NOTSET=0, DEBUG=10,
                INFO=20, WARN=30, ERROR=40, CRITICAL=50)
    """
    root = logging.getLogger()
    logging.getLogger("matplotlib.font_manager").disabled = True
    if root.handlers:
        for handler in root.handlers:
            root.removeHandler(handler)
    if config or log_file is not None or log_stream is not None:
        if log_file is not None:
            logging.basicConfig(
                filename=log_file,
                format="%(message)s",
                encoding="utf-8",
                level=log_level,
            )
        elif log_stream is not None:
            logging.basicConfig(
                stream=log_stream,
                format="%(message)s",
                encoding="utf-8",
                level=log_level,
            )
        else:
            logging.basicConfig(format="%(message)s", encoding="utf-8", level=log_level)

    # Retrieve the logger instance
    logger = logging.getLogger()
    return logger


def get_module_param(params, var, required=True):
    """Extract a variable from a parameter module.

    Keyword arguments:
    param -- [required] the parameter module
    var -- [required] the variable to extract from the param object.
    required [default True] -- is var required? If True, the code will
    exit if not found.

    """

    # Set up the logger.
    logger = get_logger()
    if var not in dir(params):
        logger.error(
            f"[ERR] parameter '{var}' is required but it is missing from the parameter module."
        )
        if required:
            sys.exit(2)
        else:
            return None

    return getattr(params, var)


def get_param(params, var, required=True):
    """Extract a variable from a parameter dictionary.

    Keyword arguments:
    param -- [required] the parameter dictionary
    var -- [required] the variable to extract from the param object.
    required [default True] -- is var required? If True, the code will
    exit if not found.

    """

    # Set up the logger.
    logger = get_logger()
    if var not in params:
        if required:
            logger.error(
                f"[ERR] parameter '{var}' is required but it is missing from the parameter file.\n{list(params.keys())}"
            )
            sys.exit(2)
        else:
            return None

    return params[var]


def merge_dictionaries(parent_dict, dict_list):
    """Merge a list of dictionaries to the parent dictionary using update() method.

    Keyword arguments:
    parent_dict -- [required] The parent dictionary.
    dict_list -- a list of dictionaries to merge with the parent dictionary.
    """
    for _dic in dict_list:
        parent_dict.update(_dic)

    return parent_dict


def insert_after(parent_dict, dict_item, after_key, after=True):
    """Insert an item into a dictionary before or after a given item.

    Keyword arguments:
    parent_dict -- [required] The parent dictionary
    dict_item -- [required] The dictionary item to insert
    after_key -- [required] The dictionary key to insert the item after
    after -- [default True] if True, will insert after the given key. Otherwise, insert before it
    """
    delta = 0
    if after:
        delta = 1

    pos = list(parent_dict.keys()).index(after_key) + delta
    items = list(parent_dict.items())
    items.insert(pos, dict_item)

    return dict(items)


def set_variable_keys(variable, var):
    updated_var = dict()
    for key in variable:
        updated_var[f"{var}_{key}"] = variable[key]


def project_lonlat_utm(
    longitude, latitude, utm_zone, ellipsoid, xy_to_latlon=False, preserve_units=False
):
    """
    Performs cartographic transformations. Converts from longitude, latitude to UTM x,y coordinates
    and vice versa using PROJ (https://proj.org).

     Keyword arguments:
    longitude (scalar or array) – Input longitude coordinate(s).
    latitude (scalar or array) – Input latitude coordinate(s).
    xy_to_latlon (bool, default=False) – If inverse is True the inverse transformation from x/y to lon/lat is performed.
    preserve_units (bool) – If false, will ensure +units=m.
    """
    utm_zone = int(utm_zone)

    P = Proj(
        proj="utm",
        zone=utm_zone,
        ellps=ellipsoid,
    )
    # preserve_units=preserve_units,

    x, y = P(
        longitude,
        latitude,
        inverse=xy_to_latlon,
    )
    if np.isinf(x) or np.isinf(y):
        print(
            f"Invalid result for converting ({longitude}, {latitude}, inf values) with xy_to_latlon set to '{xy_to_latlon}' for the given UTM zone {utm_zone}."
        )
    return x, y


def time_it(t0):
    """Compute the elapsed time since last call

    Keyword arguments:
    t0 -- [required] Time of the last call
    """
    t1 = time.time()
    dt = t1 - t0
    if dt >= 60.0:
        text = f"Elapsed time: {dt / 60.0:0.1f}m"
    else:
        text = f"Elapsed time: {dt:0.2f}s"

    return time.time(), text
