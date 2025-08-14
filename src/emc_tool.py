#!/usr/bin/env python

import sys
import os
import getopt
import traceback
from pathlib import Path
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from textwrap import fill

from pyproj import Proj, Geod

"""Extract data from an EMC netCDF file. Allow user to inspect the metadata, slice the data, plot, and 
    save the sliced data. The slicing can be performed along the existing coordinate planes or as an interpolated cross-sectional slice through gridded data.

    Call arguments:
        -h, --help: this message.
        -i, --input: [required] the input nefiletCDF filename.
        -v, --verbose [optional] provide informative information during the run.
"""

script = os.path.basename(sys.argv[0])
script = Path(script).stem
# Get the directory paths.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop
import tool_prop as slicer_prop

LIB_DIR = os.path.join(ROOT_DIR, prop.LIB_DIR)
sys.path.append(LIB_DIR)
import shared_lib as lib
import tool_lib as slicer_lib

dash = 10 * "-"

# Set up the logger.
logger = lib.get_logger()
valid_file_types = {"csv": ".csv", "geocsv": ".csv", "netcdf": ".nc"}


def usage():
    logger.info(
        f"""
            A tool to extract data from an EMC netCDF file. It allows to interactively inspect metadata, s
            lice data, plot results, and save the sliced output.
            
            Currently, slicing is supported only along existing coordinate planes.

            Command-line arguments:
                -h, --help  Display this help message.
                -v, --verbose [Optional] Enable verbose output.
                -i, --input  [Required] Path to the input netCDF file.

        """
    )


import os


def custom_formatter(x):
    """
    Custom formatter function
    If the value is close enough to zero (including -0.0), format it as '0'.
    Otherwise, use the default formatting.
    """
    if abs(x) < 1e-12:  # 1e-12 is used as a threshold for floating-point comparison
        return "0"
    else:
        return f"{x}"


def format_numbers(numbers, digits=2):
    """
    Formats a single number or a list/tuple of numbers to a specified number of decimal places.

    Parameters:
        numbers: A single number, or a list/tuple of numbers to be formatted.
        digits: The number of decimal places to format to. Default is 2.

    Returns:
        A string with the formatted number(s).
    """
    if isinstance(numbers, (int, float)):
        # If it's a single number, format it directly
        formatted_numbers = f"{numbers:.{digits}f}"
    elif isinstance(numbers, (list, tuple)):
        # If it's a list or tuple, format each number individually
        formatted_numbers = f"[{', '.join(f'{num:.{digits}f}' for num in numbers)}]"
    else:
        raise TypeError("Input must be a single number, a list, or a tuple.")

    return formatted_numbers


def get_netcdf_engine(filename):
    """Finds the netcdf engine for reading a netcdf file."""
    for extension in prop.netcdf_engines:
        if filename.endswith(extension):
            return prop.netcdf_engines[extension]
    return None


def display_range(ds, ds_type):
    """Output ranges for the given dataset."""
    logger.info(f"\n{ds_type} netCDF file:")
    logger.info(f"\n[Range] Coordinate Variables:")
    for var in list(ds.dims):
        if "units" not in ds[var].attrs:
            logger.error(f"Missing units for variable '{var}'")
            sys.exit(2)

        logger.info(
            f"\t{var}: {np.nanmin(ds[var].data):0.2f} to  {np.nanmax(ds[var].data):0.2f} {ds[var].attrs['units']}"
        )
    logger.info(f"\n[Range] Data Variables:")
    for var in list(ds.data_vars):
        if "units" not in ds[var].attrs:
            logger.error(f"Missing units for variable '{var}'")
            sys.exit(2)

        logger.info(
            f"\t{var}: {np.nanmin(ds[var].data):0.2f} to  {np.nanmax(ds[var].data):0.2f} {ds[var].attrs['units']}"
        )
    logger.info("\n")


def _to_python_type(val):
    """Convert NumPy scalars/arrays to Python native types."""
    if isinstance(val, np.generic):
        return val.item()  # scalar → Python float/int
    elif isinstance(val, np.ndarray):
        return val.tolist()  # array → Python list
    return val


def _print_coord_values(ds, names, indent):
    """Print coordinate variable metadata with Python-native floats."""
    for name in names:
        logger.info(f"\t{name}: ")
        # Print attributes
        for attr_name, attr_val in ds[name].attrs.items():
            logger.info(f"{' ' * indent}{attr_name}: {attr_val}")
        # Convert values to plain floats
        arr = ds[name].values
        if isinstance(arr, np.ndarray):
            vals = np.asarray(arr).astype(float).ravel().tolist()
        else:
            vals = [float(arr)]
        # Nicely format the values
        val_str = ", ".join(f"{v:g}" for v in vals)
        wrapped = fill(val_str, width=88, subsequent_indent=" " * indent)
        logger.info(f"{' ' * indent}Values:")
        logger.info(f"{' ' * indent}{wrapped}")


def display_metadata(ds, ds_type):
    """Output metadata for the given dataset."""
    logger.info(f"\n[Metadata] {ds_type} netCDF file\n")

    # Global attributes
    logger.info(f"\n[Metadata] Global attributes:\n")
    for row in ds.attrs:
        if "geospatial" in row:
            val = _to_python_type(ds.attrs[row])
            logger.info(f"\t{row}: {val}")

    # Coordinate Variables (converted to Python floats)
    logger.info(f"\n[Metadata] Coordinate Variables:")
    indent = 14
    _print_coord_values(ds, list(ds.dims), indent)

    # Data Variables (original slicer_lib behavior, no values printed)
    logger.info(f"\n[Metadata] Data Variables:")
    slicer_lib.display_var_meta(ds, list(ds.data_vars), indent, values=False)
    logger.info("\n")


def output_messages(messages):
    """Output messages stored in a message list."""
    if messages:
        new_line = "\n"
        logger.info(f"{new_line}NOTES: {new_line.join(messages)}")
    return list()


def output_data(subset_volume, filename, messages):

    # Define the desired precision
    precision = 3

    # Data variables.

    data_vars = list(subset_volume.data_vars)

    # The main coordinate variables.
    main_coords_list = [coord for coord in subset_volume.dims]

    # Round the auxiliary coordinates to the desired precision
    aux_coords_list = [
        coord for coord in subset_volume.coords if coord not in main_coords_list
    ]

    for var in aux_coords_list:
        if var in ["latitude", "longitude"]:
            subset_volume[var].values = np.round(subset_volume[var].values, precision)
        else:
            subset_volume[var].values = np.round(subset_volume[var].values, precision)

    if filename.endswith(".gcsv"):
        meta = lib.get_geocsv_metadata_from_ds(subset_volume)
        meta = f"# dataset: GeoCSV2.0\n# delimiter: ,\n{meta}"
        with open(filename, "w") as fp:
            fp.write(meta)
        subset_volume.to_dataframe().to_csv(filename, mode="a")
        messages.append(f"[INFO] Saved as {filename}")
    elif filename.endswith(".csv"):
        meta = lib.get_geocsv_metadata_from_ds(subset_volume)
        with open(filename, "w") as fp:
            subset_volume.to_dataframe().to_csv(filename, mode="w")
        messages.append(f"[INFO] Saved as {filename}")
    elif filename.endswith(".nc"):
        # Generate the encoding dictionary for all variables to apply compression
        encoding = {var: {"zlib": True, "complevel": 4} for var in data_vars}

        # To prevent _FillValue from being automatically assigned to dimension variables.
        for coord in main_coords_list:
            if coord not in encoding:
                encoding[coord] = {"_FillValue": None}
            else:
                encoding[coord]["_FillValue"] = None

        subset_volume.to_netcdf(
            filename, mode="w", format="NETCDF4_CLASSIC", encoding=encoding
        )
        messages.append(f"[INFO] Saved as {filename}")
    elif filename.endswith(".h5"):
        # Generate the encoding dictionary for all variables to apply compression
        encoding = {var: {"zlib": True, "complevel": 4} for var in data_vars}

        # To prevent _FillValue from being automatically assigned to dimension variables.
        for coord in main_coords_list:
            if coord not in encoding:
                encoding[coord] = {"_FillValue": None}
            else:
                encoding[coord]["_FillValue"] = None

        subset_volume.to_netcdf(
            filename, mode="w", engine=prop.netcdf_engines[".h5"], encoding=encoding
        )
        messages.append(f"[INFO] Saved as {filename}")
    else:
        messages.append(f"[ERR] invalid file: {filename}")
    return messages


def determine_dataset_type_by_data_vars(dataset):
    """
    Determine whether the dataset is 2D or 3D based on the number of dimensions of its data variables.

    Parameters:
        dataset (xarray.Dataset): The xarray dataset to evaluate.

    Returns:
        str: "2D" if the dataset contains primarily two-dimensional data variables,
             "3D" if three-dimensional data variables are present,
             or "Unknown" if no clear determination can be made.
    """
    # Initialize counters for 2D and 3D variables
    is_2d = False
    is_3d = False

    # Iterate over data variables and check their dimensionality
    for var_name, data_array in dataset.data_vars.items():
        num_dims = len(
            data_array.dims
        )  # Count the number of dimensions for the variable

        if num_dims == 2:
            is_2d = True
        elif num_dims == 3:
            is_3d = True

        # If both 2D and 3D are found, prioritize 3D as the dataset type
        if is_3d:
            return "3D"

    # Determine dataset type based on findings
    if is_2d:
        return "2D"
    else:
        return "Unknown"


def main():
    vmin = None
    vmax = None
    messages = list()
    cmap = prop.cmap
    interpolation_method = slicer_prop.interpolation_method
    xsection_steps = slicer_prop.steps
    vertical_exaggeration = 0

    # Capture the input parameters.
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(
            argv, "hvmi:o:", ["help", "verbose", "meta", "input=", "output="]
        )
    except getopt.GetoptError as err:
        # Print the error, and help information and exit:
        logger.error(err)
        usage()
        sys.exit(2)
    # Initialize the variables.
    input_file = None
    verbose = False
    meta = None
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-v", "--verbose"):
            verbose = True
        elif o in ("-m", "--meta"):
            meta = True
        elif o in ("-i", "--input"):
            input_file = a
            if not os.path.isfile(input_file):
                logger.error(
                    f"[ERR] Invalid input netCDF file: [{input_file}]. File not found!"
                )
                sys.exit(2)
        else:
            assert False, "unhandled option"

    # The input filename is required.
    if input_file is None:
        usage()
        logger.error("[ERR] missing -i or --input.")
        sys.exit(1)

    logger.info(f"[INFO] Working on input {input_file}")

    # Figure out the input's file type.
    file_type = lib.check_file_type(input_file)

    option = "start"
    coordinate_values = dict()

    # Convert netCDF to GeoCSV.
    logger.info("\n")
    if file_type["engine"] == "netcdf" and file_type["valid"] == True:
        if verbose:
            logger.info(
                f"\n\n{script}\n{dash}\n"
                f"Tool for interactively extracting data from an EMC netCDF file. Users can inspect metadata, slice data, \n"
                f"plot results, and save the output. Slicing can be done along existing coordinate planes or for 3D models, as an \n"
                f"interpolated cross-sectional slice through gridded data.\n\n"
            )

        engine = get_netcdf_engine(input_file)
        with xr.open_dataset(input_file, engine=engine) as ds:
            base_title = ""
            if "model" in ds.attrs:
                base_title = ds.attrs["model"]
            elif "id" in ds.attrs:
                base_title = ds.attrs["id"]

            # Determine if the dataset is 2D or 3D
            dataset_type = determine_dataset_type_by_data_vars(ds)
            messages.append(
                f"[INFO] Loaded {input_file} and it is {dataset_type} netCDF file"
            )

            data_var = list(ds.data_vars)
            coordinates = list(ds.coords)
            try:
                for var in coordinates:
                    if ds[var].dims == ():
                        logger.warning(
                            f"[WARN] Missing data array for variable {var}: {ds[var]}\nData: {ds[var].data},\nskipped"
                        )
                        continue
                    if var not in ds:
                        logger.error(f"[ERR] Missing data for variable {var}")
                        sys.exit(3)
                    elif ds[var].data.size < 1:
                        logger.warning(
                            f"[WARN] Missing data for variable {var}, will not use"
                        )
                        # In spacial cases, like information on the axis-rotation, we may have variables with no data.
                        coordinates.remove(var)
                        ds = ds.drop_dims(var)
                        continue

                    coordinate_values[var] = list(ds[var].data)
            except Exception as ex:
                logger.error(f"[ERR] failed to get data for var:{var}\n{ex}")
                sys.exit(2)

            # Identify main coordinates (x, y, latitude, longitude)
            main_coords_list = ["x", "y", "latitude", "longitude"]

            # Identify main and auxiliary coordinates
            main_coords = [coord for coord in ds.coords if coord in ds.dims]
            aux_coords_list = [coord for coord in ds.coords if coord not in ds.dims]
            aux_coords = aux_coords_list.copy()
            # Make sure other coordinates not in the main list appear in aux_coords.
            aux_coords += [
                coord for coord in main_coords_list if coord not in main_coords
            ]

            # Make sure the aux coord items are not included in data_var.
            for _aux in aux_coords:
                if _aux in data_var:
                    data_var.remove(_aux)

            while option != "exit":
                if verbose:
                    logger.info(
                        f"\nThe available options are:\n\tmeta - to view file's metadata\n\trange - to display value ranges for variables\n\tsubset - to subset the data\n\t{dash}\n\thelp\n\texit "
                    )
                messages = output_messages(messages)
                option = input(
                    "[data] select option [meta, range, subset, help, exit]? "
                )

                # Done, back!
                if (
                    option.strip() == "exit"
                    or option.strip() == "back"
                    or not option.strip()
                ):
                    sys.exit()

                # Need help.
                elif option == "help":
                    usage()

                # Variables value range.
                elif option == "range":
                    display_range(ds, dataset_type)

                # Display the metadata.
                elif option == "meta":
                    display_metadata(ds, dataset_type)

                # Subset the model.
                elif option == "subset":
                    subset_type = "volume"
                    # Actions.
                    while subset_type != "back":
                        # f"\nSubset type [volume, slice, xsection, back, exit]? "

                        if dataset_type == "2D":
                            if verbose:
                                logger.info(
                                    f"\nYou can subset the data as a:"
                                    f"\n\tvolume - a subvolume of data"
                                    f"\n\tsurface - a surface plot of a 2D variable"
                                    f"\n\t{dash}"
                                    f"\n\tback - takes you to the previous step\n\texit "
                                )
                            selection_list = ["volume", "surface", "back", "exit"]
                            subset_type = input(
                                f"[subset] select [volume, surface, back, exit]? "
                            )
                        else:
                            if verbose:
                                logger.info(
                                    f"\nYou can subset the data as a:"
                                    f"\n\tvolume - a subvolume of data"
                                    f"\n\tslice - a slice along a coordinate axes"
                                    f"\n\txsection - a vertical slice in an arbitrary direction"
                                    f"\n\t{dash}"
                                    f"\n\tback - takes you to the previous step\n\texit "
                                )
                            selection_list = [
                                "volume",
                                "slice",
                                "xsection",
                                "back",
                                "exit",
                            ]
                            subset_type = input(
                                f"[subset] select [volume, slice, xsection, back, exit]? "
                            )

                        # Did user make a valid selection?
                        if subset_type not in selection_list:
                            continue

                        # Done, back!
                        if subset_type.strip() == "exit":
                            sys.exit()
                        elif subset_type.strip() == "back" or not option.strip():
                            break

                        if subset_type == "volume":
                            subtitle = "volume"
                            # Get the volume limits.
                            subset = ds.copy()
                            subset_limits = dict()
                            for dim in subset.dims:
                                subset_limits[dim] = (
                                    np.nanmin(subset[dim].data),
                                    np.nanmax(subset[dim].data),
                                )
                                if verbose:
                                    logger.info(
                                        f"\nEnter the {dim} range for the volume as minimum,maximum.\n\tPress return to accept the default values (full range) \n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                    )
                                messages = output_messages(messages)
                                # Format each number in the list to two decimal places and join them with a comma
                                _limits = input(
                                    f"[subset-volume] {dim} range [default values are: {subset_limits[dim]}, back, exit]? "
                                )

                                # Done, back!
                                if _limits.strip() == "exit":
                                    sys.exit()
                                elif _limits.strip() == "back":
                                    break
                                elif _limits and "," not in _limits:
                                    messages.append(
                                        f"[WARN] Invalid limits for {dim} ({_limits}), using the full range of {subset_limits[dim]}."
                                    )
                                    logger.warning(messages[-1])
                                    _limits = ""
                                if _limits:
                                    values = _limits.split(",")
                                    if lib.is_numeric(values[0]) and lib.is_numeric(
                                        values[1]
                                    ):
                                        values[0] = float(values[0].strip())
                                        values[1] = float(values[1].strip())
                                        if lib.is_in_range(
                                            values[0], dim
                                        ) and lib.is_in_range(values[1], dim):
                                            subset_limits[dim] = (
                                                min(values),
                                                max(values),
                                            )
                                            messages.append(
                                                f"[INFO] {dim} limits: {values}."
                                            )
                                        else:
                                            messages.append(
                                                f"[WARN] Invalid limits of {_limits} for {dim}, using the full range of {subset_limits[dim]}."
                                            )
                                            logger.warning(messages[-1])
                                            _limits = ""
                                    else:
                                        messages.append(
                                            f"[WARN] Invalid limits of {_limits} for {dim} (both limits must be provided as numbers). Will be using the full range of {subset_limits[dim]}."
                                        )
                                        logger.warning(messages[-1])
                                        _limits = ""
                                else:
                                    messages.append(
                                        f"[Info] {dim} will be set to full range: {subset_limits[dim]}."
                                    )
                                    logger.warning(messages[-1])
                                    _limits = ""
                            if _limits.strip() == "back":
                                break
                            # warnings: Checking for all NaN values
                            subset_volume, warnings = slicer_lib.subsetter(
                                ds, subset_limits, dataset_type
                            )

                            if not warnings:
                                # Subset the data.
                                messages.append(
                                    f"\nThe selected volume information:\n{subset_volume}"
                                )
                                subset_action = "continue"
                            else:
                                messages.append(
                                    f"Invalid limits - {warnings} - please try again!"
                                )
                                subset_action = "back"
                            # Actions.
                            messages = output_messages(messages)
                            while subset_action != "back":
                                # Save the subset data.
                                subset_action == "save"
                                if verbose:
                                    logger.info(
                                        f"\nSave the data\n\tFilename -- The name of the file to\n\tFiletype -- The output file type must be one of: {list(valid_file_types.keys())}.\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                    )

                                messages = output_messages(messages)
                                file_info = input(
                                    f"[subset-volume] Output filename and type (filename,type)? [valid file types are: {list(valid_file_types.keys())}; back, exit]: "
                                )

                                # Done, back!
                                if file_info.strip() == "exit":
                                    sys.exit()
                                elif file_info.strip() == "back":
                                    break
                                elif "," not in file_info:
                                    messages.append(f"[ERR] Bad input: {file_info}")
                                    continue
                                else:
                                    try:
                                        filename, file_type = file_info.split(",")
                                        filename = filename.strip()
                                        file_type = file_type.strip()
                                    except Exception as ex:
                                        messages.append(
                                            f"[ERR] Bad input: {file_info}\n{ex}"
                                        )
                                        continue

                                    if file_type not in valid_file_types:
                                        messages.append(f"[ERR] Bad input: {file_info}")
                                        continue
                                    else:
                                        filename = f"{filename.strip()}{valid_file_types[file_type]}"

                                messages = output_data(
                                    subset_volume, filename, messages
                                )
                                messages = output_messages(messages)

                        elif subset_type == "xsection":
                            # A cross-section of the model.
                            subtitle = f"cross-section"
                            if verbose:
                                logger.info(
                                    f"[INFO] Plotting an interpolated cross-sectional slice through gridded data\nPlease provide the cross-section limits"
                                )
                            logger.info(
                                "\nThe model coordinate ranges for the cross-section are:"
                            )
                            xvar = "longitude"
                            yvar = "latitude"
                            zvar = "depth"
                            coords_list = list(ds.coords.keys())

                            # Remove the auxiliary coordinates.
                            if "northing" in coords_list:
                                coords_list.remove("northing")
                            if "easting" in coords_list:
                                coords_list.remove("easting")

                            do_warn = True
                            if xvar not in ds.coords or yvar not in ds.coords:
                                logger.err(
                                    f"[ERR] x-section is only only available from model geographic coordinates!"
                                )
                                continue

                            coords_list.remove(xvar)
                            coords_list.remove(yvar)

                            if zvar not in ds.coords:
                                if do_warn:
                                    logger.info(
                                        f"[WAR] Plotting a non-geospatial cross-section"
                                    )
                                    do_warn = False
                                _z = zvar
                                while _z not in coords_list:
                                    if _z == "back":
                                        continue
                                    elif _z == "exit":
                                        sys.exit()

                                    _z = input(
                                        f"[xsection] Z-variable? [valid variables are: {coords_list}; back, exit]: "
                                    )
                                zvar = _z
                            coords_list.remove(zvar)
                            for var in (xvar, yvar, zvar):
                                if "units" not in ds[var].attrs:
                                    logger.error(f"Missing units for variable '{var}'")
                                    sys.exit(2)

                                logger.info(
                                    f"\t{var}: {np.nanmin(ds[var].data):0.2f} to  {np.nanmax(ds[var].data):0.2f} {ds[var].attrs['units']}"
                                )

                            # Get the unit of the depth
                            _z_unit = ds[zvar].attrs.get(
                                "units", "No unit attribute found"
                            )
                            z_unit = lib.standard_units(_z_unit)
                            if zvar == "depth":
                                units = "cgs"
                                if z_unit == "km":
                                    units = "mks"

                                logger.info(
                                    f"[INFO] The {zvar} unit is detected as: {_z_unit}, will use: {z_unit}. Units set distance units based on: {units}"
                                )
                            else:
                                units = "mks"

                            start = slicer_lib.get_point("start", (xvar, yvar, zvar))
                            if start == "back":
                                subset_type = "back"
                                break
                            end = slicer_lib.get_point("end", (xvar, yvar, zvar))
                            if end == "back":
                                subset_type = "back"
                                break
                            depth_range = slicer_lib.get_point("z", (xvar, yvar, zvar))
                            if depth_range == "back":
                                subset_type = "back"
                                break

                            try:
                                depth = [
                                    min(float(depth_range[0]), float(depth_range[1])),
                                    max(float(depth_range[0]), float(depth_range[1])),
                                ]
                                if zvar == "depth":
                                    input_depth = depth
                                else:
                                    input_depth = [
                                        float(depth_range[0]),
                                        float(depth_range[1]),
                                    ]

                            except Exception as ex:
                                logger.error(f"[ERR] Bad depth values {depth}\n{ex}")
                                subset_type = "back"
                                break

                            plot_data = ds.copy()
                            plot_data = plot_data.where(
                                (plot_data[zvar] >= float(depth[0]))
                                & (plot_data[zvar] <= float(depth[1])),
                                drop=True,
                            )
                            utm_zone = None
                            meta = ds.attrs
                            if "grid_ref" not in meta:
                                messages.append(
                                    f"[WARN] The 'grid_ref' attribute not found. Assuming geographic coordinate system"
                                )
                                grid_ref = "latitude_longitude"
                            else:
                                grid_ref = meta["grid_ref"]
                                messages.append(f"[INFO] grid_ref is {grid_ref}")

                            # Cross-section interpolation type.
                            interp_type = interpolation_method[0]
                            if verbose:
                                logger.info(
                                    f"\nSelect the interpolation method for creating the cross-section. The {', '.join(interpolation_method)} methods are available.\n\tPress return to select the default {interp_type} method\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                )
                            messages = output_messages(messages)
                            _interp_type = input(
                                f"[xsection] Interpolation Method [{', '.join(interpolation_method+['back', 'exit'])}, default: {interp_type}]? "
                            )
                            if _interp_type.strip() == "exit":
                                sys.exit()
                            elif _interp_type.strip() == "back":
                                break
                            elif _interp_type.strip():
                                interp_type = _interp_type
                            elif interp_type not in (interpolation_method):
                                messages.append(
                                    f"[WARN] Invalid interpolation method of {interp_type}. Will use {interpolation_method[0]}"
                                )
                                interp_type = interpolation_method[0]

                            # Steps in the cross-section.
                            steps = xsection_steps
                            if verbose:
                                logger.info(
                                    f"\nThe number of points along the geodesic between the start and the end point (including the end points) to use in the cross section.\n\tPress enter to select the default of 100.\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                )
                            messages = output_messages(messages)
                            steps = input(
                                f"[xsection] Number of points ['back', 'exit', default: {steps}]? "
                            )
                            if not steps.strip():
                                steps = xsection_steps
                            elif steps.strip() == "exit":
                                sys.exit()
                            elif steps.strip() == "back":
                                break
                            elif not steps.isnumeric():
                                messages.append(
                                    f"[WARN] Invalid number of steps, expected an integer. Will use {xsection_steps}"
                                )
                                steps = xsection_steps
                            else:
                                steps = int(steps)

                            if verbose:
                                logger.info(
                                    f"\nVertical exaggeration for the cross section.\n\tPress enter to select the default of {vertical_exaggeration}.\n\t{dash}\n\tback - takes you to the previous step\n\texit "
                                )
                            messages = output_messages(messages)

                            _exagg = 0
                            while _exagg != "back":
                                _exagg = input(
                                    f"[xsection] Vertical exaggeration (0: dynamic aspect ratio)['back', 'exit', default: {vertical_exaggeration}]? "
                                )
                                relabel_y = False
                                if not _exagg.strip():
                                    exaggeration = vertical_exaggeration
                                elif _exagg.strip() == "exit":
                                    sys.exit()
                                elif _exagg.strip() == "back":
                                    _exagg = "back"
                                    continue
                                elif not _exagg.replace(".", "").isnumeric():
                                    messages.append(
                                        f"[WARN] Invalid vertical exaggeration, expected a number. Will use {vertical_exaggeration}"
                                    )
                                    exaggeration = vertical_exaggeration
                                else:
                                    exaggeration = float(_exagg)
                                    relabel_y = True
                                # Extract the cross-section.
                                try:
                                    xsection_data, latitudes, longitudes = (
                                        slicer_lib.interpolate_path(
                                            plot_data,
                                            start,
                                            end,
                                            num_points=steps,
                                            method=interp_type,
                                            grid_ref=grid_ref,
                                            utm_zone=meta["utm_zone"],
                                            ellipsoid=meta["ellipsoid"],
                                        )
                                    )
                                except Exception as ex:
                                    message = f"[ERR] interpolate_path failed: {ex}\n{traceback.print_exc()}"
                                    logger.error(message)
                                    break

                                # Calculate distances between consecutive points
                                geod = Geod(ellps=meta["ellipsoid"])
                                _, _, distances = geod.inv(
                                    longitudes[:-1],
                                    latitudes[:-1],
                                    longitudes[1:],
                                    latitudes[1:],
                                )

                                # Compute cumulative distance, starting from 0
                                cumulative_distances = np.concatenate(
                                    ([0], np.cumsum(distances))
                                )

                                if "grid_ref" not in meta:
                                    messages.append(
                                        f"[WARN] The 'grid_ref' attribute not found. Assuming geographic coordinate system"
                                    )
                                    grid_ref = "latitude_longitude"

                                else:
                                    grid_ref = meta["grid_ref"]

                                if units == "mks":
                                    cumulative_distances = cumulative_distances / 1000.0
                                    distance_label = "distance (km)"
                                else:
                                    distance_label = "distance (m)"

                                # Create a new coordinate 'distance' based on the cumulative distances
                                xsection_data = xsection_data.assign_coords(
                                    distance=("points", cumulative_distances)
                                )

                                # If you want to use 'distance' as a dimension instead of 'index',
                                # you can swap the dimensions (assuming 'index' is your current dimension)
                                xsection_data = xsection_data.swap_dims(
                                    {"points": "distance"}
                                )

                                # Iterate through the model variables and plot each cross-section.
                                plot_var = data_var[0]
                                while plot_var:

                                    if len(data_var) > 1:
                                        if verbose:
                                            logger.info(
                                                f"\nThe model variable:\n\tback - takes you to the previous step\n\texit "
                                            )
                                        messages = output_messages(messages)
                                        plot_var = input(
                                            f"[xsection] Variable to plot {data_var}, back, exit]: "
                                        )
                                    # Done, back!
                                    if plot_var.strip() == "exit":
                                        sys.exit()
                                    elif (
                                        plot_var.strip() == "back" or not option.strip()
                                    ):
                                        break
                                    elif plot_var not in data_var:
                                        messages.append(
                                            f"[ERR] Invalid input {plot_var}"
                                        )
                                        continue
                                    """
                                    _action = "plot"
                                    while _action != "back":

                                        if verbose:
                                            logger.info(
                                                f"\nWhat to do with the data:\n\tback - takes you to the previous step\n\texit "
                                            )
                                        messages = output_messages(messages)
                                        _action = input(
                                            f"[plot, save, back, exit]: "
                                        )
                                        # Done, back!
                                        if _action.strip() == "exit":
                                            sys.exit()
                                        elif _action.strip() == "back" or not option.strip():
                                            _action="break"
                                            continue
                                        elif _action == "plot":
                                        """
                                    pdata = xsection_data.copy()
                                    if zvar == "depth":
                                        positive = ds[zvar].attrs.get(
                                            "positive", "No positive attribute found"
                                        )

                                        if positive not in ("down", "up"):
                                            messages.append(
                                                f"[WARN] For {zvar} the 'positive' attribute is set to {positive}. It must be set to 'down' or 'up'. We set it here to 'down'"
                                            )
                                            z_factor = -1
                                        elif positive == "down":
                                            messages.append(
                                                f"[INFO] For {zvar} the 'positive' attribute is set to {positive}"
                                            )
                                            z_factor = -1
                                        else:
                                            messages.append(
                                                f"[INFO] For {zvar} the 'positive' attribute is set to {positive}"
                                            )
                                            z_factor = 1

                                        pdata[zvar] = z_factor * pdata[zvar]
                                        messages = output_messages(messages)

                                    if vmin and vmax:
                                        pdata[plot_var].plot.contourf(
                                            cmap=cmap, vmin=vmin, vmax=vmax
                                        )
                                    else:
                                        pdata[plot_var].plot.contourf(
                                            cmap=cmap,
                                        )
                                    # Get the current axes
                                    ax = plt.gca()
                                    ax.set_xlabel(distance_label)
                                    if relabel_y:
                                        # Retrieve the current y-axis label
                                        current_label = ax.get_ylabel()

                                        # Append new text to the existing label
                                        new_label = (
                                            f"{current_label}; VE x{exaggeration}"
                                        )

                                        # Set the new label
                                        ax.set_ylabel(new_label)
                                        if exaggeration > 0:
                                            plt.gca().set_aspect(exaggeration)

                                    # Set the depth limits for display.
                                    if zvar == "depth":
                                        if z_factor < 0:
                                            plt.ylim(
                                                z_factor * depth[1], z_factor * depth[0]
                                            )
                                        else:
                                            plt.ylim(input_depth[0], input_depth[1])
                                    else:
                                        z_factor = 1
                                        plt.ylim(input_depth[0], input_depth[1])

                                    # plt.gca().invert_yaxis()  # Invert the y-axis to show depth increasing downwards
                                    # Getting current y-axis tick labels
                                    labels = [
                                        item.get_text()
                                        for item in plt.gca().get_yticklabels()
                                    ]
                                    y_ticks = plt.gca().get_yticks()

                                    # Assuming the labels are numeric, convert them to float, multiply by -1, and set them back.
                                    # If labels are not set or are custom, you might need to adjust this part.
                                    new_labels = [
                                        (
                                            custom_formatter(
                                                z_factor
                                                * float(label.replace("−", "-"))
                                            )
                                            if label
                                            else 0.0
                                        )
                                        for label in labels
                                    ]  # Handles empty labels as well

                                    # Setting new labels ( better to explicitly set both the locations of the ticks and their labels using
                                    # set_yticks along with set_yticklabels.)
                                    plt.gca().set_yticks(y_ticks)
                                    plt.gca().set_yticklabels(new_labels)

                                    plt.xticks(
                                        rotation=45
                                    )  # Rotating the x-axis labels to 45 degrees

                                    # Adding vertical text for start and end locations
                                    plt.text(
                                        cumulative_distances[0],
                                        1.05,
                                        f"⟸{start}",
                                        rotation=90,
                                        transform=plt.gca().get_xaxis_transform(),
                                        verticalalignment="bottom",
                                        horizontalalignment="center",
                                        fontsize=9,
                                    )
                                    plt.text(
                                        cumulative_distances[-1],
                                        1.05,
                                        f"⟸{end}",
                                        rotation=90,
                                        transform=plt.gca().get_xaxis_transform(),
                                        verticalalignment="bottom",
                                        horizontalalignment="center",
                                        fontsize=9,
                                    )

                                    # Add title if the model name is available.
                                    plt.title(f"{base_title}\n{subtitle} of {plot_var}")

                                    # Adjust the layout
                                    plt.tight_layout()
                                    plt.show()
                                    break
                                    """
                                        # Save the slice data.
                                        elif _action == "save":

                                            filename = ""
                                            while filename != "back":

                                                messages = output_messages(messages)
                                                file_info = input(f"Output filename and type (filename,type)? [valid file types are: {list(valid_file_types.keys())}; back, exit]: ")

                                                # Done, back!
                                                if file_info.strip() == "exit":
                                                    sys.exit()
                                                elif file_info.strip() == "back":
                                                    break
                                                elif "," not in file_info:
                                                    messages.append(f"[ERR] Bad input: {file_info}")
                                                    continue
                                                else:
                                                    filename, file_type = file_info.split(",")
                                                    if file_type.strip() not in valid_file_types:
                                                        messages.append(f"[ERR] Bad input: {file_info}")
                                                        continue
                                                    else:
                                                        filename = f"{filename.strip()}{valid_file_types[file_type]}"
                                                messages = output_data(xsection_data[plot_var], filename, messages)
                                                messages = output_messages(messages)
                                        """

                        elif subset_type in ["slice", "surface"]:
                            # Slice the model.
                            slice_dir = "depth"
                            while slice_dir:

                                if subset_type == "slice":
                                    if verbose:
                                        logger.info(
                                            f"\nA slice cuts the model along one of the coordinate axes."
                                            f"\n\tdirection - direction of the slice, the coordinate to cut the model along"
                                            f"\n\t{dash}\n\tback - takes you to the previous step\n\texit  "
                                        )
                                    messages = output_messages(messages)
                                    slice_dir = input(
                                        f"[subset-slice] direction [{', '.join(main_coords+['back', 'exit'])}]? "
                                    )
                                    subtitle = f"{slice_dir} slice"
                                else:
                                    if verbose:
                                        logger.info(
                                            f"\nA surface displays a surface plot of a 2D variable."
                                            f"\n\t{dash}\n\tback - takes you to the previous step\n\texit  "
                                        )
                                    subtitle = "surface plot"

                                # Decide on the coordinate system to work with.
                                # If in both, choose the geographic one.
                                coords = main_coords
                                if subset_type == "slice":
                                    if (
                                        slice_dir in main_coords
                                        and slice_dir in aux_coords
                                    ):
                                        if "latitude" not in coords:
                                            coords = aux_coords
                                    elif slice_dir not in coords:
                                        coords = aux_coords
                                # Done, back!
                                if slice_dir.strip() == "exit":
                                    sys.exit()
                                elif slice_dir.strip() == "back":
                                    # messages.append(f"[ERR] invalid variable {slice_dir}")
                                    slice_dir = "back"
                                    break
                                elif (
                                    slice_dir not in coords or not slice_dir.strip()
                                ) and subset_type != "surface":
                                    messages.append(
                                        f"[ERR] invalid selection {slice_dir}"
                                    )
                                    continue
                                # Explore what to do with the slice?
                                else:
                                    if subset_type == "slice":
                                        slice_value = None
                                        while slice_value is None:
                                            messages = output_messages(messages)
                                            slice_value = input(
                                                f"[subset-slice-{slice_dir}] {slice_dir} [{np.nanmin(ds[slice_dir].data)} to {np.nanmax(ds[slice_dir].data)}, back, exit]? "
                                            )

                                            # Exit.
                                            if slice_value.strip() == "exit":
                                                sys.exit()
                                            # A step back.
                                            elif (
                                                slice_value.strip() == "back"
                                                or not slice_value.strip()
                                            ):
                                                slice_value = "back"
                                                break
                                            else:
                                                try:
                                                    slice_value = float(slice_value)
                                                except Exception as ex:
                                                    messages.append(
                                                        f"[ERR] invalid value {slice_value}\n{ex}"
                                                    )
                                                    slice_value = None
                                            if slice_value == "back":
                                                break
                                        # Actions.
                                        # while slice_value:
                                        if slice_value == "back":
                                            break
                                        closest_slice_value = lib.closest(
                                            coordinate_values[slice_dir],
                                            slice_value,
                                        )
                                        messages.append(
                                            f"Slicing at the closest {slice_dir} of {closest_slice_value}"
                                        )

                                    slice_dims = coords.copy()
                                    if subset_type == "slice":
                                        slice_dims.remove(slice_dir)
                                    slice_limits = dict()
                                    slice_input_limits = dict()
                                    slice_input_order = list()
                                    gmap_limits = dict()
                                    gmap_option = ""
                                    # Geographic coordinates?
                                    if (
                                        "latitude" in coordinates
                                        and "longitude" in coordinates
                                    ):
                                        # Make sure the geographical coordinates are independent of the slice direction.
                                        if (
                                            slice_dir not in ds["latitude"].dims
                                            and slice_dir not in ds["longitude"].dims
                                        ):
                                            gmap_option = ", gmap"
                                            for dim in ["longitude", "latitude"]:
                                                gmap_limits[dim] = (
                                                    np.nanmin(ds[dim].data),
                                                    np.nanmax(ds[dim].data),
                                                )
                                    # Get the slice limits.
                                    for dim in slice_dims:
                                        slice_limits[dim] = (
                                            np.nanmin(ds[dim].data),
                                            np.nanmax(ds[dim].data),
                                        )

                                        if verbose:
                                            messages.append(
                                                f"\nSelect slice limits in the {dim} direction.\n\tlimits - provide minimum,maximum or press enter to accept the full range default \n\t{dash}\n\tback - takes you to the previous step\n\texit  "
                                            )
                                        messages = output_messages(messages)
                                        _limits = input(
                                            f"[slice-{slice_dir}] {dim} limits [default values {format_numbers(slice_limits[dim], digits=2)}, back, exit]: "
                                        )
                                        # Done, back!
                                        if _limits.strip() == "exit":
                                            sys.exit()
                                        elif _limits.strip() == "back":
                                            break
                                        elif _limits.strip() and "," not in _limits:
                                            messages.append(
                                                f"[WARN] Invalid limits, using the full range of {slice_limits[dim]}."
                                            )
                                            _limits = ""

                                        if _limits:
                                            values = _limits.split(",")
                                            try:
                                                slice_input_limits[dim] = [
                                                    float(values[0]),
                                                    float(values[1]),
                                                ]
                                                slice_input_order.append(dim)
                                                slice_limits[dim] = (
                                                    min(
                                                        float(values[0]),
                                                        float(values[1]),
                                                    ),
                                                    max(
                                                        float(values[0]),
                                                        float(values[1]),
                                                    ),
                                                )

                                                if dim in gmap_limits:
                                                    gmap_limits[dim] = (
                                                        slice_input_limits[dim]
                                                    )

                                            except Exception as ex:
                                                messages.append(
                                                    f"[ERR] Invalid limits {_limits}.\n{ex}"
                                                )
                                                _limits = "back"
                                                break
                                        else:
                                            _limits = slice_limits[dim]
                                            slice_input_order.append(dim)
                                            slice_input_limits[dim] = slice_limits[dim]
                                    # Slice the data.
                                    if _limits == "back":
                                        break
                                    if subset_type == "slice":
                                        sliced_data = slicer_lib.slicer(
                                            ds,
                                            slice_dir,
                                            closest_slice_value,
                                            slice_limits,
                                        )
                                    else:
                                        sliced_data = slicer_lib.slicer(
                                            ds,
                                            slice_dir,
                                            None,
                                            slice_limits,
                                        )
                                    messages.append(
                                        f"\nThe sliced data summary: \n{sliced_data}"
                                    )
                                    slice_action = "continue"
                                    # Actions.
                                    while slice_action != "back":

                                        if verbose:
                                            if gmap_option:
                                                logger.info(
                                                    f"\nWhat to do with the slice.\n\tplot2d - a 2D plot of {slice_dims}\n\tplot3d - a 3D plot of {slice_dims} and the model variable on the 3rd axis.\n\t\tThe plot is interactive and can be rotated.\n\tgmap - a 2D plot of {slice_dims} in geographical coordinate system\n\tChange the colormap for the plots. You can also provide cmap, vmin, and vmax, with vmin and vmax defining the minimum and maximum values for the colormap.\n\tsave - save the slice data\n\t{dash} \n\tback - takes you to the previous step\n\texit  "
                                                )
                                            else:
                                                logger.info(
                                                    f"\nWhat to do with the slice.\n\tplot2d - a 2D plot of {slice_dims}\n\tplot3d - a 3D plot of {slice_dims} and the model variable on the 3rd axis.\n\tThe plot is interactive and can be rotated.\n\tcmap - change the color map for the plots (NOTE: You may also provide cmap,vmin,vmax)\n\tsave - save the slice data\n\t{dash}\n\tback - takes you to the previous step\n\texit  "
                                                )
                                        messages = output_messages(messages)
                                        slice_action = input(
                                            f"[slice-{slice_dir}] Action [plot2d, plot3d{gmap_option}, cmap, save, back, exit]: "
                                        )
                                        # Done, back!
                                        if slice_action.strip() == "exit":
                                            sys.exit()
                                        elif (
                                            slice_action.strip() == "back"
                                            or not option.strip()
                                        ):
                                            continue
                                        # 2D plot.
                                        if slice_action == "plot2d":
                                            # Plot each model variable.
                                            plot_var = data_var[0]
                                            while plot_var:
                                                if len(data_var) > 1:
                                                    messages = output_messages(messages)
                                                    plot_var = input(
                                                        f"[slice-{slice_dir}-plot2d] Variable to plot {data_var}, back, exit]: "
                                                    )
                                                # Done, back!
                                                if plot_var.strip() == "exit":
                                                    sys.exit()
                                                elif (
                                                    plot_var.strip() == "back"
                                                    or not option.strip()
                                                ):
                                                    break
                                                elif plot_var not in data_var:
                                                    messages.append(
                                                        f"[ERR] Invalid input {plot_var}"
                                                    )
                                                    continue
                                                plot_data = sliced_data[plot_var].copy()

                                                # Set the depth axis (if exists) direction.
                                                if "depth" in plot_data.dims:
                                                    if (
                                                        "geospatial_vertical_positive"
                                                        in sliced_data.attrs
                                                    ):
                                                        if (
                                                            sliced_data.attrs[
                                                                "geospatial_vertical_positive"
                                                            ]
                                                            == "down"
                                                        ):
                                                            plot_data["depth"] = (
                                                                -plot_data["depth"]
                                                            )
                                                            slice_input_limits[
                                                                "depth"
                                                            ] = [
                                                                -1 * element
                                                                for element in slice_limits[
                                                                    "depth"
                                                                ]
                                                            ]
                                                            slice_input_limits[
                                                                "depth"
                                                            ].sort()
                                                if vmin and vmax:
                                                    plot_data.plot(
                                                        cmap=cmap, vmin=vmin, vmax=vmax
                                                    )
                                                else:
                                                    plot_data.plot(cmap=cmap)

                                                if subset_type == "slice":
                                                    plt.title(
                                                        f"{base_title}\n{subtitle} of {plot_var} at {sliced_data[slice_dir].values}  {sliced_data[slice_dir].attrs['units']}"
                                                    )
                                                else:
                                                    plt.title(
                                                        f"{base_title}\n{plot_var} Surface Plot"
                                                    )
                                                # Adjust the layout
                                                dims = list(sliced_data.dims.keys())
                                                plt.xlim(slice_input_limits[dims[1]])
                                                plt.ylim(slice_input_limits[dims[0]])
                                                plt.tight_layout()
                                                plt.show()
                                                if len(data_var) <= 1:
                                                    plot_var = "back"
                                        # 3D plot.
                                        elif slice_action == "plot3d":
                                            plot_var = data_var[0]
                                            while plot_var:
                                                if len(data_var) > 1:
                                                    messages = output_messages(messages)
                                                    plot_var = input(
                                                        f"[slice-{slice_dir}-plot3d] Variable to plot {data_var}, back, exit]: "
                                                    )
                                                # Done, back!
                                                if plot_var.strip() == "exit":
                                                    sys.exit()
                                                elif (
                                                    plot_var.strip() == "back"
                                                    or not option.strip()
                                                ):
                                                    break
                                                elif plot_var not in data_var:
                                                    messages.append(
                                                        f"[ERR] Invalid input {plot_var}"
                                                    )
                                                    continue

                                                if vmin and vmax:
                                                    sliced_data[plot_var].plot.surface(
                                                        cmap=cmap, vmin=vmin, vmax=vmax
                                                    )
                                                else:
                                                    sliced_data[plot_var].plot.surface(
                                                        cmap=cmap,
                                                    )

                                                if subset_type == "slice":
                                                    plt.title(
                                                        f"{base_title}\n{subtitle} of {plot_var} at {sliced_data[slice_dir].values} {sliced_data[slice_dir].attrs['units']}"
                                                    )
                                                else:
                                                    plt.title(
                                                        f"{base_title}\n{plot_var} Surface Plot"
                                                    )
                                                # Adjust the layout
                                                plt.tight_layout()
                                                plt.show()
                                                if len(data_var) <= 1:
                                                    plot_var = "back"
                                        # Geographic map.
                                        elif slice_action == "gmap" and gmap_option:
                                            plot_var = data_var[0]
                                            while plot_var:
                                                if len(data_var) > 1:
                                                    messages = output_messages(messages)
                                                    plot_var = input(
                                                        f"[slice-{slice_dir}-gmap] Variable to plot {data_var}, back, exit]: "
                                                    )
                                                # Done, back!
                                                if plot_var.strip() == "exit":
                                                    sys.exit()
                                                elif (
                                                    plot_var.strip() == "back"
                                                    or not option.strip()
                                                ):
                                                    plot_var = "back"
                                                    break
                                                elif plot_var not in data_var:
                                                    messages.append(
                                                        f"[ERR] Invalid input {plot_var}"
                                                    )
                                                    continue

                                                if subset_type == "slice":
                                                    title = (
                                                        f"{base_title}\n{subtitle} of {plot_var} at {sliced_data[slice_dir].values} {sliced_data[slice_dir].attrs['units']}",
                                                    )

                                                else:
                                                    title = f"{base_title}\n{plot_var} Surface Plot"

                                                slicer_lib.gmap(
                                                    plot_var,
                                                    cmap,
                                                    gmap_limits,
                                                    sliced_data,
                                                    vmin=vmin,
                                                    vmax=vmax,
                                                    title=title,
                                                )
                                                if len(data_var) <= 1:
                                                    plot_var = "back"
                                        # Colormap.
                                        elif slice_action == "cmap":
                                            messages.append(
                                                f"\n[cmaps] {plt.colormaps()}"
                                            )
                                            _cmap = prop.cmap
                                            while _cmap != "back":
                                                messages = output_messages(messages)
                                                _cmap = input(
                                                    f"[slice-{slice_dir}-cmap] Select a color map for the plot [default cmap {cmap} (NOTE: You may also provide cmap,vmin,vmax), back, exit]: "
                                                )
                                                # Done, back!
                                                if _cmap.strip() == "exit":
                                                    sys.exit()
                                                elif _cmap.strip() == "back":
                                                    _cmap = "back"
                                                    continue
                                                elif not _cmap.strip():
                                                    cmap = prop.cmap
                                                    messages.append(
                                                        f"[INFO] Using the default cmap: {cmap}"
                                                    )
                                                    break
                                                items = _cmap.split(",")
                                                if len(items) == 1:
                                                    _cmap = items[0]
                                                    if not _cmap:
                                                        _cmap = cmap
                                                    if _cmap not in plt.colormaps():
                                                        messages.append(
                                                            f"[ERR] invalid cmap: {_cmap}"
                                                        )
                                                        continue
                                                    cmap = _cmap
                                                    break
                                                elif len(items) == 3:
                                                    _cmap = items[0]
                                                    if not _cmap:
                                                        _cmap = cmap
                                                    if _cmap not in plt.colormaps():
                                                        messages.append(
                                                            f"[ERR] invalid cmap: {_cmap}"
                                                        )
                                                        continue
                                                    cmap = _cmap
                                                    vmin = float(items[1])
                                                    vmax = float(items[2])
                                                    messages.append(
                                                        f"[INFO] CMAP:{cmap}, vmin:{vmin}, vmax:{vmax}"
                                                    )
                                                    break
                                                else:
                                                    cmap = prop.cmap
                                                    messages.append(
                                                        f"[ERR] invalid cmap selection input: {','.join(items)}"
                                                    )
                                                    continue
                                        # Save the slice data.
                                        elif slice_action == "save":
                                            filename = ""
                                            while filename != "back":

                                                messages = output_messages(messages)
                                                file_info = input(
                                                    f"[slice-{slice_dir}-save] Output filename and type (filename,type)? [valid file types are: {list(valid_file_types.keys())}; back, exit]: "
                                                )

                                                # Done, back!
                                                if file_info.strip() == "exit":
                                                    sys.exit()
                                                elif file_info.strip() == "back":
                                                    break
                                                elif "," not in file_info:
                                                    messages.append(
                                                        f"[ERR] Bad input: {file_info}"
                                                    )
                                                    continue
                                                else:
                                                    try:
                                                        filename, file_type = (
                                                            file_info.split(",")
                                                        )
                                                    except Exception as ex:
                                                        messages.append(
                                                            f"[ERR] Bad input: {file_info}\n{ex}"
                                                        )
                                                        continue
                                                    if (
                                                        file_type.strip()
                                                        not in valid_file_types
                                                    ):
                                                        messages.append(
                                                            f"[ERR] Bad input: {file_info}"
                                                        )
                                                        continue
                                                    else:
                                                        filename = f"{filename.strip()}{valid_file_types[file_type]}"
                                                messages = output_data(
                                                    sliced_data, filename, messages
                                                )
                                                messages = output_messages(messages)

    else:
        logger.error(f"[ERR] {input_file} is not a valid netCDF file")
        sys.exit(1)


if __name__ == "__main__":
    main()
