import os
import sys
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import numpy as np
import xarray as xr

import warnings

warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt

# Get the root _directory.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop

LIB_DIR = os.path.join(ROOT_DIR, prop.LIB_DIR)
sys.path.append(LIB_DIR)
import shared_lib as lib

# Set up the logger.
logger = lib.get_logger()


def get_point(location, vars):
    f"""
    Get the coordinates for a point.

    Call arguments:
        location: point's location.
    """
    point = None
    while point is None:
        if location in ("start", "end"):
            point = input(
                f"[slice-xsection] Cross-section {location} point as {vars[0]},{vars[1]} ['back', 'exit']? "
            )
        else:
            point = input(
                f"[slice-xsection] Cross-section {vars[2]} range as start, end ['back', 'exit']? "
            )

        # Done, back!
        if point.strip() == "exit":
            sys.exit()
        elif point.strip() == "back" or not point.strip():
            return "back"
        elif "," not in point:
            if location in ("start", "end"):
                logger.error(
                    f"[ERR] invalid {location} coordinates '{point}' input {location} as {vars[0]},{vars[1]}"
                )
            else:
                logger.error(
                    f"[ERR] invalid {location} range '{point}' input {location} as start, end {vars[2]}"
                )
            point = None
        else:
            try:
                point = point.split(",")
                if len(point) != 2:
                    logger.error(
                        f"[ERR] invalid {location} range '{point}' input {location} as start, end {vars[2]}"
                    )
                    point = None
                point = [float(i) for i in point]
                if location in ("start", "end"):
                    if lib.is_in_range(point[0], vars[0]) and lib.is_in_range(
                        point[1], vars[1]
                    ):
                        break
                    logger.error(
                        f"[ERR] invalid {vars[0]},{vars[1]} pair: {point[0]},{point[1]}"
                    )
                    point = None
                # Depth.
                else:
                    break
            except Exception as ex:
                logger.error(
                    f"[ERR] invalid {location} range '{point}' input {location} as value1,value2\n{ex}"
                )
            point = None
    return point


def display_var_meta(ds, meta, indent, values=True):
    """
    Display metadata fora given variable.

    Call arguments:
        ds - [required] the xarray dataset
        meta - [required] dataset metadata
        indent - [required] display indent
        values - [optional, default: True] display the values
    """
    for var in meta:
        for attr_indx, attr in enumerate(ds[var].attrs):
            if attr_indx == 0:
                logger.info(f"\t{var}: ")
            if attr == "variable":
                ConnectionRefusedError
            else:
                logger.info(f"{indent*' '}{attr}: {ds[var].attrs[attr]}")
        if values:
            logger.info(f"{indent*' '}Values:")
            lib.view_list(ds[var].data, indent=indent)


def slicer(ds, slice_dir, slice_value, slice_limits):
    """
    Slice a dataset along a given variable.

    Call arguments:
        ds - [required] the xarray dataset
        slice_dir - [required] the variable name along which to slice (None for surface plot)
        slice_value - [required] the slice location
        slice_limits - [required] limits of the slice in other directions.
    """
    slice_limits_keys = list(slice_limits.keys())
    slice_limits_values = list(slice_limits.values())

    sliced_data = ds.where(
        (ds[slice_limits_keys[0]] >= slice_limits_values[0][0])
        & (ds[slice_limits_keys[0]] <= slice_limits_values[0][1])
        & (ds[slice_limits_keys[1]] >= slice_limits_values[1][0])
        & (ds[slice_limits_keys[1]] <= slice_limits_values[1][1]),
        drop=True,
    )
    if slice_value is not None:
        sliced_data = sliced_data.sel({slice_dir: float(slice_value)})
    return sliced_data


def subsetter(ds, limits, ds_type):
    """
    Subset a dataset as a volume.

    Call arguments:
        ds - [required] the xarray dataset
        limits - [required] limits of the volume in all directions.
        ds_type - [required] dataset type 2D or 3D
    """
    geospatial_dict = {
        "latitude": ["geospatial_lat_min", "geospatial_lat_max"],
        "longitude": ["geospatial_lon_min", "geospatial_lon_max"],
    }
    if ds_type == "3D":
        geospatial_dict["depth"] = (
            ["geospatial_vertical_min", "geospatial_vertical_max"],
        )

    # Check if the array has any zero-sized dimensions
    warnings = ""
    try:
        limit_keys = list(limits.keys())
        limit_values = list(limits.values())
        if ds_type == "3D":
            sliced_data = ds.where(
                (ds[limit_keys[0]] >= limit_values[0][0])
                & (ds[limit_keys[0]] <= limit_values[0][1])
                & (ds[limit_keys[1]] >= limit_values[1][0])
                & (ds[limit_keys[1]] <= limit_values[1][1])
                & (ds[limit_keys[2]] >= limit_values[2][0])
                & (ds[limit_keys[2]] <= limit_values[2][1]),
                drop=True,
            )
        else:
            sliced_data = ds.where(
                (ds[limit_keys[0]] >= limit_values[0][0])
                & (ds[limit_keys[0]] <= limit_values[0][1])
                & (ds[limit_keys[1]] >= limit_values[1][0])
                & (ds[limit_keys[1]] <= limit_values[1][1]),
                drop=True,
            )

        for dim in limit_keys:
            if dim in geospatial_dict:
                #  The dropna method is used to remove coordinates with all NaN values along the specified dimensions
                sliced_data = sliced_data.dropna(dim=dim, how="all")
                if geospatial_dict[dim][0] in sliced_data.attrs:
                    sliced_data.attrs[geospatial_dict[dim][0]] = min(
                        sliced_data[dim].values
                    )
                if geospatial_dict[dim][1] in sliced_data.attrs:
                    sliced_data.attrs[geospatial_dict[dim][1]] = max(
                        sliced_data[dim].values
                    )
    except Exception as ex:
        warnings = ex
        return ds, warnings

    return sliced_data, warnings


def gmap(plot_var, cmap, gmap_limits, sliced_data, vmin=None, vmax=None, title=None):
    """
    Geographical display of a dataset.

    Call arguments:
        data_var - [required] the variable for which to create the plot.
        cmap - [required] colormap to use.
        gmap_limits - [required] surface geographical limits.
        slice_limits - [required] sliced data to plot.
    """
    # Color from cmap.
    color = cmap
    # Defining the figure
    fig = plt.figure(facecolor="w", edgecolor="k")
    # Axes with Cartopy projection
    ax = plt.axes(projection=ccrs.PlateCarree())
    # and extent
    ax.set_extent(
        [
            gmap_limits["longitude"][0],
            gmap_limits["longitude"][1],
            gmap_limits["latitude"][0],
            gmap_limits["latitude"][1],
        ],
        ccrs.PlateCarree(),
    )
    # Plotting using Matplotlib
    try:
        _x = sliced_data[plot_var]["longitude"]
        nan_count = np.isnan(_x).sum().item()
        message = f"Number of NaN values in  longitude: {nan_count}"
        _y = sliced_data[plot_var]["latitude"]
        nan_count = np.isnan(_x).sum().item()
        message = f"{message}\nNumber of NaN values in longitude: {nan_count}"
        if vmin and vmax:
            cf = sliced_data[plot_var].plot(
                x="longitude",
                y="latitude",
                transform=ccrs.PlateCarree(),
                cmap=color,
                add_colorbar=True,
                vmin=vmin,
                vmax=vmax,
            )
        else:
            cf = sliced_data[plot_var].plot(
                x="longitude",
                y="latitude",
                transform=ccrs.PlateCarree(),
                cmap=color,
                add_colorbar=True,
            )

        # Plot lat/lon grid
        gl = ax.gridlines(
            crs=ccrs.PlateCarree(),
            draw_labels=True,
            linewidth=0.1,
            color="k",
            alpha=1,
            linestyle="--",
        )
        gl.top_labels = False
        gl.right_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        # Add map features with Cartopy
        # Set environment variables to suppress GDAL/PROJ debug output
        os.environ["CPL_DEBUG"] = "OFF"
        os.environ["PROJ_DEBUG"] = "OFF"
        ax.add_feature(
            cfeature.NaturalEarthFeature(
                "physical",
                "land",
                "10m",
                edgecolor="face",
                facecolor="none",
            )
        )
        ax.coastlines(linewidth=1)
        if title is not None:
            plt.title(title)
        plt.tight_layout()
        plt.show()
    except Exception as ex:
        logger.error(
            f"Failed to create the plot. \n{ex}\nLongitude:{_x}\nLatitude: {_y}\n{message}"
        )


def interpolate_path(
    ds,
    start,
    end,
    num_points=100,
    method="linear",
    grid_ref="latitude_longitude",
    utm_zone=None,
    ellipsoid=None,
):
    """
    Interpolates a dataset along a path defined by start and end coordinates on an irregular grid.

    Parameters:
        ds (xarray.Dataset): The input dataset containing 'latitude' and 'longitude' as coordinates.
        start (tuple): A tuple (latitude, longitude) of the starting point.
        end (tuple): A tuple (latitude, longitude) of the ending point.
        num_points (int): Number of points to interpolate along the path.
        method (str): Interpolation method to use ('linear', 'nearest').

    Returns:
        xarray.Dataset: The interpolated dataset along the path.
    """
    # Create linearly spaced points between start and end
    lat_points = np.linspace(start[1], end[1], num_points)
    lon_points = np.linspace(start[0], end[0], num_points)

    # Define a path dataset for interpolation
    path = xr.Dataset(
        {"latitude": ("points", lat_points), "longitude": ("points", lon_points)}
    )

    # Interpolate the dataset to these points using the specified method
    if grid_ref == "latitude_longitude":
        interpolated_ds = ds.interp(
            latitude=path.latitude, longitude=path.longitude, method=method
        )
    else:
        if None in (utm_zone, ellipsoid):
            message = f"[ERR] for grid_ref: {grid_ref}, utm_zone and ellipsoid are required. Current values: {utm_zone}, {ellipsoid}!"
            logger.error(message)
            raise
        x_points = list()
        y_points = list()
        for index, lat_value in enumerate(lat_points):
            x, y = lib.project_lonlat_utm(
                lon_points[index], lat_points[index], utm_zone, ellipsoid=ellipsoid
            )
            x_points.append(x)
            y_points.append(y)

        # Define a path dataset for interpolation
        path = xr.Dataset({"x": ("points", x_points), "y": ("points", y_points)})
        interpolated_ds = ds.interp(x=path.x, y=path.y, method=method)

    return interpolated_ds, lat_points, lon_points
