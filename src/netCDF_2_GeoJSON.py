import os
import sys
import xarray as xr
import json
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
import getopt
import logging
from pathlib import Path

# Configure basic logging
logging.basicConfig(level=logging.DEBUG)  # DEBUG, INFO, WARNING, ERROR, CRITICAL

NETCDF_FILE_NAME = None
SCRIPT = sys.argv[0]
this_app = os.path.splitext(SCRIPT)[0]
logger = logging.getLogger(this_app)

"""Create geojson file from a given netcdf file's metadata"""


def usage():
    print(f"\n{SCRIPT} Create geojson file from a given netcdf file's metadata.\n\n")
    print(
        f" USAGE:\n\n",
        f"{SCRIPT} -i FILE\n"
        f"\t-i [required] FILE is the input netCDF Earth model file\n",
    )


def to_serializable(val):
    if isinstance(val, (np.integer, np.floating)):
        return val.item()
    elif isinstance(val, (np.ndarray, list, tuple)):
        return [to_serializable(v) for v in val]
    else:
        return str(val) if not isinstance(val, (int, float, str)) else val


# Check if a filename is provided on the command line
try:
    options, remainder = getopt.getopt(
        sys.argv[1:],
        "hi:",
        [
            "help",
            "input=",
        ],
    )
    for opt, arg in options:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            NETCDF_FILE_NAME = arg
except getopt.GetoptError:
    usage()
    sys.exit(2)

if not NETCDF_FILE_NAME:
    usage()
    logger.error("Missing the input file:")
    sys.exit(1)

try:
    # Detect NetCDF file format using netCDF4.Dataset
    try:
        with Dataset(NETCDF_FILE_NAME, "r") as nc:
            file_format = nc.file_format
        nc_file = Path(NETCDF_FILE_NAME)  # now a Path object
    except Exception:
        file_format = "unknown"

    ds = xr.open_dataset(NETCDF_FILE_NAME)
    attrs = ds.attrs
    props = {
        "model": to_serializable(attrs.get("model", "")),
        "title": to_serializable(attrs.get("title", "")),
        "summary": to_serializable(attrs.get("summary", "")),
        "reference": to_serializable(attrs.get("reference", "")),
        "reference_pid": to_serializable(
            attrs.get("reference_pid", "").replace(" ", "")
        ),
        "repository_doi": to_serializable(
            attrs.get("repository_pid", "").replace(" ", "")
        ),
        "year": to_serializable(attrs.get("year", "")),
        "data_revision": to_serializable(attrs.get("data_revision", "r0.0")),
        "version": to_serializable(attrs.get("version", "v0.0")),
        "model_type": to_serializable(attrs.get("model_type", "")),
        "model_subtype": to_serializable(attrs.get("model_subtype", "")),
        "filename": to_serializable(nc_file.name),
        "size_mb": round(nc_file.stat().st_size / (1024 * 1024), 2),
        "format": file_format,
        "grid_ref": to_serializable(attrs.get("grid_ref", "latitude_longitude")),
        "projection_type": to_serializable(
            attrs.get("emc:projection_type", "geographic_wgs84")
        ),
        "projection_notes": to_serializable(
            attrs.get(
                "emc:projection_notes",
                "Unprojected geographic coordinates, WGS 84 datum",
            )
        ),
        "grid_dim": to_serializable(attrs.get("grid_dim", "3D")),
        "geospatial_lat_min": to_serializable(
            float(attrs.get("geospatial_lat_min", "0"))
        ),
        "geospatial_lat_max": to_serializable(
            float(attrs.get("geospatial_lat_max", "0"))
        ),
        "geospatial_lat_resolution": to_serializable(
            float(attrs.get("geospatial_lat_resolution", "0"))
        ),
        "geospatial_lat_units": to_serializable(attrs.get("geospatial_lat_units", "")),
        "geospatial_lon_min": to_serializable(
            float(attrs.get("geospatial_lon_min", "0"))
        ),
        "geospatial_lon_max": to_serializable(
            float(attrs.get("geospatial_lon_max", "0"))
        ),
        "geospatial_lon_resolution": to_serializable(
            float(attrs.get("geospatial_lon_resolution", "0"))
        ),
        "geospatial_lon_units": to_serializable(attrs.get("geospatial_lon_units", "")),
        "geospatial_vertical_min": to_serializable(
            float(attrs.get("geospatial_vertical_min", "0"))
        ),
        "geospatial_vertical_max": to_serializable(
            float(attrs.get("geospatial_vertical_max", "0"))
        ),
        "geospatial_vertical_units": to_serializable(
            attrs.get("geospatial_vertical_units", "")
        ),
        "variables": {
            to_serializable(vname): to_serializable(var.attrs.get("long_name", ""))
            for vname, var in ds.data_vars.items()
        },
    }

    feature = {
        "type": "Feature",
        "geometry": {
            "type": "Polygon",
            "coordinates": [
                [
                    [props["geospatial_lon_min"], props["geospatial_lat_min"]],
                    [props["geospatial_lon_max"], props["geospatial_lat_min"]],
                    [props["geospatial_lon_max"], props["geospatial_lat_max"]],
                    [props["geospatial_lon_min"], props["geospatial_lat_max"]],
                    [props["geospatial_lon_min"], props["geospatial_lat_min"]],
                ]
            ],
        },
        "properties": props,
    }

    # Save GeoJSON file per model
    output_path = f"{nc_file.stem}.geojson"
    with open(output_path, "w") as f:
        json.dump({"type": "FeatureCollection", "features": [feature]}, f, indent=2)

    logger.info(f"Saved {output_path}")

except Exception as e:
    logger.error(f"Failed to process {nc_file.name}: {e}")
