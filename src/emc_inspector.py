#!/usr/bin/env python
import os
import sys
import netCDF4

description = f"""emc_inspector.py - checks a given netCDF file against the EMC requirements and flags any missing or incorrectly defined attributes. 
The check is not exhaustive, and visual validation is required to identify any additional errors.

"""
"""
USAGE:
         
       emc_inspector.py  

 HISTORY:
   2024-09-16 EarthScope DS Manoch: v2024.260 Created.
                               """

version = "v2024.260"
script = os.path.basename(sys.argv[0])
version_tag = f"""\n\n[INFO] {script} version {version}"""


def supports_unicode():
    # Check if the output supports Unicode \
    # and ensuring no exception occurs
    encoding = sys.stdout.encoding or ""
    return "UTF-8" in encoding.upper() or "UTF8" in encoding.upper()


# Status marker
if supports_unicode():
    CHECK = "✔️".encode("utf-8").decode("utf-8")
    FAIL = "❌".encode("utf-8").decode("utf-8")
    WARNING = "⚠️".encode("utf-8").decode("utf-8")
else:
    CHECK = "[OK]"
    FAIL = "[X]"
    WARNING = "[!]"


# Preferred order of dimensions.
preferred_dimensions_order = {
    "latitude": ["latitude", "longitude"],
    "longitude": ["latitude", "longitude"],
    "x": ["x", "y"],
    "y": ["x", "y"],
}

# Dictionary of required global attributes and their start string
required_global_attributes = {
    "author_institution": None,
    "author_name": None,
    "Conventions": "CF-1.0",
    "data_revision": "r",
    "id": None,
    "Metadata_Conventions": "Unidata Dataset Discovery v1.0",
    "model": None,
    "model_type": None,
    "model_subtype": None,
    "reference_pid": "doi",
    "reference": None,
    "repository_institution": None,
    "repository_name": None,
    "repository_pid": "doi:",
    "summary": None,
    "title": None,
    "year": None,
}

# Dictionary of optional global attributes and their start string
optional_global_attributes = {
    "author_email": None,
    "author_url": None,
    "grid_ref": None,
    "version": "v",
    "grid_dim": "3D",
    "grid_ref": "latitude_longitude",
}

metadata_summary = {}

# Required attributes based on variable type
required_attributes = {
    "coordinates": ["long_name", "units", "standard_name"],
    "depth": ["long_name", "units", "positive"],
    "model": ["long_name", "units", "display_name"],
}


def check_netcdf_file(file_name):
    try:
        file_size = os.path.getsize(file_name)  # Get the file size

        # Try opening the file with netCDF4
        dataset = netCDF4.Dataset(file_name, "r")

        try:
            # Determine the format
            if dataset.file_format == "NETCDF3_CLASSIC":
                metadata_summary["File Type"] = f"NetCDF-3 Classic {FAIL}"
            elif dataset.file_format == "NETCDF3_64BIT":
                metadata_summary["File Type"] = f"NetCDF-3 64-bit offset {FAIL}"
            elif dataset.file_format == "NETCDF4_CLASSIC":
                metadata_summary["File Type"] = f"NetCDF-4 Classic Model {CHECK}"
            elif dataset.file_format == "NETCDF4":
                metadata_summary["File Type"] = f"NetCDF-4 {FAIL}"
            else:
                metadata_summary["File Type"] = f"Unknown {FAIL}"
        except Exception as e:
            output(f"Error determining file format: {e} {FAIL}")
            metadata_summary["File Type"] = f"Error determining file format: {e} {FAIL}"

        print(
            f"\n\n\n{version_tag}\n\n=== NetCDF File Check ===\nFile: {file_name}\nSize: {file_size / (1024 ** 2):.2f} MB\nFile Type: {metadata_summary['File Type']}"
        )

        try:
            # List dimensions and their sizes
            output_header("Dimensions and Sizes")
            list_dimensions(dataset)

            # List primary and auxiliary coordinates
            output_header("Coordinate Information")
            primary_coords, auxiliary_coords = list_coordinates(dataset)

            # Latitude and Longitude Value Check
            output_header("Latitude and Longitude Value Check")
            lat_min, lat_max, lon_min, lon_max = check_lat_lon_values(
                dataset, primary_coords, auxiliary_coords
            )

            # Check geospatial attributes
            output_header("Geospatial Attributes Check")
            check_geospatial_attributes(dataset, lat_min, lat_max, lon_min, lon_max)

            # Check for required and optional global attributes
            output_header("Global Attributes Check")
            check_global_attributes(dataset)

            # Perform CF compliance check
            output_header("CF Compliance Check")
            check_cf_compliance(file_name)

            # List variables with their ranges and units, split into types
            output_header("Variables Summary")
            list_variables(dataset, primary_coords)

        except Exception as e:
            output(f"Error during processing: {e} {FAIL}")
            metadata_summary["Processing Error"] = (
                f"Error during processing: {e} {FAIL}"
            )

        dataset.close()

    except OSError as e:
        output(f"Error: Unable to open the file '{file_name}'. {e} {FAIL}")
        metadata_summary["File Access Error"] = (
            f"Unable to open the file '{file_name}'. {e} {FAIL}"
        )

    # Print metadata summary
    try:
        print(
            "\n\n\n============================\n",
            "Metadata Inspection Summary",
            "\n============================",
        )  # Changed to 'Metadata Inspection Summary'
        print_metadata_summary()
    except Exception as e:
        output(f"Error while printing metadata summary: {e} {FAIL}")


def list_dimensions(dataset):
    try:
        dimensions = {dim: len(dataset.dimensions[dim]) for dim in dataset.dimensions}
        metadata_summary["Dimensions"] = dimensions
        for dim, size in dimensions.items():
            output(f"- {dim}: {size} {CHECK}", indent=True)
    except Exception as e:
        output(f"Error listing dimensions: {e} {FAIL}")
        metadata_summary["Dimensions"] = f"Error listing dimensions: {e} {FAIL}"


def list_coordinates(dataset):
    try:
        primary_coords = []
        auxiliary_coords = set()  # Use a set to avoid duplicate entries
        third_dimension = None
        # Iterate over variables to identify coordinates
        for var_name, variable in dataset.variables.items():
            if hasattr(variable, "coordinates"):
                auxiliary_coords.update(variable.coordinates.split())
            elif var_name in dataset.dimensions:
                primary_coords.append(var_name)

        # Identify the third dimension
        for coord in primary_coords:
            if coord.lower() not in ["latitude", "longitude", "x", "y"]:
                third_dimension = coord
                break

        metadata_summary["Primary Coordinates"] = (
            primary_coords if primary_coords else None
        )
        metadata_summary["Auxiliary Coordinates"] = (
            list(auxiliary_coords) if auxiliary_coords else None
        )
        metadata_summary["Third Dimension"] = (
            third_dimension if third_dimension else None
        )

        if primary_coords:
            output(
                f"Primary Coordinates: {', '.join(primary_coords)} {CHECK}", indent=True
            )
            if third_dimension:
                output(f"Third Dimension: {third_dimension} {CHECK}", indent=True)
        else:
            output(f"No primary coordinates found. {FAIL}", indent=True)
            metadata_summary[f"Primary Coordinates"] = (
                "No primary coordinates found {FAIL}"
            )

        if auxiliary_coords:
            output(
                f"Auxiliary Coordinates: {', '.join(auxiliary_coords)} {CHECK}",
                indent=True,
            )
        else:
            output(f"No auxiliary coordinates found. {WARNING} (Optional)", indent=True)

        # Return both primary and auxiliary coordinates
        return primary_coords, auxiliary_coords

    except Exception as e:
        output(f"Error listing coordinates: {e} {FAIL}")
        metadata_summary["Coordinates"] = f"Error listing coordinates: {e} {FAIL}"
        return [], set()  # Return empty lists in case of an error


def check_lat_lon_values(dataset, primary_coords, auxiliary_coords):
    try:
        lat_name = None
        lon_name = None

        # Identify latitude and longitude variable names
        for coord in primary_coords + list(auxiliary_coords):
            if "lat" in coord.lower():
                lat_name = coord
            if "lon" in coord.lower():
                lon_name = coord

        if lat_name and lon_name:
            lat_values = dataset.variables[lat_name][:]
            lon_values = dataset.variables[lon_name][:]

            # Determine longitude format (0-360 or -180/180)
            if lon_values.min() >= 0 and lon_values.max() <= 180:
                lon_format = "Potentially 0-360 degrees or -180/180 degrees"
            elif lon_values.min() >= 0 and lon_values.max() <= 360:
                lon_format = "0-360 degrees"
            elif lon_values.min() >= -180 and lon_values.max() <= 180:
                lon_format = "-180/180 degrees"
            else:
                lon_format = "Unknown"

            metadata_summary["Longitude Format"] = lon_format
            output(f"Longitude Format: {lon_format} {CHECK}", indent=True)

            # Check that latitudes and longitudes are within valid ranges
            lat_min, lat_max = lat_values.min(), lat_values.max()
            lon_min, lon_max = lon_values.min(), lon_values.max()

            if lat_min < -90 or lat_max > 90:
                output(
                    f"Latitude values are out of range: {lat_min:.3f} to {lat_max:.3f} {FAIL}",
                    indent=True,
                )
                metadata_summary["Latitude Range"] = (
                    f"Latitude values out of range: {lat_min:.3f} to {lat_max:.3f} {FAIL}"
                )
            else:
                output(
                    f"Latitude Range: {lat_min:.3f} to {lat_max:.3f} {CHECK}",
                    indent=True,
                )
                metadata_summary["Latitude Range"] = f"{lat_min:.3f} to {lat_max:.3f}"

            if lon_format == "0-360 degrees" and (lon_min < 0 or lon_max > 360):
                output(
                    f"Longitude values are out of range: {lon_min:.3f} to {lon_max:.3f} {FAIL}",
                    indent=True,
                )
                metadata_summary["Longitude Range"] = (
                    f"Longitude values out of range: {lon_min:.3f} to {lon_max:.3f} {FAIL}"
                )
            elif lon_format == "-180/180 degrees" and (lon_min < -180 or lon_max > 180):
                output(
                    f"Longitude values are out of range: {lon_min:.3f} to {lon_max:.3f} {FAIL}",
                    indent=True,
                )
                metadata_summary["Longitude Range"] = (
                    f"Longitude values out of range: {lon_min:.3f} to {lon_max:.3f} {FAIL}"
                )
            elif lon_format == "Potentially 0-360 degrees or -180/180 degrees":
                output(
                    f"Longitude Range: {lon_min:.3f} to {lon_max:.3f} {CHECK}",
                    indent=True,
                )
                metadata_summary["Longitude Range"] = f"{lon_min:.3f} to {lon_max:.3f}"
            else:
                output(
                    f"Longitude Range: {lon_min:.3f} to {lon_max:.3f} {CHECK}",
                    indent=True,
                )
                metadata_summary["Longitude Range"] = f"{lon_min:.3f} to {lon_max:.3f}"

            return lat_min, lat_max, lon_min, lon_max  # Return computed values

        else:
            output(f"Latitude and/or Longitude not found. {FAIL}", indent=True)
            metadata_summary["Coordinates"] = (
                f"Latitude and/or Longitude not found {FAIL}"
            )
            return None, None, None, None

    except Exception as e:
        output(f"Error checking latitude and longitude values: {e} {FAIL}")
        metadata_summary["Coordinates"] = (
            f"Error checking latitude and longitude values: {e} {FAIL}"
        )
        return None, None, None, None


def check_geospatial_attributes(dataset, lat_min, lat_max, lon_min, lon_max):
    try:
        required_geospatial_attrs = [
            "geospatial_lat_min",
            "geospatial_lat_max",
            "geospatial_lon_min",
            "geospatial_lon_max",
            "geospatial_vertical_min",
            "geospatial_vertical_max",
        ]

        missing_attrs = [
            attr for attr in required_geospatial_attrs if attr not in dataset.ncattrs()
        ]

        if missing_attrs:
            for attr in missing_attrs:
                if (
                    "vertical" in attr and "depth" in dataset.dimensions
                ):  # Check if depth is present
                    output(
                        f"Error: {attr} is missing but depth dimension is present. {FAIL}",
                        indent=True,
                    )
                    metadata_summary[attr] = (
                        f"{attr} is missing but depth dimension is present. {FAIL}"
                    )
                elif "vertical" in attr:
                    output(
                        f"Warning: {attr} is missing (model does not have a depth dimension). {WARNING}",
                        indent=True,
                    )
                else:
                    output(
                        f"Error: Required attribute {attr} is missing. {FAIL}",
                        indent=True,
                    )
                    metadata_summary[attr] = (
                        f"Required attribute {attr} is missing. {FAIL}"
                    )
            return  # Skip further checks if required attributes are missing

        # Check if depth is a dimension and handle accordingly
        if "depth" in dataset.dimensions:
            depth_values = dataset.variables["depth"][:]
            depth_min = depth_values.min()
            depth_max = depth_values.max()
        else:
            depth_min = None
            depth_max = None

        # Use computed latitude and longitude values from earlier
        geospatial_metadata = {
            "geospatial_lat_min": lat_min,
            "geospatial_lat_max": lat_max,
            "geospatial_lon_min": lon_min,
            "geospatial_lon_max": lon_max,
            "geospatial_vertical_min": (
                depth_min if depth_min is not None else "Not applicable"
            ),
            "geospatial_vertical_max": (
                depth_max if depth_max is not None else "Not applicable"
            ),
        }

        # Compare with geospatial metadata
        for key, value in geospatial_metadata.items():
            if key in dataset.ncattrs():
                metadata_value = dataset.getncattr(key)
                if value == "Not applicable":
                    output(
                        f"{key}: {value} (No depth dimension) {WARNING}", indent=True
                    )
                elif type(metadata_value) == str:
                    output(
                        f'{key}: "{metadata_value}", Cannot be a string; it must match the variable\'s type. {FAIL}',
                        indent=True,
                    )
                    metadata_summary[key] = (
                        f'{key} = "{metadata_value}", is a string {FAIL}'
                    )
                elif is_within_tolerance(metadata_value, value):
                    output(
                        f"{key}: Metadata = {metadata_value}, Computed = {value:.3f} {CHECK}",
                        indent=True,
                    )
                else:
                    output(
                        f"{key}: Metadata = {metadata_value}, Computed = {value:.3f} {FAIL}",
                        indent=True,
                    )
                    output(
                        f"Error: {key} mismatch. Metadata = {metadata_value}, Computed = {value:.3f}",
                        indent=True,
                        level=2,
                    )
                    metadata_summary[key] = (
                        f"Metadata = {metadata_value}, Computed = {value:.3f} {FAIL}"
                    )
            else:
                output(
                    f"{key} is not found in metadata. {WARNING} (Optional)", indent=True
                )
                if value is not None and value != "Not applicable":
                    output(f"Computed {key}: {value:.3f} {CHECK}", indent=True)
    except Exception as e:
        output(f"Error checking geospatial attributes: {e} {FAIL}")
        metadata_summary["Geospatial Attributes"] = (
            f"Error checking geospatial attributes: {e} {FAIL}"
        )


def is_within_tolerance(metadata_value, computed_value, tolerance=0.1):
    try:
        if metadata_value is None or computed_value is None:
            return False
        return abs(metadata_value - computed_value) <= abs(tolerance * metadata_value)
    except Exception as e:
        output(f"Error during tolerance check: {e} {FAIL}")
        return False


def check_global_attributes(dataset):
    try:
        missing_required_attributes = []
        empty_required_attributes = []
        checked_required_attributes = []
        missing_optional_attributes = []
        empty_optional_attributes = []
        checked_optional_attributes = []

        # Required
        for attribute in required_global_attributes:
            if attribute not in dataset.ncattrs():
                missing_required_attributes.append(attribute)
            elif not str(dataset.getncattr(attribute)).strip():
                empty_required_attributes.append(attribute)
            else:
                checked_required_attributes.append(attribute)

        if missing_required_attributes:
            output(f"Missing required global attributes: {FAIL}", indent=True)
            for attr in missing_required_attributes:
                output(f"- {attr}", indent=True, level=2)
            metadata_summary["Missing required global attributes"] = (
                missing_required_attributes
            )

        if empty_required_attributes:
            output(f"Required global attributes NOT assigned: {FAIL}", indent=True)
            for attr in empty_required_attributes:
                output(f"- {attr}", indent=True, level=2)
            metadata_summary["Required global attributes present but not assigned"] = (
                empty_required_attributes
            )

        output(f"Global attributes (values not checked): {CHECK}", indent=True)
        for attr in checked_required_attributes:
            output(f"- {attr}", indent=True, level=2)
        metadata_summary["Present global attributes (values not checked)"] = (
            checked_required_attributes
        )

        # Optional
        for attribute in optional_global_attributes:
            if attribute not in dataset.ncattrs():
                missing_optional_attributes.append(attribute)
            elif not str(dataset.getncattr(attribute)).strip():
                empty_optional_attributes.append(attribute)
            else:
                checked_optional_attributes.append(attribute)

        if missing_optional_attributes:
            output(f"Missing optional global attribute(s): {WARNING}", indent=True)
            for attr in missing_optional_attributes:
                output(f"- {attr}", indent=True, level=2)
            metadata_summary["Missing optional global attributes"] = (
                missing_optional_attributes
            )

        if empty_optional_attributes:
            output(f"Optional global attributes NOT assigned: {WARNING}", indent=True)
            for attr in empty_optional_attributes:
                output(f"- {attr}", indent=True, level=2)
            metadata_summary["Optional global attributes present but not assigned"] = (
                empty_optional_attributes
            )

        output(f"Optional global attributes (values not checked): {CHECK}", indent=True)
        for attr in checked_optional_attributes:
            output(f"- {attr}", indent=True, level=2)
        metadata_summary["Optional global attributes (values not checked)"] = (
            checked_optional_attributes
        )
    except Exception as e:
        output(f"Error checking global attributes: {e} {FAIL}")
        metadata_summary["Global Attributes"] = (
            f"Error checking global attributes: {e} {FAIL}"
        )


def check_cf_compliance(file_name):
    try:
        # Use the cfchecker command-line tool to check CF compliance
        result = os.popen(f"cfchecks -v CF-1.0 {file_name}").read()

        # Initialize status variables
        errors_detected = 0
        warnings_given = 0
        information_messages = 0

        # Parse the output for specific lines indicating status
        for line in result.splitlines():
            if "ERRORS detected" in line:
                errors_detected = int(line.split(":")[-1].strip())
            elif "WARNINGS given" in line:
                warnings_given = int(line.split(":")[-1].strip())
            elif "INFORMATION messages" in line:
                information_messages = int(line.split(":")[-1].strip())

        # Determine the compliance result based on the parsed information
        if errors_detected == 0 and warnings_given == 0:
            output(f"CF Compliance: Compliant {CHECK}", indent=True)
            metadata_summary["CF Compliance"] = "Compliant"
        elif errors_detected == 0 and warnings_given > 0:
            output(f"CF Compliance: with warning {WARNING}", indent=True)
            output(f"Warnings: {warnings_given} {WARNING}", indent=True, level=2)
            metadata_summary["CF Compliance"] = f"with warning {WARNING}"
        elif errors_detected > 0:
            output(f"CF Compliance: Non-compliant {FAIL}", indent=True)
            output(f"Issues:\n{result}", indent=True, level=2)
            metadata_summary[f"CF Compliance"] = "Non-compliant {FAIL}"
    except Exception as e:
        output(f"CF compliance check failed for '{file_name}'. Error: {e} {FAIL}")
        metadata_summary["CF Compliance"] = f"Check failed: {e} {FAIL}"


def list_variables(dataset, primary_coords):
    try:
        coordinate_variables = {}
        auxiliary_coordinates = set()  # Set to store auxiliary coordinates
        model_variables = {}

        # Retrieve the third dimension from the metadata summary
        third_dimension = metadata_summary.get("Third Dimension")

        # Step 1: Identify all auxiliary coordinates from model variable attributes
        for var_name, variable in dataset.variables.items():
            if "coordinates" in variable.ncattrs():
                coordinates = variable.getncattr("coordinates").split()
                auxiliary_coordinates.update(coordinates)

        for var_name, variable in dataset.variables.items():
            var_type = None
            if third_dimension and "depth" in var_name.lower():
                var_type = "depth"
            elif var_name in dataset.dimensions:
                var_type = "coordinates"
            elif var_name in auxiliary_coordinates:
                var_type = "auxiliary"

            else:
                var_type = "model"

            # Check for required attributes
            missing_attributes = []
            if var_type == "model":
                required_attributes_model = ["long_name", "units", "display_name"]
                for attr in required_attributes_model:
                    if attr not in variable.ncattrs():
                        missing_attributes.append(attr)
            else:
                required_attributes_coordinates = [
                    "long_name",
                    "units",
                    "standard_name",
                ]
                if var_type in ["depth", "elevation"]:
                    required_attributes_coordinates = [
                        "long_name",
                        "units",
                        "positive",
                    ]
                for attr in required_attributes_coordinates:
                    if attr not in variable.ncattrs():
                        missing_attributes.append(attr)

            var_min = variable[:].min()
            var_max = variable[:].max()
            var_units = variable.units if "units" in variable.ncattrs() else "No units"
            var_dims = list(variable.dimensions)
            var_summary = {
                "Range": f"{var_min:.3f} to {var_max:.3f}",
                "Units": var_units,
                "Dimensions": ", ".join(var_dims),
            }
            if missing_attributes:
                var_summary["Missing Attributes"] = missing_attributes

            # Determine expected order of dimensions based on coordinates and third dimension
            if var_type == "model" and len(variable.dimensions) == 3:
                # Ensure the third dimension is not already in primary_coords
                for coord in primary_coords:
                    if coord != third_dimension:
                        expected_order = preferred_dimensions_order[coord]
                        break

                expected_order = [third_dimension] + expected_order

                if var_dims != expected_order:
                    output(
                        f"Error: Variable '{var_name}' has incorrect dimension order: {', '.join(var_dims)}. Expected order is '{', '.join(expected_order)}'. {FAIL}",
                        indent=True,
                    )
                    metadata_summary[var_name] = (
                        f"Incorrect dimension order: {', '.join(var_dims)} {FAIL}"
                    )
                else:
                    output(
                        f"Variable '{var_name}' has correct dimension order: {', '.join(var_dims)}. {CHECK}",
                        indent=True,
                    )

            # Classify variables
            if var_type == "coordinates" or (var_type == "depth" and third_dimension):
                coordinate_variables[var_name] = var_summary
            elif var_type == "auxiliary" and var_type != "depth":
                coordinate_variables[var_name] = (
                    var_summary  # Auxiliary coordinates are also part of coordinates
                )
            else:
                model_variables[var_name] = var_summary

        output("Coordinate Variables:", indent=True)
        for var_name, summary in coordinate_variables.items():
            output(
                f"- {var_name}: Range = {summary['Range']}, Units = {summary['Units']}, Dimensions = {summary['Dimensions']} {CHECK}",
                indent=True,
                level=2,
            )
            if "Missing Attributes" in summary:
                output(
                    f"Missing Attributes: {', '.join(summary['Missing Attributes'])} {FAIL}",
                    indent=True,
                    level=3,
                )
                metadata_summary[var_name] = (
                    f"Missing Attributes: {', '.join(summary['Missing Attributes'])} {FAIL}"
                )

        output("Model Variables:", indent=True)
        for var_name, summary in model_variables.items():
            output(
                f"- {var_name}: Range = {summary['Range']}, Units = {summary['Units']}, Dimensions = {summary['Dimensions']} {CHECK}",
                indent=True,
                level=2,
            )
            if "Missing Attributes" in summary:
                output(
                    f"Missing Attributes: {', '.join(summary['Missing Attributes'])} {FAIL}",
                    indent=True,
                    level=3,
                )
                metadata_summary[var_name] = (
                    f"Missing Attributes: {', '.join(summary['Missing Attributes'])} {FAIL}"
                )

        metadata_summary["Coordinate Variables"] = coordinate_variables
        metadata_summary["Model Variables"] = model_variables
    except Exception as e:
        output(f"Error while listing variables: {e} {FAIL}")
        metadata_summary["Variables"] = f"Error while listing variables: {e} {FAIL}"


def print_metadata_summary():
    try:
        for key, value in metadata_summary.items():
            if value is not None:
                if isinstance(value, dict):
                    output(f"{key}:", indent=True)
                    for sub_key, sub_value in value.items():
                        if isinstance(sub_value, str) and FAIL in sub_value:
                            output(
                                f"- {sub_key}: {sub_value} {FAIL}", indent=True, level=2
                            )
                        else:
                            output(f"- {sub_key}: {sub_value}", indent=True, level=2)
                else:
                    if isinstance(value, str) and FAIL in value:
                        output(f"{key}: {value} {FAIL}", indent=True)
                    else:
                        output(f"{key}: {value}", indent=True)
    except Exception as e:
        output(f"Error printing metadata summary: {e} {FAIL}")


def output(text, indent=False, level=1):
    indentation = "  " * (level - 1) if indent else ""
    print(f"{indentation}{text}")


def output_header(header):
    print(f"\n=== {header} ===\n")


if __name__ == "__main__":
    # Prompt the user to enter the file name
    file_name = input(f"\n\n{description}\nPlease enter the netCDF file name: ").strip()
    check_netcdf_file(file_name)
