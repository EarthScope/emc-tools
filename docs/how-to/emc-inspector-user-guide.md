# EMC Inspector User Guide

------------------------------------------------------------------------

## Introduction

`emc_inspector.py` validates a **NetCDF** file against EarthScope's EMC
conventions and common CF metadata expectations. It prints clear,
emoji-marked results (✔️ / ❌ / ⚠️), plus a compact **Metadata
Inspection Summary** at the end.

## 1. Quick Start

``` bash
# From the EMC_Tools root (recommended)
python src/emc_inspector.py

# or, if you are inside that folder already
python emc_inspector.py
```

When prompted:

    Please enter the NetCDF file name: /path/to/model.nc

The script will:

1.  Open the file and detect the NetCDF format.
2.  Report **dimensions**, **coordinate variables** (primary &
    auxiliary), and identify a third dimension if present.
3.  Check **latitude/longitude value ranges** and report the longitude
    system (0--360 or --180/180).
4.  Validate **geospatial global attributes** against computed values.
5.  Check **required/optional global attributes**.
6.  Run a **CF compliance** check via `cfchecks` (CF-1.0) if available.
7.  Summarize everything in **Metadata Inspection Summary**.

------------------------------------------------------------------------

## 2. What the Script Checks

### - File Type

Reports the underlying NetCDF container: - **NETCDF-4 Classic Model**:
✔️ (Preferred) - **NETCDF-4**: ❌ (not preferred for EMC) - **NETCDF-3
Classic / 64-bit offset**: ❌ (not preferred) - Unknown/Errors: ❌

### - Dimensions

Lists each dimension and its size, e.g.:

    === Dimensions and Sizes ===
    - longitude: 721 ✔️
    - latitude: 361 ✔️
    - depth: 100 ✔️

### - Coordinate Information

Identifies **primary** coordinates (dimensions) and **auxiliary**
coordinates (via the `coordinates` attribute found on variables).\
It also tries to recognize a "third dimension" (e.g., `depth`).

### - Latitude/Longitude Range Check

-   Detects whether longitude values look like **0--360** or
    **--180/180**.
-   Validates ranges:
    -   Latitude must be within **\[-90, 90\]**.
    -   Longitude must match the detected system.
-   Prints min/max for lat/lon and flags out-of-range values.

### - Geospatial Global Attributes

Confirms presence/values of:

    geospatial_lat_min, geospatial_lat_max,
    geospatial_lon_min, geospatial_lon_max,
    geospatial_vertical_min, geospatial_vertical_max

**Note**: Values assigned to **geospatial\_**\* attributes must match
the data type of their corresponding coordinate variables.

-   If `depth` exists as a dimension, vertical attributes are
    **required**.
-   If values exist but are strings, that's flagged (they must be
    numeric, see the **Note** above).
-   Metadata values are compared to **computed** values from the data
    with a default tolerance of **10% of the metadata value**.

### - Global Attributes

**Required** (must exist and often must start with a given prefix):

    author_institution, author_name, Conventions (CF-1.0), data_revision (r*),
    id, Metadata_Conventions ("Unidata Dataset Discovery v1.0"),
    model, model_type, model_subtype, reference_pid ("doi"),
    reference, repository_institution ("EarthScope DS"),
    repository_name ("EMC"), repository_pid ("doi:"),
    summary, title, year

**Optional** (nice to have):

    author_email, author_url, grid_ref ("latitude_longitude"),
    version ("v"), grid_dim ("3D")

### - Variables Summary

Variables are grouped into: - **Coordinate Variables**
(incl. auxiliary); expected attributes depend on type: - Coordinates:
`long_name`, `units`, `standard_name` - Depth/Elevation: `long_name`,
`units`, `positive` - **Model Variables** (data variables): require
`long_name`, `units`, `display_name`.

For 3D model variables, it also validates **dimension order**. The guide
encodes a preferred order based on the two primary spatial dims and the
detected third dimension.

------------------------------------------------------------------------

## Example Session

    $ python emc_inspector.py


    emc_inspector.py - checks a given NetCDF file against the EMC requirements and flags any missing or incorrectly defined attributes. 
    The check is not exhaustive, and visual validation is required to identify any additional errors.


    Please enter the NetCDF file name: /data/models/ak135_emc_example.nc


    [INFO] emc_inspector.py version v2024.260

    === NetCDF File Check ===
    File: /data/models/ak135_emc_example.nc
    Size: 12.34 MB
    File Type: NetCDF-4 Classic Model ✔️

    === Dimensions and Sizes ===
    - longitude: 721 ✔️
    - latitude: 361 ✔️
    - depth: 100 ✔️

    === Coordinate Information ===
    Primary Coordinates: longitude, latitude, depth ✔️
    Third Dimension: depth ✔️
    Auxiliary Coordinates: time ✔️

    === Latitude and Longitude Value Check ===
    Longitude Format: -180/180 degrees ✔️
    Latitude Range: -89.75 to 89.75 ✔️
    Longitude Range: -180.00 to 180.00 ✔️

    === Geospatial Attributes Check ===
    geospatial_lat_min: Metadata = -89.75, Computed = -89.75 ✔️
    geospatial_lat_max: Metadata = 89.75, Computed = 89.75 ✔️
    geospatial_lon_min: Metadata = -180.0, Computed = -180.0 ✔️
    geospatial_lon_max: Metadata = 180.0, Computed = 180.0 ✔️
    geospatial_vertical_min: Metadata = 0.0, Computed = 0.0 ✔️
    geospatial_vertical_max: Metadata = 700.0, Computed = 700.0 ✔️

    === Global Attributes Check ===
    Global attributes included (values not checked): ✔️
    - author_name: Jane Doe
    - author_institution: Example University
    - Conventions: CF-1.0
    ...

    === CF Compliance Check ===
    CF Compliance: Compliant ✔️

    === Variables Summary ===
    Coordinate Variables:
    - latitude: Range = -89.750 to 89.750, Units = degrees_north, Dimensions = latitude ✔️
    - longitude: Range = -180.000 to 180.000, Units = degrees_east, Dimensions = longitude ✔️
    - depth: Range = 0.000 to 700.000, Units = km, Dimensions = depth ✔️

    Model Variables:
    - VS: Range = 3.000 to 4.900, Units = km/s, Dimensions = depth, latitude, longitude ✔️

    ============================
     Metadata Inspection Summary
    ============================
    File Type: NetCDF-4 Classic Model ✔️
    Dimensions:
    - longitude: 721
    - latitude: 361
    - depth: 100
    Longitude Format: -180/180 degrees
    Latitude Range: -89.750 to 89.750
    Longitude Range: -180.000 to 180.000
    ...

------------------------------------------------------------------------

## Exit Codes & Errors

-   The script currently **prints** errors and does not use structured
    exit codes for each failure.\
    For CI, you can detect failure patterns (`❌`, "FAIL", or "Error")
    in stdout.
-   If a file cannot be opened you'll see:
    `Error: Unable to open the file '...'.`

------------------------------------------------------------------------

## Interpreting Results

-   **✔️ OK** --- the item satisfies the EMC convention or an acceptable
    default.
-   **⚠️ Warning** --- optional or recommended item missing/ambiguous.
-   **❌ Fail** --- required item missing/invalid, or mismatch between
    metadata and data values.

> The tool is **not exhaustive**. Always perform visual inspection and
> domain‑specific checks for scientific validity.

------------------------------------------------------------------------

## 3. Troubleshooting

-   **Unicode icons print as \[OK\]/\[X\]/\[!\]**: Your terminal likely
    doesn't support UTF‑8; the script falls back automatically.
-   **"CF compliance check failed"**: Ensure `cfchecks` is installed and
    on your `PATH`.
-   **Longitude wrap confusion**: The script reports a "Potentially
    0--360 or --180/180" message when the range is ambiguous (e.g.,
    0--180). Inspect your data and set/verify the correct system.

------------------------------------------------------------------------

## 4. Contact / Feedback

For questions, comments, or feedback, please email: 📧 **data-help@earthscope.org**
