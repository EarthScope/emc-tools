# GeoCSV to NetCDF Converter User Guide

## Introduction

The `GeoCSV_2_netCDF.py` script converts a  [GeoCSV](https://ds.iris.edu/files/geocsv/GeoCSV.pdf) Earth model (2D or 3D, including projected grids) into a compressed [NetCDF-4 Classic](https://docs.unidata.ucar.edu/netcdf-c/4.9.2/file_format_specifications.html) (`NETCDF4_CLASSIC`) file that follows EMC conventions.  This ensures compatibility with EMC standards, improves data sharing, and supports long-term usability of Earth model files.  

---

## Usage Examples

```bash
# From the repository root
python src/GeoCSV_2_netCDF.py -i data/my_model.csv
# → writes data/my_model.nc (NETCDF4_CLASSIC)
```

Show headers only (no output .nc created):
```bash
python src/GeoCSV_2_netCDF.py -i data/my_model.csv -H
```

Verbose debugging:
```bash
python src/GeoCSV_2_netCDF.py -i data/my_model.csv -d
```

> **Output**: One NetCDF file named like the input (same base name), e.g., `my_model.nc`.

---

## Functionality Overview

- Parses a **GeoCSV** file that includes a metadata header (lines beginning with `#`), a header row, and data rows.  
- Detects **2D** (Y, X) or **3D** (Z, Y, X) grids and sorts input accordingly.  
- Creates **coordinate variables** (`y`, `x`, and optionally `z`) and **model variables** with attributes from the GeoCSV header.  
- Writes a **NETCDF4_CLASSIC** file with **zlib compression** (no chunking).  
- Preserves and extends provenance via the **`history`** global attribute.  

---

## System Requirements

- Python 3.8+  
- Packages:  
  - `netCDF4`  
  - `numpy`  

Install (example):  
```bash
pip install netCDF4 numpy
```

---

## Command-Line Usage

```
GeoCSV_2_netCDF.py -i FILE [-d] [-H]
```

| Option | Meaning |
|-------:|--------|
| `-i, --input` | **Required.** Path to the input GeoCSV file. |
| `-d, --debug` | Debug/verbose mode (prints details while processing). |
| `-H, --header` | **Header-only** mode. Parses headers and prints both GeoCSV and derived NetCDF headers, but **does not** write a `.nc` file. |

Exit codes: `0` on success; non‑zero for input/validation errors.  

---

## Input: GeoCSV Expectations

A **GeoCSV** file is a text file with two major parts:  

1. **Metadata header**: lines starting with `#` containing key–value pairs, e.g.:  
   ```
   # key: value
   # key_subkey: value
   ```  

2. **Tabular data**: a header row with column names, followed by rows of values.  

### Header Templates

When preparing a GeoCSV file, you should use one of the **header templates** provided in the `samples/` folder as a starting point. These templates include the required coordinate definitions, variable attributes, and formatting conventions expected by EMC-Tools:

- `samples/header_2D.csv` — template for 2D models (`y`, `x`)  
- `samples/header_3D.csv` — template for 3D models (`z`, `y`, `x`)  
- `samples/header_3D_projected.csv` — template for projected 3D grids (e.g., UTM with latitude/longitude auxiliary coordinates)  

These templates ensure that all required metadata keys (coordinates, delimiter, and variable attributes) are defined correctly and consistently.  


---

## 2D vs 3D Models

- The script infers **max dimensions** from variable blocks.  
- If any variable declares `dimensions: 3`, the dataset is treated as **3D** and expects a `z` column.  

Sorting:  
- **3D**: `(z, y, x)`  
- **2D**: `(y, x)`  

---

## Output: NetCDF Structure

- **File format**: `NETCDF4_CLASSIC`  
- **Compression**: zlib enabled, `complevel=4`, `shuffle=True` (no chunking)  
- **Dimensions**:  
  - 2D: `(y, x)`  
  - 3D: `(z, y, x)`  
- **Coordinate variables**: derived from `y`, `x`, and (if 3D) `z`.  
- **Model variables**: created from variable sections; attributes copied from GeoCSV header.  
- **Global attributes**: copied from `# global_*` entries; `history` always updated.  

> **Data type**: defaults to `f4` unless `VAR_DTYPE` is changed to `f8` in the script.  

---

## Typical GeoCSV Example (2D)

```text
# delimiter: ,
# x_column: longitude
# y_column: latitude

# vp_column: Vp
# vp_variable: vp
# vp_dimensions: 2
# vp_units: km/s
# vp_long_name: P-wave velocity
# vp_display_name: Vp
# vp__FillValue: -9999.0
# vp_missing_value: -9999.0

longitude,latitude,Vp
-123.0,45.0,6.5
-122.0,45.0,6.6
-123.0,46.0,6.4
-122.0,46.0,6.7
```

Run:  
```bash
python tools/GeoCSV_2_netCDF.py -i data/example_2d.csv
```

---

## Typical GeoCSV Example (3D)

```text
# delimiter: ,
# x_column: x
# y_column: y
# z_column: depth_km

# vp_column: Vp
# vp_variable: vp
# vp_dimensions: 3
# vp_units: km/s
# vp_long_name: P-wave velocity
# vp_display_name: Vp

x,y,depth_km,Vp
0,0,0,6.0
0,0,50,7.0
0,1,0,6.1
0,1,50,7.1
```

Run:  
```bash
python tools/GeoCSV_2_netCDF.py -i data/example_3d.csv
```

---

## Variable Mapping

For each variable tag (e.g., `vp`), the script uses:  

- `tag_column`: column name in the CSV  
- `tag_variable`: output NetCDF variable name  
- `tag_dimensions`: `2` or `3`  
- Additional attributes copied directly (e.g., `units`, `long_name`)  
- Special handling:  
  - `tag__FillValue` → `fill_value` in NetCDF  
  - `tag_missing_value` → attribute only  
  - `tag__range` → numeric array attribute  

---

## Coordinate Variables

- Created as 1D arrays: `y[y]`, `x[x]`, `z[z]` if 3D  
- Assigned from sorted unique values in data columns  
- Attributes may include `long_name`, `units`, `standard_name`, `_range`, `missing_value`, `_FillValue`  

---

## Global Attributes & Provenance

- All `# global_*: value` entries are written as NetCDF global attributes.  
- A `history` entry is appended:  
  ```
  YYYY-mm-dd HH:MM:SS UTC Converted to NetCDF by GeoCSV_2_netCDF.py vYYYY.DDD from my_model.csv
  ```  

---

## Output Location & Naming

- Input `/path/to/foo.csv` → Output `/path/to/foo.nc`  
- With `-H`, no file is written; headers are displayed only.  

---

## Additional Resources

- [**EMC User Guide**](../emc-user-guide.md)  
- [**EMC Model Files Standards and Conventions**](../reference/emc-standards-conventions.md)  
- [**NetCDF to GeoCSV Converter Guide**](emc-netcdf-2-geocsv-user-guide.md)  
- [**NetCDF to GeoJSON Converter Guide**](emc-netcdf-2-geojson-user-guide.md)  

---

**Comments or Questions?**  
For any questions or feedback about EMC Earth models or EMC-Tools,  
please email: **[data-help@earthscope.org](mailto:data-help@earthscope.org)**  
