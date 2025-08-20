# NetCDF to GeoCSV Converter User Guide

## Introduction

The `netCDF_2_GeoCSV.py` script converts a **2D or 3D NetCDF Earth model** into **GeoCSV 2.0**. It supports writing a single consolidated file or one file per depth level for 3D datasets.  
This ensures EMC models can be distributed in a lightweight, tabular format while preserving key metadata and attributes.  

---

## Usage Examples

### Basic command

```bash
python src/netCDF_2_GeoCSV.py -i samples/CSEM2-Africa.v2024.12.01.r0.0-n4c.nc -x longitude -y latitude -z depth
```

- Coordinate variable flags must match variables in your NetCDF file:
  - `-x <X_VARIABLE>` (e.g., `longitude` or `x`)  
  - `-y <Y_VARIABLE>` (e.g., `latitude` or `y`)  
  - `-z <Z_VARIABLE>` (e.g., `depth`) or `-z none` for 2D models  

### Command-line flags

| Flag | Long form | Required | Default | Description |
|------|-----------|----------|---------|-------------|
| `-i` | `--input` | ✅ | — | Path to input NetCDF model. Looks in CWD, then `data/` if not found. |
| `-x` | `--xvar` | ✅ | — | Name of the x-axis variable (e.g., `longitude`, `x`). |
| `-y` | `--yvar` | ✅ | — | Name of the y-axis variable (e.g., `latitude`, `y`). |
| `-z` | `--zvar` | ✅ | — | Name of the z-axis variable (e.g., `depth`). For 2D, may be set to `none`. |
| `-m` | `--mode` |  | `single` | Output mode: `single` (all depths in one file) or `depth` (one file per depth, 3D only). |
| `-d` | `--debug` |  | `False` | Enable verbose debug output. |
| `-H` | `--header` |  | `False` | Print NetCDF and GeoCSV headers only, then exit. |

Exit codes: `0` on success; non-zero for input or validation errors.  

---

## Examples

### 1. Single-file output (3D)

```bash
python netCDF_2_GeoCSV.py -i samples/KEA20-Moho.r0.0.nc -x longitude -y latitude -z depth -m single
```
**Output:** `models/MODEL_3D.csv`

### 2. One file per depth (3D)

```bash
python netCDF_2_GeoCSV.py -i samples/KEA20-Moho.r0.0.nc -x longitude -y latitude -z depth -m depth
```
**Output:** `models/MODEL_3D_<DEPTH>_km.csv` for each depth level

### 3. 2D model

```bash
python src/netCDF_2_GeoCSV.py -i samples/KEA20-Moho.r0.0_1.nc -x latitude -y longitude -m single
```
**Output:** `models/MODEL_2D.csv`

### 4. Show headers only (2D)

```bash
python src/netCDF_2_GeoCSV.py -i samples/KEA20-Moho.r0.0.nc -x longitude -y latitude -H
```
Prints both the NetCDF header (dimensions, variables, globals) and the GeoCSV header that would be written, then exits.  

---

## Output Structure

### GeoCSV Header (2.0)

The script writes a full header that includes:  

- Dataset metadata:  
  ```text
  # dataset: GeoCSV2.0
  # created: <UTC time> (netCDF_2_GeoCSV.py)
  # netCDF_file: <basename>
  # delimiter: |
  ```  
- Global attributes: written as `# global_<attr>: <value>` (with a `global_history` entry noting the conversion).  
- For each coordinate (`x`, `y`, `z` if 3D) and each model variable, the header includes:  
  - `_column`, `_variable`, `_dimensions`  
  - All NetCDF attributes (e.g., `long_name`, `units`, `standard_name`, `_FillValue`, `_range` if present).  

> `_range` values are written as `min,max`.  

### GeoCSV Data Section

- **3D / single mode:**  
  ```
  <y>|<x>|<z>|<var_1>|...|<var_N>
  ```
- **3D / depth mode:**  
  ```
  # depth: <value>
  <y>|<x>|<var_1>|...|<var_N>
  ```
- **2D:**  
  ```
  <y>|<x>|<var_1>|...|<var_N>
  ```  

- Delimiter is a **pipe** (`|`) to avoid conflicts with commas in float values.  
- All values are written as strings to preserve precision. Missing values are written as `nan`.  

---

## Additional Resources

- [**EMC User Guide**](../index.md)  
- [**EMC Model Files Standards and Conventions**](../reference/emc-standards-conventions.md)  
- [**GeoCSV to NetCDF Converter Guide**](emc-geocsv-2-netcdf-user-guide.md)  
- [**NetCDF to GeoJSON Converter Guide**](emc-netcdf-2-geojson-user-guide.md)  

---

**Comments or Questions?**  
For any questions or feedback about EMC Earth models or EMC-Tools,  
please email: **[data-help@earthscope.org](mailto:data-help@earthscope.org)**  
