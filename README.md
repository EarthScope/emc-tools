# EarthScope

Data Services (DS)  
Data Products Unit  
EarthScope Earth Model Collaboration (EMC) - Release 3.0

### COMMENTS/QUESTIONS:
Please contact [data-help@earthscope.org](mailto:data-help@earthscope.org)

---

## Release Notes

- **2025-08-20 (Release 3.0):**  
  - All netCDF output now uses **NetCDF4 Classic** format for compatibility and long-term support.  
  - Introduces **emc_explorer.py**, a tool for interactive browsing and visualization of EMC models.  
  - Enhanced **metadata validation** in `emc_inspector.py`, ensuring stricter compliance with EMC and CF standards.  

- **2024-09-16 (Release 2.1):**  
  Introduced **emc_inspector.py** (v2024.260) to help identify potential issues with netCDF metadata.  

- **Release 2.0:**  
  The R2 release included changes to the GeoCSV fields. R1 GeoCSV files **are not** compatible with R2 tools.  
  However, R1-generated netCDF files can be used to produce R2-style GeoCSV files via **EMC Tools Release 2.0**.  

---

## 1. DESCRIPTION

The **EMC Tools** provide Python utilities for:  

- Converting EMC netCDF Earth model files **to and from GeoCSV**.  
- **Validating metadata** against EMC requirements.  
- **Exploring and visualizing** model content interactively.  

### Available Tools
- **netCDF_2_GeoCSV.py** – Converts a 2D or 3D netCDF Earth model file (geographic or projected coordinates) to GeoCSV format or displays header info.  
- **GeoCSV_2_netCDF.py** – Converts a 2D or 3D GeoCSV Earth model file to netCDF4 Classic format or displays header info.  
- **emc_inspector.py** – Checks a netCDF file against EMC requirements and flags missing or incorrectly defined attributes.  
- **emc_explorer.py** – Interactive exploration of EMC netCDF files, including variable inspection, coordinate validation, and quick visual previews.  

### References
- [EMC (Earth Model Collaboration)](http://ds.iris.edu/ds/products/emc/)  
- [netCDF Documentation](https://www.unidata.ucar.edu/software/netcdf/)  
- [GeoCSV Format Specification](http://geows.ds.iris.edu/documents/GeoCSV.pdf)  

---

## 2. USAGE

Each tool can be run independently. Use the `-h` option for help:  

- `netCDF_2_GeoCSV.py -h`  
- `GeoCSV_2_netCDF.py -h`  
- `emc_inspector.py` (run without arguments for metadata validation prompt)  
- `emc_explorer.py -h`  

---

### 2.1 netCDF → GeoCSV
Convert a netCDF file into GeoCSV.  
```bash
netCDF_2_GeoCSV.py -i FILE -x longitude -y latitude -z depth -m mode -d -H
```

**Options:**  
- `-i [required]`: Input netCDF Earth model file  
- `-x [required]`: x-axis variable  
- `-y [required]`: y-axis variable  
- `-z [required]`: z-axis variable (`none` for 2D models)  
- `-m [default: single]`: Output mode (`single` = one GeoCSV file, `depth` = one file per depth)  
- `-d [default: False]`: Debug mode  
- `-H [default: False]`: Display headers only  

_Note: Only `.csv` files converted with `-m single` can be converted back to netCDF._  

---

### 2.2 GeoCSV → netCDF
Convert a GeoCSV file back into netCDF4 Classic format.  
```bash
GeoCSV_2_netCDF.py -i FILE -d -H
```

**Options:**  
- `-i [required]`: Input GeoCSV Earth model file  
- `-d [default: False]`: Debug mode  
- `-H [default: False]`: Display headers only  

_Note: Only a single `.csv` file can be converted back to netCDF._  

---

### 2.3 Metadata Inspector
Run the metadata inspector:  
```bash
emc_inspector.py
```
It will prompt:  
```
Please enter the netCDF file name:
```
Then prints detailed metadata validation results.  

---

### 2.4 EMC Explorer
Launch the EMC Explorer for interactive browsing:  
```bash
emc_explorer.py -i FILE
```

Features include:  
- Listing model dimensions, coordinates, and variables.  
- Quick metadata summaries.  
- Optional visualization of slices or grids (where supported).  

---

### 2.5 Creating EMC GeoCSV Model Files
To create EMC GeoCSV Earth model files:  

1. Follow EMC’s [Contribution Guide](http://ds.iris.edu/ds/products/emc-contributionguide/).  
2. Review the [GeoCSV format](http://geows.ds.iris.edu/documents/GeoCSV.pdf).  
3. Use the **samples** directory as templates:  
   - `samples/header_2D.csv`  
   - `samples/header_3D.csv`  
   - `samples/header_3D_projected.csv`  
4. Ensure your data is sorted in the correct coordinate order:  
   - **3D:** depth, latitude, longitude  
   - **2D:** latitude, longitude  
