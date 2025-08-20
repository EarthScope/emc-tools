# Convert CSV to GeoCSV (EMC)

This guide shows how to turn a plain **CSV** table into a valid **GeoCSV** file that EMC tools can ingest — using the sample header templates in `samples/`.

> **GeoCSV = header (lines starting with `#`) + your CSV data**.  
> You can usually make one by concatenating a prepared header file and your CSV body.

---

## What is GeoCSV (in EMC)?

GeoCSV is a simple, text-based format:
- A **header block** of metadata lines beginning with `# key: value`
- Followed by a **data table** (comma-separated by default)
- For gridded Earth models, the header declares which columns are coordinates (`x`, `y`, optional `z`) and which columns are **model variables** (e.g., `vp`, `vs`, `rho`).

This repository provides three header templates:

- `samples/header_2D.csv` – 2D models (y, x + variables)
- `samples/header_3D.csv` – 3D models (y, x, z + variables)
- `samples/header_3D_projected.csv` – 3D with a projected CRS (e.g., UTM) **and** geographic lat/lon variables

Each template already includes the required **GeoCSV keys**, for example:
```
# dataset: GeoCSV2.0
# delimiter: ,
# y_column: latitude
# x_column: longitude
# z_column: depth             <-- only in 3D
# vp_column: vp               <-- example model variable
# vp_units: km/s
...
latitude,longitude,depth,vp,vs,rho    <-- FIRST non-# line is the CSV header row
```

---

## Prerequisites

- Your CSV data file (no `#` lines) — call it `my_data.csv`
- One of the header templates from `samples/`
- UTF‑8 encoding, Unix line endings (`\n`) recommended
- **Column names in the first row must match the header’s `*_column:` values**

Optional, for validation:
- `GeoCSV_2_netCDF.py` (in EMC_Tools) to convert GeoCSV → NetCDF for quick checks
- `emc_inspector.py` to verify NetCDF compliance

---

## Quick recipe (concatenate header + data)

Pick the right header and concatenate with your data:

### 2D model
```bash
# Inspect the header template
sed -n '1,80p' samples/header_2D.csv

# Ensure your CSV header row matches the template’s *_column names
head -1 my_data.csv

# Build GeoCSV
cat samples/header_2D.csv my_data.csv > out_2D.geocsv.csv
```

### 3D model (geographic lat/lon/depth)
```bash
sed -n '1,80p' samples/header_3D.csv
head -1 my_data.csv
cat samples/header_3D.csv my_data.csv > out_3D.geocsv.csv
```

### 3D model (projected CRS + lat/lon)
```bash
sed -n '1,120p' samples/header_3D_projected.csv
head -1 my_data.csv
cat samples/header_3D_projected.csv my_data.csv > out_3D_proj.geocsv.csv
```

> **Tip:** The output filename extension is arbitrary (`.csv` is fine). EMC tools look at the **header keys**, not the file extension.

---

## Mapping your CSV columns

Your CSV’s **first row** must contain column names. The header’s `*_column:` keys must point to those exact names. Examples:

### 2D
- `# y_column: latitude`
- `# x_column: longitude`
- Model variables:
  - `# vp_column: vp`
  - `# vp_units: km/s`
  - `# vp_long_name: P-wave velocity`

### 3D
- `# y_column: latitude`
- `# x_column: longitude`
- `# z_column: depth` (positive **down** usually; make sure units are clear)
- Variables as above

### 3D projected
- **Projected axes**:
  - `# y_column: y` (easting or northing depending on your convention)
  - `# x_column: x`
  - `# z_column: depth`
- **Geographic variables also present** so the tools can convert/plot:
  - `# latitude_column: latitude`
  - `# longitude_column: longitude`
- **CRS metadata** (in the header):
  - `# grid_ref: projected` (and `# utm_zone: 11N`, `# ellipsoid: WGS84`, etc.)
  - `# geospatial_*_units: meters` for x/y if they are projected

Edit the template to match your column names, units, and descriptions.

---

## Validate your GeoCSV

### 1) Convert to NetCDF (fast structural check)

```bash
# Show header only (no file written)
python GeoCSV_2_netCDF.py -i out_2D.geocsv.csv -H

# Write a NetCDF4 Classic file
python GeoCSV_2_netCDF.py -i out_2D.geocsv.csv
# → produces out_2D.geocsv.nc (or <input>.nc)
```

### 2) Inspect the resulting NetCDF (optional)

```bash
# Quick look using ncdump (if available)
ncdump -h out_2D.geocsv.nc

# EMC inspector (checks EMC/CF requirements)
python emc_inspector.py
# (You'll be prompted for the file name: enter: out_2D.geocsv.nc)
```

If `emc_inspector.py` flags issues (e.g., missing units), update the **GeoCSV header** and re-run the conversion.

---

## Example end-to-end workflows

### A) 2D model (lat/lon)

1. Prepare `my_data.csv` with columns:
   ```text
   latitude,longitude,vp,vs
   34.5,-117.0,6.10,3.50
   34.6,-117.1,6.05,3.48
   ...
   ```
2. Edit `samples/header_2D.csv` to ensure:
   ```text
   # y_column: latitude
   # x_column: longitude
   # vp_column: vp
   # vs_column: vs
   # vp_units: km/s
   # vs_units: km/s
   # delimiter: ,
   ```
3. Build and validate:
   ```bash
   cat samples/header_2D.csv my_data.csv > my_model_2D.csv
   python GeoCSV_2_netCDF.py -i my_model_2D.csv
   python emc_inspector.py   # then enter: my_model_2D.nc
   ```

### B) 3D model (lat/lon/depth)

1. `my_data.csv`:
   ```text
   latitude,longitude,depth,vp
   34.5,-117.0,0.0,5.90
   34.5,-117.0,5.0,6.05
   ...
   ```
2. Edit `samples/header_3D.csv`:
   ```text
   # z_column: depth
   # depth_units: km
   # geospatial_vertical_positive: down
   ```
3. Concatenate and convert:
   ```bash
   cat samples/header_3D.csv my_data.csv > my_model_3D.csv
   python GeoCSV_2_netCDF.py -i my_model_3D.csv
   python emc_inspector.py   # enter: my_model_3D.nc
   ```

### C) 3D projected (e.g., UTM 11N, meters)

1. `my_data.csv`:
   ```text
   x,y,depth,latitude,longitude,vp
   379000,3768000,0.0,34.50,-117.00,6.10
   379000,3768000,5.0,34.50,-117.00,6.18
   ...
   ```
2. Edit `samples/header_3D_projected.csv`:
   ```text
   # grid_ref: projected
   # utm_zone: 11N
   # ellipsoid: WGS84
   # x_units: meters
   # y_units: meters
   # latitude_column: latitude
   # longitude_column: longitude
   ```
3. Concatenate and convert:
   ```bash
   cat samples/header_3D_projected.csv my_data.csv > my_model_3D_proj.csv
   python GeoCSV_2_netCDF.py -i my_model_3D_proj.csv
   python emc_inspector.py   # enter: my_model_3D_proj.nc
   ```

---

## Common pitfalls & fixes

- **Header/data delimiter mismatch**: If your CSV uses `;` or `|`, set `# delimiter: ;` (or `|`) in the header.
- **Column names don’t match**: The header’s `*_column` values **must** match the CSV header row exactly.
- **Units missing**: Provide `<var>_units` key for every variable and coordinate.
- **Depth sign/origin**: Clarify with `geospatial_vertical_positive: down` (commonly for depth), and set `depth_units`.
- **Projected models without lat/lon**: For EMC tools, include **both** projected coordinates **and** latitude/longitude.
- **Encoding/BOM**: Save in UTF‑8 **without BOM**.
- **Sorting**: Not required for GeoCSV, but `GeoCSV_2_netCDF.py` will sort when building gridded arrays.


---

## Next steps

- Convert GeoCSV → NetCDF: `python GeoCSV_2_netCDF.py -i my_model.csv`
- Validate the NetCDF: `python emc_inspector.py` (enter the produced `.nc` file)
- Contribute the model following: `docs/how-to/model-contribution-guide.md`

