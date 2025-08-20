# NetCDF to GeoJSON Converter User Guide

## Introduction

The `netCDF_2_GeoJSON.py` script converts a NetCDF Earth model file into a **GeoJSON footprint** containing both the geographic bounding box and key metadata.  
This enables quick visualization of model coverage and inspection of metadata in a lightweight, widely supported format.  

---

## Usage Examples

### Basic command

```bash
python src/netCDF_2_GeoJSON.py -i /path/to/model.nc
```

### Command-line flags

| Flag | Long form   | Required | Value                 | Description                       |
|------|-------------|----------|-----------------------|-----------------------------------|
| `-i` | `--input`   | ✅        | Path to a NetCDF file | Input NetCDF Earth model file     |
| `-h` | `--help`    | ❌        | —                     | Show usage information and exit   |

### Example

```bash
python src/netCDF_2_GeoJSON.py -i samples/KEA20-Moho.r0.0.nc
# → writes KEA20-Moho.r0.0.geojson in the current directory
```

---

## Output Structure

- **Filename**: `{input_stem}.geojson`  
  Example: `samples/KEA20-Moho.r0.0.geojson`  
- **Format**: GeoJSON FeatureCollection containing one Feature.  

### Example (abridged)

```json
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "geometry": {
        "type": "Polygon",
        "coordinates": [[
          [lon_min, lat_min],
          [lon_max, lat_min],
          [lon_max, lat_max],
          [lon_min, lat_max],
          [lon_min, lat_min]
        ]]
      },
      "properties": {
        "model": "...",
        "title": "...",
        "summary": "...",
        "reference": "...",
        "reference_pid": "...",
        "repository_doi": "...",
        "year": "...",
        "data_revision": "r0.0",
        "version": "v0.0",
        "model_type": "Tomography Earth Model",
        "model_subtype": "Radially anisotropic S velocities (km/s) and Moho depth (km)",
        "filename": "KEA20-Moho.r0.0.nc",
        "size_mb": 0.19,
        "format": "NETCDF4_CLASSIC",
        "grid_ref": "latitude_longitude",
        "projection_type": "geographic_wgs84",
        "projection_notes": "Unprojected geographic coordinates, WGS 84 datum",
        "grid_dim": "2D",
        "geospatial_lat_min": 15.0,
        "geospatial_lat_max": 51.0,
        "geospatial_lat_resolution": 0.20,
        "geospatial_lat_units": "degrees_north",
        "geospatial_lon_min": 70.0,
        "geospatial_lon_max": 150.0,
        "geospatial_lon_resolution": 0.20,
        "geospatial_lon_units": "degrees_east",
        "geospatial_vertical_min": 0.0,
        "geospatial_vertical_max": 0.0,
        "geospatial_vertical_units": "",
        "variables": {
          "moho": "Moho depth relative to a mean Earth radius of 6371 km"
        }
      }
    }
  ]
}
```

**Note:** The footprint polygon is assumed to be in **geographic (WGS84)** coordinates. The script does not perform reprojection; it uses metadata attributes directly.  

---

## Additional Resources

- [**EMC User Guide**](../index.md)  
- [**EMC Model Files Standards and Conventions**](../reference/emc-standards-conventions.md)  
- [**GeoCSV to NetCDF Converter Guide**](emc-geocsv-2-netcdf-user-guide.md)  
- [**NetCDF to GeoCSV Converter Guide**](emc-netcdf-2-geocsv-user-guide.md)  

---

**Comments or Questions?**  
For any questions or feedback about EMC Earth models or EMC-Tools,  
please email: **[data-help@earthscope.org](mailto:data-help@earthscope.org)**  