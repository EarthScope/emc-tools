# EMC Model Files Standards and Conventions

## Introduction

The **EMC Standards and Conventions** define the supported file formats, metadata requirements, and coordinate systems to ensure a consistent approach to creating, distributing, and using EMC model files. Adhering to these conventions improves the **consistency, interoperability, and reliability** of contributed models, enabling more effective collaboration and data sharing across the community.  

This guide highlights the critical role of metadata, specifying the required attributes, their formats, and guidelines for coordinate system usage. Following these practices helps minimize distortion, improve accuracy, and maintain the overall integrity and usability of EMC data and tools.  

Contributed Earth models are available for download from the [EarthScope Data Repository](https://data.dev.earthscope.org/archive/seismology/products/emc/netcdf/). For interactive filtering and selection of models, visit the [EMC Model Filtering & Access API](https://data.dev.earthscope.org/archive/seismology/products/emc/README.html).  

---

## File Format Standards

This section describes the required formats for EMC model files and the benefits of adopting the current standards.

### Required Format
- All EMC models must use **NetCDF-4 Classic** as the required container format. Earlier versions of the standard required **NetCDF-3 Classic**.  
- Apply **compression** for large models to reduce file size and improve data transfer efficiency.  
- Ensure all models are **CF-compliant** (see [CF Metadata Conventions](https://cfconventions.org/)) and compatible with EMC-Tools.  

### Format Conversion
- Use the [**EMC-Tools**](https://github.com/EarthScope/emc-tools) utilities to convert Earth model files between [GeoCSV](https://ds.iris.edu/files/geocsv/GeoCSV.pdf) and [NetCDF-4 Classic](https://docs.unidata.ucar.edu/netcdf-c/4.9.2/file_format_specifications.html).  
- The transition from **NetCDF-3 Classic** to **NetCDF-4 Classic** provides key benefits:  
    - Enable internal compression and chunking to reduce storage and transfer costs.  
    - Maintain backward compatibility with the NetCDF-3 data model to ensure broad software support.  
    - Adopt a modern container format that improves long-term sustainability and usability of EMC model files.  

---

## Metadata Requirements

All EMC NetCDF files must include metadata following the **CF Metadata Conventions**. In addition, EMC introduces the following required attributes:

### `source` (Variable-Level Attribute) 
Each model data variable must include a `source` attribute describing the origin of the data:

- If derived from observations:  

```  
vp:source = data-derived  
```  

- If derived from an empirical formula:  
  Include the formula or reference. Example:  

```  
vs:source = assumed to follow vp as vs = vp / sqrt(3)  
```  

### `data_layout` (Global Attribute) 
Specifies the layout of data in the file:

- `vertex` — values specified at vertices  
- `cell` — values specified at the centers of grid cells  

If `data_layout` is not defined, the default `vertex` is assumed.  

---

## Coordinate Systems (`grid_ref`)

The `grid_ref` attribute specifies the coordinate system used in the model file:

- The default coordinate system is geographic:  

```  
grid_ref = latitude_longitude  
```  

- Models using a projected coordinate system (e.g., UTM) must include variables for both the projected coordinates and the geographic latitude/longitude coordinates.  
- The value of `grid_ref` should indicate the projection used, and metadata must clearly document the projection, including EPSG codes where applicable.  

---

## Additional Resources

The following resources provide additional standards, guidelines, and workflows for contributing and maintaining EMC models:

- [**EMC User Guide**](../emc-user-guide.md) &mdash; overview of EMC, how to access model files, and how to use EMC-Tools.  
- [**Contribute Models to EMC**](emc-model-contribution-guide.md) &mdash; guidelines for contributing new models to EMC.  
- [**Convert CSV Files to GeoCSV**](emc-csv-to-geocsv.md) &mdash; instructions for converting a plain CSV table into a valid GeoCSV file that EMC tools can ingest.  
- [**Update an Existing Model**](emc-update-existing-model.md) &mdash; explain how to update an existing Earth model in EMC.  

---

**Comments or Questions?**  
For any questions or feedback about EMC Earth models or EMC-Tools,  
please email: **[data-help@earthscope.org](mailto:data-help@earthscope.org)**  
