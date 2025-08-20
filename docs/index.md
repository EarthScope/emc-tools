# EMC User Guide

## Introduction

The [Earth Model Collaboration (EMC)](https://ds.iris.edu/ds/products/emc/) is a **community-supported repository of Earth models** that enables researchers to share, preview, and access diverse geophysical Earth models. To support this effort, we provide [EMC-Tools](https://github.com/EarthScope/emc-tools), a suite of **Python 3 utilities**, that streamline the **storage**, **extraction**, and **visualization** of these models. By adopting a **uniform, self-contained data format**, the tools ensure consistency across contributions and promote the long-term usability of Earth models within the community.  

---

## Model Files 

This section describes the structure of EMC model files and the ways they can be accessed.

### Standards and Conventions

All EMC models follow a set of **metadata and format standards** to ensure compatibility across tools. See the [Model Files Standards and Conventions](reference/emc-standards-conventions.md) for complete details.

### Access

EMC model files are available through both direct downloads and interactive tools. The resources below make it easy to obtain NetCDF files from the EarthScope data archive and to explore models dynamically with a filtering and access API.

- [**EMC Model File Server**](https://data.dev.earthscope.org/archive/seismology/products/emc/netcdf/) &mdash; access NetCDF model files directly from the EarthScope data archive.  
- [**EMC Model Filtering & Access API**](https://data.dev.earthscope.org/archive/seismology/products/emc/README.html) &mdash; interactively filter, explore, visualize, and download Earth model files.  

---

## EMC-Tools

The EMC-Tools suite provides utilities for validating metadata, exploring models, and converting between common geophysical data formats.

### Tool Guides

The guides below provide step-by-step instructions for using the core EMC-Tools utilities.

- [**EMC Metadata Inspector**](how-to/emc-inspector-user-guide.md) &mdash; validate a NetCDF file against EarthScope’s EMC conventions and common CF metadata expectations.   
- [**EMC Model Explorer**](how-to/emc-explorer-user-guide.md) &mdash; an interactive command-line tool for exploring, visualizing, and exporting Earth model data from **local** NetCDF files.  

### Format Converter Guides

The following guides cover utilities for format conversion:

- [**GeoCSV to NetCDF Converter**](how-to/emc-geocsv-2-netcdf-user-guide.md) &mdash; convert a GeoCSV Earth model (2D or 3D, including projected grids) into a compressed NetCDF-4 Classic file that follows EMC conventions.  
- [**NetCDF to GeoCSV Converter**](how-to/emc-netcdf-2-geocsv-user-guide.md) &mdash; convert a 2D or 3D NetCDF Earth model into GeoCSV 2.0.  
- [**NetCDF Metadata to GeoJSON Converter**](how-to/emc-netcdf-2-geojson-user-guide.md) &mdash; convert a NetCDF Earth model file into a GeoJSON footprint containing both the geographic bounding box and key metadata.  

### Getting Started

Follow these steps to begin using EMC-Tools:

1. **Install** Python 3.8+ and dependencies from `requirements.txt`.  
2. **Explore** the Tool Guides listed above.  
3. **Test** with sample datasets in the `samples/` directory.  
4. **Follow** EMC standards when preparing or submitting models.  

> **Tip:** EMC-Tools is actively maintained.  
> Always check the [GitHub repository](https://github.com/EarthScope/emc-tools) for the latest updates before using the tools.

---

## Additional Resources

The following resources provide standards, guidelines, and workflows for contributing and maintaining EMC models.

- [**Model Files Standards and Conventions**](reference/emc-standards-conventions.md) &mdash; define supported file formats, metadata requirements, and coordinate systems to ensure consistency across EMC model files.  
- [**Contribute Models to EMC**](reference/emc-model-contribution-guide.md) &mdash; provide guidelines for contributing new models to EMC.  
- [**Convert CSV Files to GeoCSV**](reference/emc-csv-to-geocsv.md) &mdash; describe how to convert a plain CSV table into a valid GeoCSV file that EMC tools can ingest.  
- [**Update an Existing Model**](reference/emc-update-existing-model.md) &mdash; explain how to update an existing Earth model in EMC.  

---

**Comments or Questions?**  
For any questions or feedback about EMC Earth models or EMC-Tools,  
please email: **[data-help@earthscope.org](mailto:data-help@earthscope.org)**  
