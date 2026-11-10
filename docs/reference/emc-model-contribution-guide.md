# EMC Model Contribution Guidelines

## Introduction

The **Earth Model Collaboration (EMC)** is a **community-supported repository for Earth models**, created to facilitate the sharing, previewing, and access of diverse geophysical Earth models. Authors of Earth models are **strongly encouraged** to contribute their models to the repository, supporting greater visibility and accessibility within the research community.  

If you are interested in contributing, please follow the guidelines below to ensure your model can be seamlessly integrated into EMC.

---

## Submission Checklist

**All contributed models must be described in a peer-reviewed publication.**  
Please include DOI of the associated peer-reviewed article.  

Before submitting your model to EMC, review the [EMC Standards and Conventions](emc-standards-conventions.md) and confirm that your submission meets the following requirements:

!!! tip "Using EMC-Tools to create your netCDF file"
    Model netCDF files created with the EMC-Tools [GeoCSV to netCDF converter](../how-to/emc-geocsv-2-netcdf-user-guide.md) automatically follow the required EMC data and metadata format (items 1–5 below). If you use the converter, you can skip directly to validating and previewing your file (items 6–7).

1. Save model data in **NetCDF-4 Classic** format (use compression for large files).  
2. Ensure that metadata is **CF-compliant**.  
3. Include all **required global attributes** (Note: See the header templates under **samples/** for the required attributes).  
4. Document the **coordinate system** (geographic or projected).  
5. Verify that **geospatial bounds** are correctly defined in metadata.  
6. Validate metadata compliance using the [EMC Metadata Inspector](../how-to/emc-inspector-user-guide.md).  
7. Preview and verify model data integrity with the [EMC Model Explorer](../how-to/emc-explorer-user-guide.md).  
8. Complete the **EMC Model Page Template** ([EMC-ModelPageTemplate.toml](../samples/EMC-ModelPageTemplate.toml)) describing your model. EMC uses this file to build your model's page in the catalog. Instructions are included as comments in the template. See [Model Page Template](#model-page-template) below.  
9. Submit your final model file(s), completed model page template, figure files, and publication reference to:  
   📧 [help@earthscope.org](mailto:help@earthscope.org).  

---

## Model Page Template

Every model in EMC has its own page in the catalog showing the model description, authors, coverage, figures, download links, and citation. To create this page, EMC needs a completed **model page template** from you.

📥 **[Download EMC-ModelPageTemplate.toml](../samples/EMC-ModelPageTemplate.toml)**

The template is a plain-text [TOML](https://toml.io/) file that you can edit in any text editor. It is pre-filled with example values for a fictional model; replace them with your own information. Every field is labeled in the file:

- **[REQUIRED]**: you must fill this in.
- **[OPTIONAL]**: leave as `""` if it does not apply.
- **[EMC STAFF]**: completed by EMC; leave unchanged.

When completing the template:

- Save your copy using your model name, for example `Doe2024-Vs.toml`.
- Do not rename, delete, or reorder section headers (`[model]`, `[[authors]]`, etc.) or field names.
- Name your figure files to match the `file` entries in the `[[figures]]` blocks (for example, `Doe2024-Vs_1.png`) and include them with your submission. PNG or JPG, roughly 1000–2000 pixels wide, is preferred.
- Give DOIs without the `https://doi.org/` prefix and without a trailing period.

If anything in the template is unclear, add a comment (a line starting with `#`) next to the field and EMC staff will follow up.
