# EMC Model Contribution Guidelines

## Introduction

The **Earth Model Collaboration (EMC)** is a **community-supported repository for Earth models**, created to facilitate the sharing, previewing, and access of diverse geophysical Earth models. Authors of Earth models are **strongly encouraged** to contribute their models to the repository, supporting greater visibility and accessibility within the research community.  

If you are interested in contributing, please follow the guidelines below to ensure your model can be seamlessly integrated into EMC.

---

## Submission Checklist

**All contributed models must be described in a peer-reviewed publication.**  
Please include DOI of the associated peer-reviewed article.  

Before submitting your model to EMC, review the [EMC Standards and Conventions](emc-standards-conventions.md) and confirm that your submission meets the following requirements:

1. Save model data in **NetCDF-4 Classic** format (use compression for large files).  
2. Ensure that metadata is **CF-compliant**.  
3. Include all **required global attributes** (Note: See the header templates under **samples/** for the required attributes). 
1. Document the **coordinate system** (geographic or projected).  
2. Verify that **geospatial bounds** are correctly defined in metadata.  
3. Validate metadata compliance using the [EMC Metadata Inspector](../how-to/emc-inspector-user-guide.md).  
4. Preview and verify model data integrity with the [EMC Model Explorer](../how-to/emc-explorer-user-guide.md).  
5. Submit your final model file and publication reference to:  
   📧 [data-help@earthscope.org](mailto:data-help@earthscope.org).  
