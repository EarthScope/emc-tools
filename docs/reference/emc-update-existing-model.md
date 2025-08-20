# How to Update an Existing Model in the EarthScope Earth Model Collaboration (EMC)

This guide explains the process for updating an existing Earth model in EMC.

---

## 1. Determine the Version/Revision Number
If the model data has changed, you must decide on an appropriate version or revision number before submission.  
- **Version**: Use when significant scientific changes are made (e.g., new methodology, added datasets).  
- **Revision**: Use for corrections, minor updates, or formatting changes.

Example:
- Previous: `r0.0`
- Minor revision: `r0.1`
- Major update: `r1.0`

---

## 2. Document the Changes
Provide a short paragraph summarizing what has changed in the new version. This will be included in the EMC model description page.

Example:
> *This update corrects a sign error in the Vs parameter at depths greater than 100 km and includes new measurements for the western region based on 2024 seismic data.*

---

## 3. Validate Metadata with the EMC Inspector
Run the **EMC Inspector** tool to ensure the new NetCDF file's metadata follows EMC and CF conventions.  
- Check that all required global and variable attributes are present.  
- Verify that `data_layout` and `source` attributes are properly set.  
- Ensure `geospatial` attributes match the coordinate ranges in the data.

Command example:
```bash
python emc_inspector.py your_model_file.nc
```

---

## 4. Validate Data Content with the EMC Explorer
Use the **EMC Explorer** tool to visually inspect and verify that the data values are correct.  
- Plot key variables.  
- Check for unexpected gaps or anomalies.  
- Verify coordinate grids and projections.

Command example:
```bash
python emc_explorer.py your_model_file.nc
```

---

## 5. Send the Update to EarthScope
Once validated:

1. Attach the **updated NetCDF file**.
2. Include your **change description** paragraph.
3. State the **new version/revision number**.
4. Send everything to: **data-help@earthscope.org**

**Important:**  
If you have already communicated about this model via `data-help@earthscope.org`, **always reply to one of the existing email threads** instead of starting a new one. This ensures your update is tracked under the same ticket.

---

**Following these steps will help ensure smooth and accurate updates to your EMC model.**
