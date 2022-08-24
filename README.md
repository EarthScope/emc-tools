 
Incorporated Research Institutions for Seismology (IRIS):

    Data Management Center (DMC)
    Data Products Team
    IRIS Earth Model Collaboration (EMC) - Tools, R2

COMMENTS/QUESTIONS:

    Please contact manoch@iris.washington.edu


 2022-08-25
----------------------------------------------------------------------------------------------------------------------------------------------------------------------

 DESCRIPTION:

 EMC Tools contains two Python scripts for converting the EMC's netCDF Earth model files (in netCDF 3 format) to 
 and from GeoCSV format:
 
    - netCDF_2_GeoCSV.py - a Python script to read a 2D or 3D netCDF Earth model file in geographic or projected 
      coordinate system and display its header information or convert it to GeoCSV format.

    - GeoCSV_2_netCDF.py - a Python script to read a 2D or 3D GeoCSV Earth model file in geographic or projected 
      coordinate system and display its header information or convert it to netCDF format.

 * For information on EMC (IRIS Earth Model Collaboration) visit: http://ds.iris.edu/ds/products/emc/
 * For information on netCDF (Network Common Data Form) visit: https://www.unidata.ucar.edu/software/netcdf/
 * For information on GeoCSV (tabular text formatting for geoscience data) visit: 
       http://geows.ds.iris.edu/documents/GeoCSV.pdf

 **USAGE:**
You can work with each script directly and independently. To view the usage message, run the script with the 
help (**-h**) option:

* netCDF_2_GeoCSV.py -h
* GeoCSV_2_netCDF.py -h

**How to run netCDF\_2\_GeoCSV.py:**

    - use the sample 3D netCDF file under the "samples" directory or download a 3D netCDF Earth model file from the EMC 
      repository:
   http://ds.iris.edu/ds/products/emc-earthmodels/
   
    - [optional] create a "data" directory with a path relative to the script as "../data/" and place the netCDF file 
      there
    - run the script to convert the netCDF file. The converted file will have the same name as the netCDF file but with 
      extension ".csv" and will be placed under the same directory where the netCDF file was found.
            
     netCDF_2_GeoCSV.py -i FILE -x longitude -y latitude -z depth -m mode -d -H
        -i [required] FILE is the input GeoCSV Earth model file
        -x [required] the netCDF variable representing the x-axis (must match the netCDF file's x variable)
        -y [required] the netCDF variable representing the y-axis (must match the netCDF file's y variable)
        -z [required] the netCDF variable representing the z-axis (must match the netCDF file's z variable OR set to 'none' for 2D models)
        -m [default:single] output mode [single: single GeoCSV file OR depth: one GeoCSV file per depth
        -d [default:False] debug mode
        -H [default:False] display headers only
                 
    NOTE: only a .csv files converted with "-m single" option can be converted back to netCDF using the 
          GeoCSV_2_netCDF_3D.py script

**How to run GeoCSV\_2\_netCDF.py:**

    - use the sample GeoCSV file under the "samples" directory as a template to create your own GeoCSV file with proper 
      header, or convert a netCDF model file to GeoCSV using the netCDF_2_GeoCSV.py script above
    - run the script to convert the GeoCSV file to netCDF. The converted file will have the same name as the GeoCSV file
      but with extension ".nc" and will be placed under the same directory where the GeoCSV file was found.

     GeoCSV_2_netCDF.py -i FILE -d -H
        -i [required] FILE is the input GeoCSV Earth model file
        -d [default:False] debug mode
        -H [default:False] display headers only

    NOTE: only a single ".csv" files can be converted back to netCDF using this script

            
 **HEADERS:**
 
    GeoCSV Headers are similar to the netCDF headers with the following exceptions:
        - GeoCSV-specific header lines:
            # dataset: GeoCSV2.0
            # created: 2022-08-18 13:33:30 UTC (netCDF_2_GeoCSV.py)
            # netCDF_file: SAW642ANB_kmps.nc
            # delimiter: |
            # {VAR}_column: {VAR} where {VAR} represents a model or coordinate variable such as latitude, vs, etc.
            
         - The global attributes in netCDF header start with a ":", like    :geospatial_lat_min = "-90.00" ; 
           while the same global attribute under GeoCSV starts with global_, like # global_geospatial_lat_min: -90.00

    EMC Earth models netCDF header follow the EMC netCDF 3 guidelines see: 
   http://ds.iris.edu/ds/products/emc-contributionguide/

    To view the file header information for both netCDF and GeoCSV, run the above scripts with -H option:
            netCDF_2_GeoCSV.py -H -i SAW642ANB_kmps.nc
            GeoCSV_2_netCDF.py -H -i SAW642ANB_kmps.csv
            netCDF_2_GeoCSV.py -H -i Crustal_Thickness_Error.nc
            GeoCSV_2_netCDF.py -H -i Crustal_Thickness_Error.csv

**Creating EMC GeoCSV Earth Model Files**

   The following highlights the steps you need to follow in order to create an EMC Earth model file in GeoCSV format:
       
       Note: The GeoCSV Earth model files in general should follow "EMC's guidelines" at:
   http://ds.iris.edu/ds/products/emc-contributionguide/
   
       Note: Familiarize yourself with the "GeoCSV format" at: 
   http://geows.ds.iris.edu/documents/GeoCSV.pdf
   
      Note: The following steps are for converting 3D models. For 2D models use the corresponding scripts but note that
            in 2D domain, each point is defined by a set of (latitude,longitude) values. 
            Your data must be sorted in the same order (latitude column first and then the longitude column).
       
     - download the "SAW642ANB_kmps.csv" from the samples folder and inspect it. The file has three parts:
          * metadata header - starts at the beginning of the file with each line starting with a "#" 
          * data header - a single line following metadata header that identifies columns
          * data - lines following the header
          
     - download the header.csv file from the samples directory as a header template and edit it to include information 
       about your model (do not change the first line: # dataset: GeoCSV2.0). This file contains the metadata and data 
       headers of the "SAW642ANB_kmps.csv" file. Also edit the data header (last line with no starting # to include all 
       variables that you have defined under the mtadata header.
     - create the data section by writing data to a file using delimiter identified in the metadata header section and 
       in the same order as the defined in the data header line.
       **NOTE:** For the EMC models with a 3D grid, each point is defined by a set of (depth,latitude,longitude) values. 
                 Your data must be sorted in the same order (depth column first, latitude column second and the 
                 longitude column third).
     - add the data section under the "data header line" of the header file
     - test your GeoCSV file by extracting its header using the "GeoCSV_2_netCDF_3D.py" script:
               GeoCSV_2_netCDF_3D.py -H -i {YOUR_FILE_NAME}.csv
     - convert your file to netCDF format using the "GeoCSV_2_netCDF_3D.py" script:
               GeoCSV_2_netCDF_3D.py -i {YOUR_FILE_NAME}.csv
          
     
**HEADER EXAMPLE:**

netCDF_2_GeoCSV.py -i KEA20-Moho.r0.0.nc -H

[INFO] netCDF_2_GeoCSV.py version v2022.237
[INFO] working on the 3D netCDF file: KEA20-Moho.r0.0.nc

[INFO] Input netCDF File: KEA20-Moho.r0.0.nc


**netCDF header information:**


	dimensions:
		latitude 181
		longitude 401

	variables:
		moho:
			variable = moho
			long_name = Moho depth relative to a mean Earth radius of 6371 km
			display_name = Moho depth (km)
			units = km

	global attributes:
			title = Moho depth model for East Asia
			id = KEA20_Moho
			data_revision = r.0.0
			Conventions = CF-1.0
			Metadata_Conventions = Unidata Dataset Discovery v1.0
			keywords = seismic, surface wave, tomography, East Asia, Moho, receiver functions, reflection
			acknowledgment = Model was provided by Michael Witek, Department of Geophysics, Kangwon National University
			history = 2022-08-22 19:48:38 UTC Converted to netCDF by GeoCSV_2_netCDF.py v2022.230 from KEA20-Moho.r0.0.csv; 2022-08-22 19:46:46 UTC Converted to GeoCSV by netCDF_2_GeoCSV.py ,v2022.230 from KEA20-Moho.r0.0.nc; v1.0.0 created 2021-06-02
			summary = Here we provide depths to the Moho boundary in KEA20, relative to a mean Earth radius of 6371 km.; These Moho depths result from an initial inversion of point constraints from teleseismic; receiver function and seismic reflection studies, plus additional Moho perturbations from the; surface wave dispersion data inversion. The reference model Moho is taken from CRUST1.0.  ; For additional information please see the supplementary information of Witek et al. (2021).
			reference = ; @article{witek-2021-kea20,; title={Radial anisotropy in East Asia from multimode surface wave tomography},; author={M. Witek and S.-J. Chang and D.Y. Lim and S. Ning and J. Ning},; journal={Journal of Geophysical Research},; year={2021},; volume={},; pages={}; };
			reference_pid = 
			author_name = M. Witek
			author_email = mwitek2@gmail.com
			author_institution = Department of Geophysics, Kangwon National University
			geospatial_lat_min = 15.0
			geospatial_lat_max = 51.0
			geospatial_lat_resolution = 0.20000000298023224
			geospatial_lat_units = degrees_north
			geospatial_lon_min = 70.0
			geospatial_lon_max = 150.0
			geospatial_lon_resolution = 0.20000000298023224
			geospatial_lon_units = degrees_east
			netcdf_file = KEA20_Moho.nc
			repository_name = EMC
			repository_institution = IRIS EMC
			repository_pid = doi:10.17611/dp/emc.2021.kea20.1


**GeoCSV header information:**

    # dataset: GeoCSV2.0
    # created: 2022-08-24 16:09:10 UTC (netCDF_2_GeoCSV.py)
    # netCDF_file: KEA20-Moho.r0.0.nc
    # delimiter: |
    # global_title: Moho depth model for East Asia
    # global_id: KEA20_Moho
    # global_data_revision: r.0.0
    # global_Conventions: CF-1.0
    # global_Metadata_Conventions: Unidata Dataset Discovery v1.0
    # global_keywords: seismic, surface wave, tomography, East Asia, Moho, receiver functions, reflection
    # global_acknowledgment: Model was provided by Michael Witek, Department of Geophysics, Kangwon National University
    # global_history: 2022-08-24 16:09:10 UTC Converted to GeoCSV by netCDF_2_GeoCSV.py ,v2022.237 from KEA20-Moho.r0.0.nc; 2022-08-22 19:48:38 UTC Converted to netCDF by GeoCSV_2_netCDF.py v2022.230 from KEA20-Moho.r0.0.csv; 2022-08-22 19:46:46 UTC Converted to GeoCSV by netCDF_2_GeoCSV.py ,v2022.230 from KEA20-Moho.r0.0.nc; v1.0.0 created 2021-06-02
    # global_summary: Here we provide depths to the Moho boundary in KEA20, relative to a mean Earth radius of 6371 km.; These Moho depths result from an initial inversion of point constraints from teleseismic; receiver function and seismic reflection studies, plus additional Moho perturbations from the; surface wave dispersion data inversion. The reference model Moho is taken from CRUST1.0.  ; For additional information please see the supplementary information of Witek et al. (2021).
    # global_reference: ; @article{witek-2021-kea20,; title={Radial anisotropy in East Asia from multimode surface wave tomography},; author={M. Witek and S.-J. Chang and D.Y. Lim and S. Ning and J. Ning},; journal={Journal of Geophysical Research},; year={2021},; volume={},; pages={}; };
    # global_reference_pid: 
    # global_author_name: M. Witek
    # global_author_email: mwitek2@gmail.com
    # global_author_institution: Department of Geophysics, Kangwon National University
    # global_geospatial_lat_min: 15.0
    # global_geospatial_lat_max: 51.0
    # global_geospatial_lat_resolution: 0.20000000298023224
    # global_geospatial_lat_units: degrees_north
    # global_geospatial_lon_min: 70.0
    # global_geospatial_lon_max: 150.0
    # global_geospatial_lon_resolution: 0.20000000298023224
    # global_geospatial_lon_units: degrees_east
    # global_netcdf_file: KEA20_Moho.nc
    # global_repository_name: EMC
    # global_repository_institution: IRIS EMC
    # global_repository_pid: doi:10.17611/dp/emc.2021.kea20.1


NOTE: For information on variables, the -x, -y, and -z parameters are required			


**Models with projected coordinate systems***

The extended EMC netCDF format includes both the projected coordinate system variables (x and y) and the geographic 
latitude and longitude variables that are defined as a function of the primary x and y coordinate variables 
(Table 1, right). These two-dimensional latitude and longitude coordinate variables allow users to quickly switch 
between the model’s projected coordinates and the geographic latitudes and longitudes. In this extended netCDF format, 
the x and y of the projection coordinate system become the primary coordinates, and the geographic latitude and 
longitude, which depend on x and y, represent the auxiliary coordinate variables. The directions of the primary 
coordinate variables are defined by adding the **axis** attribute to the x and y coordinate variables. The model 
variables (thickness and vp in Table 1, right) are now defined based on the primary x and y coordinate variables, 
and are tied to the auxiliary coordinate variables via the **coordinates** attribute.

* for more information see: https://ds.iris.edu/ds/newsletter/vol99/no1/543/emc-two-dimensional-latitude-longitude-coordinate-variables/

*Table 1*. A comparison of dimensions and variables in a classic EMC model with geographic coordinates (top), 
and in an extended EMC model with a projected coordinate system (bottom).

* Dimensions and variables in a classic EMC model

**dimensions**:

    depth = 5;
    latitude = 10;
    longitude = 20;
**variables**:

    double longitude(longitude);
        .
        .
    double latitude(latitude);
        .
        .
    double depth(depth);
        .
        .
    double thickness(latitude, longitude);
        .
        .
    double vp(depth, latitude, longitude);
        .
* Dimensions and variables in an extended EMC model with a projected coordinate system 

**dimensions**:

    depth = 5;
    y = 10;
    x = 20;
**variables**:

    float y(y);
        y:axis = "Y";
        .
        .
    float x(x);
        x:axis = "X"
        .
        .
    double depth(depth);
        .
        .
    double latitude(y, x);
        .
        .
    double longitude(y, x);
        .
        .
    double thickness(y, x);
        .
        thickness:coordinates = "longitude latitude";
      .
    double vp(depth, y, x);
        .
        vp:coordinates = "longitude latitude";