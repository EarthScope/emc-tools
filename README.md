 
Incorporated Research Institutions for Seismology (IRIS):

    Data Management Center (DMC)
    Data Products Team
    IRIS Earth Model Collaboration (EMC) - Tools

COMMENTS/QUESTIONS:

    Please contact manoch@iris.washington.edu


 2018-10-25
----------------------------------------------------------------------------------------------------------------------------------------------------------------------

 DESCRIPTION:

 EMC Tools contains a set of Python scripts for converting the EMC's netCDF Earth model files (in netCDF 3 format) to 
 and from GeoCSV format:
 
    - netCDF_2_GeoCSV_3D.py - a Python script to read a 3D netCDF Earth model file and display its header information or
      convert it to GeoCSV format.
    - GeoCSV_2_netCDF_3D.py - a Python script to read a 3D GeoCSV Earth model file and display its header information or
      convert it to netCDF format.
    - netCDF_2_GeoCSV_2D.py - a Python script to read a 2D netCDF Earth model file and convert it to GeoCSV format.
    - GeoCSV_2_netCDF_2D.py - a Python script to read a 2D GeoCSV Earth model file and convert it to netCDF format.

 * For information on EMC (IRIS Earth Model Collaboration) visit: http://ds.iris.edu/ds/products/emc/
 * For information on netCDF (Network Common Data Form) visit: https://www.unidata.ucar.edu/software/netcdf/
 * For information on GeoCSV (tabular text formatting for geoscience data) visit: 
       http://geows.ds.iris.edu/documents/GeoCSV.pdf

 **USAGE:**
You can work with each script directly and independently. To view the usage message, run the script with the 
help (**-h**) option:

* netCDF_2_GeoCSV_3D.py -h
* GeoCSV_2_netCDF_3D.py -h
* netCDF_2_GeoCSV_2D.py -h
* GeoCSV_2_netCDF_2D.py -h

**How to run netCDF\_2\_GeoCSV\_3D.py:**

    - use the sample 3D netCDF file under the "samples" directory or download a 3D netCDF Earth model file from the EMC 
      repository:
   http://ds.iris.edu/ds/products/emc-earthmodels/
   
    - [optional] create a "data" directory with a path relative to the script as "../data/" and place the netCDF file 
      there
    - run the script to convert the netCDF file. The converted file will have the same name as the netCDF file but with 
      extension ".csv" and will be placed under the same directory where the netCDF file was found.
            
            netCDF_2_GeoCSV_3D.py -i SAW642ANB_kmps.nc (will look for the file under the script directory and data 
            directory)

            The following options are also available:
                -d debug mode with verbose execution
                -m [single, depth] the mode option tells the script how to convert the netCDF file
                    single (default) will create a single ".csv" file
                    depth will create one ".csv" file per depth
                -H display file header information for both netCDF and GeoCSV (no output file created)
                 
            NOTE: only a .csv files converted with "-m single" option can be converted back to netCDF using the 
                  GeoCSV_2_netCDF_3D.py script

**How to run GeoCSV\_2\_netCDF\_3D.py:**

    - use the sample GeoCSV file under the "samples" directory as a template to create your own GeoCSV file with proper 
      header, or convert a netCDF model file to GeoCSV using the netCDF_2_GeoCSV_3D.py script above
    - run the script to convert the GeoCSV file to netCDF. The converted file will have the same name as the GeoCSV file
      but with extension ".nc" and will be placed under the same directory where the GeoCSV file was found.

            GeoCSV_2_netCDF_3D.py -i SAW642ANB_kmps.csv (will look for the file under the script directory and data 
            directory)

            The following options are also available:
                -d debug mode with verbose execution
                -H display file header information for both netCDF and GeoCSV 

            NOTE: only a single ".csv" files can be converted back to netCDF using this script
        
**How to run netCDF\_2\_GeoCSV\_2D.py:**

    - use the sample 2D netCDF file under the "samples" directory or download a 2D netCDF Earth model file from the EMC 
      repository:
   http://ds.iris.edu/ds/products/emc-earthmodels/
   
    - [optional] create a "data" directory with a path relative to the script as "../data/" and place the netCDF file 
      there
    - run the script to convert the netCDF file. The converted file will have the same name as the netCDF file but with 
      extension ".csv" and will be placed under the same directory where the netCDF file was found.
            
            netCDF_2_GeoCSV_2D.py -i Crustal_Thickness_Error.nc (will look for the file under the script directory and 
            data directory)

            The following options are also available:
                -d debug mode with verbose execution
                -H display file header information for both netCDF and GeoCSV (no output file created)

**How to run GeoCSV\_2\_netCDF\_2D.py:**

    - use the sample GeoCSV file under the "samples" directory as a template to create your own GeoCSV file with proper 
      header, or convert a netCDF model file to GeoCSV using the netCDF_2_GeoCSV_2D.py script above
    - run the script to convert the GeoCSV file to netCDF. The converted file will have the same name as the GeoCSV file
      but with extension ".nc" and will be placed under the same directory where the GeoCSV file was found.

            GeoCSV_2_netCDF_2D.py -i Crustal_Thickness_Error.csv (will look for the file under the script directory and 
            data directory)

            The following options are also available:
                -d debug mode with verbose execution
                -H display file header information for both netCDF and GeoCSV 
            
 **HEADERS:**
 
    GeoCSV Headers are similar to the netCDF headers with the following exceptions:
        - GeoCSV-specific header lines:
            # dataset: GeoCSV2.0
            # created: 2018-10-18 13:33:30 UTC (netCDF_2_GeoCSV_3D.py)
            # netCDF_file: SAW642ANB_kmps.nc
            # delimiter: |
            # {VAR}_column: {VAR} where {VAR} represents a model or coordinate variable such as latitude, vs, etc.
            
         - The global attributes in netCDF header start with a ":", like    :geospatial_lat_min = "-90.00" ; 
           while the same global attribute under GeoCSV starts with global_, like # global_geospatial_lat_min: -90.00

    EMC Earth models netCDF header follow the EMC netCDF 3 guidelines see: 
   http://ds.iris.edu/ds/products/emc-contributionguide/

    To view the file header information for both netCDF and GeoCSV, run the above scripts with -H option:
            netCDF_2_GeoCSV_3D.py -H -i SAW642ANB_kmps.nc
            GeoCSV_2_netCDF_3D.py -H -i SAW642ANB_kmps.csv
            netCDF_2_GeoCSV_2D.py -H -i Crustal_Thickness_Error.nc
            GeoCSV_2_netCDF_2D.py -H -i Crustal_Thickness_Error.csv

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

netCDF_2_GeoCSV_3D.py -H -i SAW642ANB_kmps.nc


[INFO] netCDF_2_GeoCSV_3D.py version R.0.0.2018.295
[INFO] Input netCDF File: SAW642ANB_kmps.nc


**netCDF header information:**


	dimensions:
		depth 46
		latitude 91
		longitude 180

	variables:
		vs:
			long_name = Shear Velocity
			display_name = S Velocity (km/s)
			units = km.s-1
			missing_value = 99999.0
			_FillValue = 99999.0
		vp:
			long_name = P Velocity
			display_name = P Velocity (km/s)
			units = km.s-1
			missing_value = 99999.0
			_FillValue = 99999.0
		rho:
			long_name = Density
			display_name = Density (kg/m^3)
			units = kg.m-3
			missing_value = 99999.0
			_FillValue = 99999.0
		Qs:
			long_name = Shear Q
			display_name = Shear Q (dimensionless)
			units = count
			missing_value = 99999.0
			_FillValue = 99999.0

	global attributes:
			title = A radially anisotropic whole mantle model with improved crustal corrections
			id = SAW642ANb_kmps
			summary = SAW642ANb is a radially anisotropic shear velocity model, parameterized in terms of isotropic S velocity (Voigt average) and the anisotropic parameter, xi (Vsh^2/Vsv^2). 
			reference = Panning, Lekic and Romanowicz (2010)
			references = http://ds.iris.edu/ds/products/emc-references
			keywords = seismic, tomography, shear wave, s wave, elastic waveform 
			Conventions = CF-1.0
			Metadata_Conventions = Unidata Dataset Discovery v1.0
			creator_name = IRIS EMC
			creator_url = http://ds.iris.edu/ds/products/emc/
			creator_email = product@iris.edu
			institution = IRIS DMC
			acknowledgment = Model was provided by Professor Mark Panning,  Department of Geological SciencesUniversity of Florida 
			history = 2011-06-21 created 2016-08-11 IRID DMC, updated variable definitions to make sure they are compatible with CF
			comment = model converted to netCDF by IRIS EMC
			geospatial_lat_min = -90.00
			geospatial_lat_max = 90.00
			geospatial_lat_units = degrees_north
			geospatial_lat_resolution = 2
			geospatial_lon_min = -180.00
			geospatial_lon_max = 178.0
			geospatial_lon_units = degrees_east
			geospatial_lon_resolution = 2
			geospatial_vertical_min = 50
			geospatial_vertical_max = 2850
			geospatial_vertical_units = km
			geospatial_vertical_positive = down
			
**GeoCSV header information:**

	    # dataset: GeoCSV2.0
        # created: 2018-10-23 21:23:32 UTC (netCDF_2_GeoCSV_3D.py)
        # netCDF_file: SAW642ANB_kmps.nc
        # delimiter: |
        # global_title: A radially anisotropic whole mantle model with improved crustal corrections
        # global_id: SAW642ANb_kmps
        # global_summary: SAW642ANb is a radially anisotropic shear velocity model, parameterized in terms of isotropic S velocity (Voigt average) and the anisotropic parameter, xi (Vsh^2/Vsv^2).
        # global_reference: Panning, Lekic and Romanowicz (2010)
        # global_references: http://ds.iris.edu/ds/products/emc-references
        # global_keywords: seismic, tomography, shear wave, s wave, elastic waveform
        # global_Conventions: CF-1.0
        # global_Metadata_Conventions: Unidata Dataset Discovery v1.0
        # global_creator_name: IRIS EMC
        # global_creator_url: http://ds.iris.edu/ds/products/emc/
        # global_creator_email: product@iris.edu
        # global_institution: IRIS DMC
        # global_acknowledgment: Model was provided by Professor Mark Panning,  Department of Geological SciencesUniversity of Florida
        # global_history: 2011-06-21 created 2016-08-11 IRID DMC, updated variable definitions to make sure they are compatible with CF
        # global_comment: model converted to netCDF by IRIS EMC
        # global_geospatial_lat_min: -90.00
        # global_geospatial_lat_max: 90.00
        # global_geospatial_lat_units: degrees_north
        # global_geospatial_lat_resolution: 2
        # global_geospatial_lon_min: -180.00
        # global_geospatial_lon_max: 178.0
        # global_geospatial_lon_units: degrees_east
        # global_geospatial_lon_resolution: 2
        # global_geospatial_vertical_min: 50
        # global_geospatial_vertical_max: 2850
        # global_geospatial_vertical_units: km
        # global_geospatial_vertical_positive: down
        # latitude_column: latitude
        # latitude_long_name: Latitude; positive north
        # latitude_units: degrees_north
        # latitude_standard_name: latitude
        # longitude_column: longitude
        # longitude_long_name: Longitude; positive east
        # longitude_units: degrees_east
        # longitude_standard_name: longitude
        # depth_column: depth
        # depth_long_name: depth below earth surface
        # depth_units: km
        # depth_positive: down
        # vs_column: vs
        # vs_long_name: Shear Velocity
        # vs_display_name: S Velocity (km/s)
        # vs_units: km.s-1
        # vs_missing_value: 99999.0
        # vs__FillValue: 99999.0
        # vp_column: vp
        # vp_long_name: P Velocity
        # vp_display_name: P Velocity (km/s)
        # vp_units: km.s-1
        # vp_missing_value: 99999.0
        # vp__FillValue: 99999.0
        # rho_column: rho
        # rho_long_name: Density
        # rho_display_name: Density (kg/m^3)
        # rho_units: kg.m-3
        # rho_missing_value: 99999.0
        # rho__FillValue: 99999.0
        # Qs_column: Qs
        # Qs_long_name: Shear Q
        # Qs_display_name: Shear Q (dimensionless)
        # Qs_units: count
        # Qs_missing_value: 99999.0
        # Qs__FillValue: 99999.0
			





