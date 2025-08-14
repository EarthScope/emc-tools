SRC_DIR = "src"
WORK_DIR = "work"
LIB_DIR = "lib"
dataset_version = "GeoCSV 2.0"
delimiter = "|"
na_rep = "nan"
extension = {
    "geocsv": ".csv",
    "ggeocsv": ".gcsv",
    "csv": ".csv",
    "netcdf": ".nc",
    "json": ".json",
    "nonlinloc": ".txt",
}
netcdf_engines = {".nc": "netcdf4"}
line_break = "\n"
cmap = "jet_r"
grid_spatial_ref = ["latitude_longitude", "transverse_mercator"]

valid_variable_list = [
    "Vs",
    "Vp",
    "Vs_perturbation",
    "Vp_perturbation",
    "Density",
    "Rayleigh_Phase",
    "Love_Phase",
    "Rayleigh_Group",
    "Love_Group",
    "Depth",
]
