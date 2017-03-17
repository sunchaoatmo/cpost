import os
import sys
from netCDF4 import Dataset
#output_folder="//homes/sunchao/lustre//Netcdf/"
r95tfile="//home-4/sunchao@umd.edu/data/OBS_R95T_hist_seasonal.nc"
if os.path.isfile(r95tfile):
  f=Dataset(r95tfile)
  r95t_hist=f.variables["R95T_hist"]
else:
  sys.exit("please provide the file to support r95 calculation")
r95t_hist=Dataset(r95tfile, 'r')
