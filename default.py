import os
import sys
from netCDF4 import Dataset
import numpy as np
#output_folder="//homes/sunchao/lustre//Netcdf/"
wrfinputfilename="///homes/sunchao/lustre/US/icbc_1999/wrfinput_d01.1998120100"
r95tfile="/homes/sunchao/data//OBS_R95T_hist_seasonal.nc"
if os.path.isfile(r95tfile):
  f=Dataset(r95tfile)
  r95t_hist=f.variables["R95T_hist"]
else:
  sys.exit("please provide the file to support r95 calculation")
r95t_hist=Dataset(r95tfile, 'r')
wrfinput=Dataset(wrfinputfilename, 'r')
# be awared the z_levs should be sorted from large to small
z_levs=[1000.0, 975.0, 950.0, 925.0, 900.0, 875.0, 850.0, 825.0, 800.0, 775.0, 750.0, 700.0, 650.0, 600.0, 550.0, 500.0, 450.0, 400.0, 350.0, 300.0, 250.0, 225.0, 200.0, 175.0, 150.0, 125.0, 100.0, 70.0]
z_levs=np.array(z_levs)
z_levs=100*z_levs
number_of_zlevs  =len(z_levs)
