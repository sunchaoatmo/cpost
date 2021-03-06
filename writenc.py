from netCDF4 import Dataset
def createsmnc(casename,vname,periods,fields,nx,ny,nz=None):
  filename            = "%s_%s_%s.nc"%(casename,vname,periods)
  rootgrp             = Dataset(filename                   , "w")
  rootgrp.createDimension("lat"      , nx )
  rootgrp.createDimension("lon"      , ny )
  nt=4 if periods=="seasonal" else 12
  rootgrp.createDimension(periods   , nt)  #hard coded but usually we consider only 4 seasons
  subtime              = rootgrp.createVariable(periods    , "i4"  , (periods , ))
  subtime[:]           = range(1,nt+1)
  subtime.units        = periods
  rootgrp.createDimension("time"     , None)
  times               = rootgrp.createVariable("time"      , "i4"  , ("time"   , ))
  times.units         = "year"
  if nz:
    nc_dim=("time",periods,"bottom_top","lat","lon",)
    rootgrp.createDimension("bottom_top"      , nz )
  else:
    nc_dim=("time",periods,"lat","lon",)
  for field in fields:
    rootgrp.createVariable(field,"f4",nc_dim)

  nc_dim=("lat","lon",)
  lat=rootgrp.createVariable("lat","f4",nc_dim,fill_value=-999)
  lat.units ="degrees_north"
  lat.long_name ="latitude"
  lat.standard_name ="latitude"

  lon=rootgrp.createVariable("lon","f4",nc_dim,fill_value=-999)
  lon.units ="degrees_east"
  lon.long_name ="longitude"
  lon.standard_name ="longitude"

  return rootgrp


#def createnc(casename,vname,periods,units,calendar,fields,nx,ny,nz=None):
def createnc(filename,vname,periods,units,calendar,fields,nx,ny,nz=None):
  #filename            = "%s_%s_%s.nc"%(casename,vname,periods)
  rootgrp             = Dataset(filename                   , "w")
  rootgrp.createDimension("west_east"        , nx )
  rootgrp.createDimension("south_north"      , ny )
  rootgrp.createDimension("time"     , None)
  T                   = rootgrp.createVariable("time"    , "f8"  , ("time" , ))
  T.units             = units
  T.calendar          = calendar
  T.axis              = "T"
  if nz:
    rootgrp.createDimension("bottom_top"      , nz )
    nc_dim=("time","bottom_top","south_north","west_east",)
  else:
    nc_dim=("time","south_north","west_east",)
  for field in fields:
    rootgrp.createVariable(field,"f4",nc_dim,fill_value=-999)

  nc_dim=("south_north","west_east",)
  lat=rootgrp.createVariable("lat","f4",nc_dim,fill_value=-999)
  lat.units ="degrees_north"
  lat.long_name ="latitude"
  lat.standard_name ="latitude"

  lon=rootgrp.createVariable("lon","f4",nc_dim,fill_value=-999)
  lon.units ="degrees_east"
  lon.long_name ="longitude"
  lon.standard_name ="longitude"

  return rootgrp 


