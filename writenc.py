from netCDF4 import Dataset
def createsmnc(casename,vname,periods,nx,ny,dimension=3):
  filename            = "%s_%s_%s.nc"%(casename,vname,periods)
  rootgrp             = Dataset(filename                   , "w")
  rootgrp.createDimension("lat"      , nx )
  rootgrp.createDimension("lon"      , ny )
  nt=4 if periods=="seasonal" else 12
  rootgrp.createDimension(periods   , nt)  #hard coded but usually we consider only 4 seasons
  subtime              = rootgrp.createVariable(periods    , "i4"  , (periods , ))
  subtime[:]           = range(1,nt+1)
  subtime.units        = periods
  if dimension==3:
    rootgrp.createDimension("time"     , None)
    times               = rootgrp.createVariable("time"      , "i4"  , ("time"   , ))
    times.units         = "year"
    nc_var              = rootgrp.createVariable(vname,"f4",("time",periods,"lat","lon",))
  else:
    nc_var              = rootgrp.createVariable(vname,"f4",(periods,"lat","lon",))
  return rootgrp


def createnc(casename,vname,periods,units,calendar,nx,ny,dimension=3):
  filename            = "%s_%s_%s.nc"%(casename,vname,periods)
  rootgrp             = Dataset(filename                   , "w")
  rootgrp.createDimension("lat"      , nx )
  rootgrp.createDimension("lon"      , ny )
  rootgrp.createDimension("time"     , None)
  T                   = rootgrp.createVariable("time"    , "f4"  , ("time" , ))
  T.units             = units
  T.calendar          = calendar
  nc_var              = rootgrp.createVariable(vname,"f4",("time","lat","lon",))
  return rootgrp 


