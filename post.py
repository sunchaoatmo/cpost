#!/usr/bin/env python
import numpy.ma as ma
import numpy as np
import calendar
from netCDF4 import Dataset
from netCDF4 import date2num,num2date
from datetime import datetime,timedelta
from constant import *   # only several constants
from POSTparameter import var_parameters
from argument import args
from cs_stat import cs_stat
import glob
import os.path
import default

# This code assumes each wrfout file contains one full day of output 
# and that all data is consecutive (no gaps)


def selectseasonaldata(yb_select,ye_select,data_beg_date,units_cur,calendar_cur,seasonnames):
# This is a function can select each season's index for standard daily data
# input yb_select the begin year you want to select the date
# input ye_select the end   year you want to select the date
# input date_beg_date the first day of your date
# output season_periods_loc is a 2d list each containing index for each season in a[0][:] a[1][:] ...
  from collections import defaultdict
  seasonList=[ [12, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]
  season_periods_loc=defaultdict(list) #this is for mdata, because in temp, it has already cutted onece and the relative starting point changes to 0 now
  first_day_num=date2num(data_beg_date,units=units_cur,calendar=calendar_cur)
  for iseason,seasonname in enumerate(seasonnames):
    for year in range(yb_select,ye_select+1):
      if iseason==0 and not year==yb_select:
        cur_y=year-1
      else:
        cur_y=year
      m=seasonList[iseason][0]
      d=1 
      cur_beg_date=datetime(cur_y,m,d,0,0,0)

      cur_y= year
      m=seasonList[iseason][-1]
      d=calendar.monthrange(cur_y,m)[1]
      cur_end_date=datetime(cur_y,m,d,0,0,0)
      beg_index=int(date2num(cur_beg_date,units=units_cur,calendar=calendar_cur)-first_day_num)
      end_index=int(date2num(cur_end_date,units=units_cur,calendar=calendar_cur)-first_day_num)
      season_periods_loc[seasonname].extend(range(beg_index,end_index+1))
  return season_periods_loc 

def createseasonalnc(casename,vname,units,nx,ny,dimension=3):
  filename            = "%s_%s_seasonal.nc"%(casename,vname)
  rootgrp             = Dataset(filename                   , "w")
  lat                 = rootgrp.createDimension("lat"      , nx )
  lon                 = rootgrp.createDimension("lon"      , ny)
  season              = rootgrp.createDimension("season"   , 4)  #hard coded but usually we consider only 4 seasons
  season              = rootgrp.createVariable("season"    , "i4"  , ("season" , ))
  season.units        = "season"
  season[:]           = range(1,5)
  if dimension==3:
    time                = rootgrp.createDimension("time"     , None)
    times               = rootgrp.createVariable("time"      , "i4"  , ("time"   , ))
    times.units         = "year"
    nc_var              = rootgrp.createVariable(vname,"f4",("time","season","lat","lon",))
    return rootgrp
  else:
    nc_var              = rootgrp.createVariable(vname,"f4",("season","lat","lon",))
    nc_var.units        = units

def daysbetween(by,bm,bd,ey,em,ed):
   tb=date(by,bm,bd)
   te=date(ey,em,ed)
   td=abs(te-tb).days
   return td;


def createnc(casename,vname,periods,units,calendar,nx,ny,dimension=3):
  filename            = "%s_%s_%s.nc"%(casename,vname,periods)
  rootgrp             = Dataset(filename                   , "w")
  lat                 = rootgrp.createDimension("lat"      , nx )
  lon                 = rootgrp.createDimension("lon"      , ny)
  time                = rootgrp.createDimension("time"     , None)
  T                   = rootgrp.createVariable("time"    , "f4"  , ("time" , ))
  T.units             = units
  T.calendar          = "gregorian"
  nc_var              = rootgrp.createVariable(vname,"f4",("time","lat","lon",))
  return rootgrp  #,nc_var,clat0,clon0

def wrftimetodate(wrfstr):
    YEAR=slice(0,4)
    MONTH=slice(5,7)  
    DAY=slice(8,10)  
    HOUR=slice(11,13)
    MINUTE=slice(14,16)
    SECOND=slice(17,19)
    Year=int(''.join(list(wrfstr[YEAR])))
    Month=int(''.join(list(wrfstr[MONTH])))
    Day=int(''.join(list(wrfstr[DAY])))
    Hour=int(''.join(list(wrfstr[HOUR])))
    Minute=int(''.join(list(wrfstr[MINUTE])))
    Second=int(''.join(list(wrfstr[SECOND])))
    wrfdate=datetime(Year, Month, Day, Hour, Minute, Second)
    return wrfdate

calendar_cur=args.calendar
periods=args.p[0]
casename=args.n
vnames=args.v

if periods=="daily":
  units_cur = 'days since 0001-01-01 00:00'
elif periods=="hourly":
  units_cur = 'hours since 0001-01-01 00:00'
Oneday=timedelta(days=1)
filenames_all=sorted(glob.glob(prefix))
filesize0=os.path.getsize(filenames_all[0])
filenames=[]
for filename in filenames_all:
  if not filesize0==os.path.getsize(filename):
    break
  filenames.append(filename)

for vname in vnames:
  shiftday=var_parameters[vname]['shiftday'] 
  compute_mode=var_parameters[vname]['compute_mode'] 
  rawfname="%s_%s_%s.nc"%(casename,vname,periods)
  ncexist=os.path.isfile(rawfname)
  lastindex=0
  if ncexist:
    rawnc=Dataset(rawfname,'a')
    lastday=num2date(rawnc.variables["time"][-1],units=units_cur,calendar=calendar_cur)
    lastwrfout="wrfout_d01_%s"%lastday.strftime(wrfout_data_fmt)
    try:
      lastindex=filenames.index(lastwrfout)
      del filenames[:lastindex]
    except:
      print("STOP! There is a GAP between the record of the last day %s and earlieast wrfout we have in this folder"%(lastday))
      import sys
      sys.exit()
  ncfile_last=Dataset(filenames[0],'r')
  var_shape=ncfile_last.variables[vname].shape
  outputdim=len(var_shape)
  if outputdim==3:
    (nstep,nx,ny)=ncfile_last.variables[vname].shape
  elif outputdim==4:
    (nstep,nlev,nx,ny)=ncfile_last.variables[vname].shape
  outputdata=np.empty([len(filenames)-shiftday,nx,ny])
  outputtime=np.empty([len(filenames)-shiftday])
  print(filenames)
  simbeg_date=wrftimetodate(Dataset(filenames[shiftday],'r').variables['Times'][0])

  for iday,filename in enumerate(filenames[shiftday:]):
    ncfile_cur=Dataset(filename,'r')
    curtime=ncfile_cur.variables['Times']
    date_curstep=wrftimetodate(curtime[0])
    if date_curstep==iday*Oneday+simbeg_date:
      if len(curtime)==nstep:
        if compute_mode==6:
          if vname=="PRAVG":
            outputdata[iday,:,:]=(ncfile_cur.variables['RAINC'][0,:,:]-ncfile_last.variables['RAINC'][0,:,:]
                                 +ncfile_cur.variables['RAINNC'][0,:,:]-ncfile_last.variables['RAINNC'][0,:,:])
          else:
            outputdata[iday,:,:]=ncfile_cur.variables[vname][0,:,:]-ncfile_last.variables[vname][0,:,:]
        elif compute_mode==1:
          outputdata[iday,:,:]=np.mean(ncfile_cur.variables[vname][:,:,:],axis=0)
        outputtime[iday]=date2num( date_curstep,units=units_cur,calendar=calendar_cur)
        print(date_curstep)
      else:
        print("STOP! one wrfout is incomplete %s ",(filename))
        import sys
        sys.exit()
    else:
      print("STOP! one day is missing bettween output %s and %s",(filenames[iday],filenames[iday-1]))
      import sys
      sys.exit()
    ncfile_last=ncfile_cur

  if not ncexist:
    rawnc=createnc(casename,vname,periods,units_cur,calendar_cur,nx,ny,dimension=outputdim)
    rawnc.variables[vname].units=ncfile_last.variables[vname].units
    rawnc.variables[vname].description=ncfile_last.variables[vname].description

  if outputdim==3:
    rawnc.variables[vname][lastindex:,:,:]=outputdata

  rawnc.variables["time"][lastindex:]=outputtime
#clat0[:]=ncfile_last.variables['XLAT']
#clon0[:]=ncfile_last.variables['XLON']


########################DIAG PART#################################
  if periods=="daily":
    if vname=="PRAVG":
      postlist=["PCT","PRAVG","RAINYDAYS","R10","R5D","CDD","R95T","SDII"]
    else:
      postlist=[vname]

    nctime=rawnc.variables["time"]
    start_ymd=num2date(nctime[0],units=units_cur,calendar=calendar_cur)
    end_ymd  =num2date(nctime[-1],units=units_cur,calendar=calendar_cur)
    dodiag=False
    diagnc={}
    if datetime(start_ymd.year,12,1)+shiftday*Oneday>=start_ymd and datetime(end_ymd.year,11,1)<=end_ymd and start_ymd.year<end_ymd.year: 
      diagfname="%s_%s_seasonal.nc"%(casename,vname)
      diagexist=os.path.isfile(diagfname) 
      if diagexist:
        for postvar in postlist:
          diagfname="%s_%s_seasonal.nc"%(casename,postvar)
          diagnc[postvar]=Dataset(diagfname,'a')
        nctime_diag=diagnc[vname].variables["time"]
        lastindex=len(nctime_diag)
        if nctime_diag[-1]<end_ymd.year:
          print("do seasonal mean diagnostic analysis")
          dodiag=True
          diag_startyear=nctime_diag[-1]+1
          diag_endyear=end_ymd.year
      else:
        lastindex=0
        dodiag=True
        diag_startyear=start_ymd.year+1
        diag_endyear=end_ymd.year
        for postvar in postlist:
          diagnc[postvar]=createseasonalnc(casename,postvar,units_cur,nx,ny)

    if dodiag:
      Years=range(diag_startyear,diag_endyear+1)
      for postvar in postlist:
        diagnc[postvar].variables['time'][lastindex:]=Years
      print(Years)
      for i,year in  enumerate(Years):
        i_cur=i+lastindex
        for j,month in enumerate(seasonList):
          if 12 in month:
            byear=year-1
          else:
            byear=year
          eyear=year
          bmonth=month[0]
          emonth=month[2]
          bday=1
          if calendar_cur=="noleap":
            eday=calendar.monthrange(1999, emonth)[1]
          else:
            eday=calendar.monthrange(eyear, emonth)[1]
          ymd_datetime=datetime(int(byear),int(bmonth),1,0,0,0)
          dayb=date2num(ymd_datetime,units=units_cur,calendar=calendar_cur)
          print(ymd_datetime)
          ymd_datetime=datetime(int(eyear),int(emonth),int(eday),0,0,0)
          print(ymd_datetime)
          daye=date2num(ymd_datetime,units=units_cur,calendar=calendar_cur)
          dayb=dayb-rawnc.variables["time"][0]+shiftday
          daye=daye-rawnc.variables["time"][0]+shiftday
          data_daily_ma=ma.masked_values(rawnc.variables[vname][int(dayb):int(daye),:,:],1.e+20)
          if vname=="PRAVG":
            (diagnc["RAINYDAYS"].variables["RAINYDAYS"][i_cur,j,:,:],
             diagnc["R10"].variables["R10"][i_cur,j,:,:],
             diagnc["R5D"].variables["R5D"][i_cur,j,:,:],
             diagnc["SDII"].variables["SDII"][i_cur,j,:,:],
             diagnc["R95T"].variables["R95T"][i_cur,j,:,:])=cs_stat.precp_extrem(fields=data_daily_ma,r95t_hist=default.r95t_hist.variables["R95T_hist"][j],dry_lim=dry_lim)
            diagnc["PCT"].variables["PCT"][i_cur,j,:,:]=cs_stat.quantile_cal(data_daily_ma,dry_lim,pct)
            diagnc["CDD"].variables["CDD"][i_cur,j,:,:]=cs_stat.consective_dry(fields=data_daily_ma,dry_lim=dry_lim)
          diagnc[vname].variables[vname].units=ncfile_last.variables[vname].units
          diagnc[vname].variables[vname].description=ncfile_last.variables[vname].description
          diagnc[vname].variables[vname][i_cur,j,:,:]=np.mean(data_daily_ma,axis=0)
          print("year %s season %s"% (str(year),str(j)))
      for postvar in postlist:
        diagnc[postvar].close()

  rawnc.close() #flush out rawnc
