def diag_season_month(periods,rawnc,monthList,start_ymd,end_ymd,postlist,vname,casename,shiftday,calendar_cur,time_units,units,description,nx,ny):
  from default import r95t_hist
  import numpy.ma as ma
  import numpy as np
  from cs_stat import cs_stat
  from netCDF4 import date2num,num2date
  import calendar
  import os.path
  from netCDF4 import Dataset
  from datetime import datetime,timedelta
  from constant import *   # only several constants
  from writenc import createsmnc
  dodiag=False
  diagnc={}
  Oneday=timedelta(days=1)

  if datetime(start_ymd.year,12,1)+shiftday*Oneday>=start_ymd and datetime(end_ymd.year,11,1)<=end_ymd and start_ymd.year<end_ymd.year: 
    diagfname="%s_%s_%s.nc"%(casename,vname,periods)
    diagexist=os.path.isfile(diagfname) 
    if diagexist:
      for postvar in postlist:
        diagfname="%s_%s_%s.nc"%(casename,postvar,periods)
        diagnc[postvar]=Dataset(diagfname,'a')
      nctime_diag=diagnc[vname].variables["time"]
      lastindex=len(nctime_diag)
      if nctime_diag[-1]<end_ymd.year:
        print("do %s mean diagnostic analysis"%periods)
        dodiag=True
        diag_startyear=nctime_diag[-1]+1
        diag_endyear=end_ymd.year
    else:
      lastindex=0
      dodiag=True
      diag_startyear=start_ymd.year+1
      diag_endyear=end_ymd.year
      for postvar in postlist:
        diagnc[postvar]=createsmnc(casename,postvar,periods,nx,ny)

  if dodiag:
    Years=range(diag_startyear,diag_endyear+1)
    for postvar in postlist:
      diagnc[postvar].variables['time'][lastindex:]=Years
    for i,year in  enumerate(Years):
      i_cur=i+lastindex
      for j,month in enumerate(monthList):
        if 12 in month:
          byear=year-1
        else:
          byear=year
        eyear=year
        bmonth=month[0]
        emonth=month[1]
        bday=1
        if calendar_cur=="noleap":
          eday=calendar.monthrange(1999, emonth)[1]
        else:
          eday=calendar.monthrange(eyear, emonth)[1]
        ymd_datetime=datetime(int(byear),int(bmonth),1,0,0,0)
        dayb=date2num(ymd_datetime,units=time_units,calendar=calendar_cur)
        print(ymd_datetime)
        ymd_datetime=datetime(int(eyear),int(emonth),int(eday),0,0,0)
        print(ymd_datetime)
        daye=date2num(ymd_datetime,units=time_units,calendar=calendar_cur)
        dayb=dayb-rawnc.variables["time"][0]+shiftday
        daye=daye-rawnc.variables["time"][0]+shiftday
        data_daily_ma=ma.masked_values(rawnc.variables[vname][int(dayb):int(daye),:,:],1.e+20)
        if vname=="PRAVG":
          (diagnc["RAINYDAYS"].variables["RAINYDAYS"][i_cur,j,:,:],
           diagnc["R10"].variables["R10"][i_cur,j,:,:],
           diagnc["R5D"].variables["R5D"][i_cur,j,:,:],
           diagnc["SDII"].variables["SDII"][i_cur,j,:,:],
           diagnc["R95T"].variables["R95T"][i_cur,j,:,:])=cs_stat.precp_extrem(fields=data_daily_ma,r95t_hist=r95t_hist.variables["R95T_hist"][j],dry_lim=dry_lim)
          diagnc["PCT"].variables["PCT"][i_cur,j,:,:]=cs_stat.quantile_cal(data_daily_ma,dry_lim,pct)
          diagnc["CDD"].variables["CDD"][i_cur,j,:,:]=cs_stat.consective_dry(fields=data_daily_ma,dry_lim=dry_lim)
        diagnc[vname].variables[vname].units=units
        diagnc[vname].variables[vname].description=description
        diagnc[vname].variables[vname][i_cur,j,:,:]=np.mean(data_daily_ma,axis=0)
        print("year %s season %s"% (str(year),str(j)))
    for postvar in postlist:
      diagnc[postvar].close()

def wrftimetodate(wrfstr):
  from datetime import datetime
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
