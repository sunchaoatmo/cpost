import numpy as np

def anal_daily(iday,outputdata,wrf_o,wrf_i,taskname,fields,vert_intp,outputdim,z_levs,number_of_zlevs,compute_mode):
  from constant import RCP,p0
  from ARWpost import arwpost
  import time
  pb   =wrf_i.variables['PB'][0,:,:,:]
  ntime       =wrf_o.dimensions['Time'].size
  south_north =wrf_o.dimensions['south_north'].size
  west_east   =wrf_o.dimensions['west_east'].size
  bottom_top  =wrf_o.dimensions['bottom_top'].size
  if vert_intp=="p":
    hgt  =wrf_i.variables['HGT'][0,:,:]
    if taskname=="geopt" or taskname=="height" or taskname=="temp" :
      phb  =wrf_i.variables['PHB'][0,:,:,:]
    if taskname=="uv_met":
      cosalpha = wrf_i.variables['COSALPHA'][0]
      sinalpha = wrf_i.variables['SINALPHA'][0]
    tk   =None
    geopt=None
    qv   =None
    psfc =None
    nz,ny,nx=bottom_top,south_north, west_east
    for itime in range(ntime):
      p    =wrf_o.variables['P'][itime,:,:,:]
      pres =p+pb
      if taskname=="geopt" or taskname=="height" or taskname=="temp" :
        ph   =wrf_o.variables['PH'][itime,:,:,:]
        psfc =wrf_o.variables['PSFC'][itime,:,:]
        qv   =wrf_o.variables['QVAPOR'][itime,:,:,:]
        geopt_w=ph+phb
        geopt=(geopt_w[1:,:,:]+geopt_w[:-1,:,:])/2.0
        t  =wrf_o.variables['T'][itime,:,:,:]
        tk = (t+300.) * ( pres / p0 )**RCP
        theta = (t+300.)

      for field in fields:
        if field in wrf_o.variables:
          metfield=wrf_o.variables[field][itime,:,:,:]
          (ntime,nz,ny,nx)=wrf_o.variables[field].shape
        elif field=="tk":
          metfield=tk
        elif field=="theta":
          metfield=theta
        elif field=="geopt" :
          metfield=geopt
        elif field=="height":
          metfield=geopt/9.8
        elif field=="u_met":
          metfield=wrf_o.variables["U"][itime,:,:,:]
          (ntime,nz,ny,nx)=wrf_o.variables["U"].shape
        elif field=="v_met":
          metfield=wrf_o.variables["V"][itime,:,:,:]
          (ntime,nz,ny,nx)=wrf_o.variables["V"].shape

        """
        elif taskname=="td":
        elif taskname=="cldfr":
        elif taskname=="dbz":
        """
#       outputdata[field][iday,:,:,:]+= arwpost.interp(
        outputdata[field][iday,:,:,:]+= arwpost.interp(
#             data_out=temp,
             cname=field                     , vertical_type=vert_intp         , 
             data_in=metfield                , 
             z_data=pres                     , 
             number_of_zlevs=number_of_zlevs , z_levs=z_levs                   , 
             psfc=psfc                       , hgt=hgt                         , pres=pres                    , 
             geopt=geopt                     , tk=tk                           , qv=qv                        , 
             nx =nx                          , ny =ny                          , nz=nz                        , 
             bottom_top_dim=bottom_top       , south_north_dim=south_north     , west_east_dim=west_east)
    for field in fields:
      outputdata[field][iday,:,:,:]=outputdata[field][iday,:,:,:]/ntime
    if taskname=="uv_met":
      ur=outputdata["u_met"][iday,:,:,:]
      vr=outputdata["v_met"][iday,:,:,:]
      ue = ur * cosalpha - vr * sinalpha
      ve = vr * cosalpha + ur * sinalpha
      outputdata["u_met"][iday,:,:,:]=ue
      outputdata["v_met"][iday,:,:,:]=ve
  else:
    if taskname=="conv":
      cape=np.zeros((ntime,south_north,west_east))
      cin =np.zeros((ntime,south_north,west_east))
      for itime in range(ntime):
        p    =wrf_o.variables['P'][itime,:,:,:]
        pres =p+pb
        qv   =wrf_o.variables['QVAPOR'][itime,:,:,:]
        t  =wrf_o.variables['T'][itime,:,:,:]
        tc =-273.15+(t+300.) * ( pres / p0 )**RCP
        arwpost.getcape_3d( pres=pres , tc=tc, qv=qv, cape=cape[itime,:,:] , cin=cin[itime,:,:],
                            bottom_top_dim= bottom_top, south_north_dim=south_north , west_east_dim=west_east)
    for field in fields:
      if field=="CAPE":
        metfield=cape
      elif field=="CIN":
        metfield=cin
      else:
        metfield=wrf_o.variables[taskname]
      if outputdim==3:
        if compute_mode==1:
          outputdata[field][iday,:,:]=np.mean(metfield,axis=0)
        elif compute_mode==8:
          outputdata[field][iday,:,:]=np.max(metfield,axis=0)
        elif compute_mode==9:
          outputdata[field][iday,:,:]=np.min(metfield,axis=0)
      elif  outputdim==4:
        if compute_mode==1:
          outputdata[field][iday,:,:,:]=np.mean(metfield,axis=0)
        elif compute_mode==8:
          outputdata[field][iday,:,:,:]=np.max(metfield,axis=0)
        elif compute_mode==9:
          outputdata[field][iday,:,:,:]=np.min(metfield,axis=0)
      else:
        print("can only output to 3 or 4 dim")

  return 

def anal_sea_mon(periods,rawnc,monthList,fields,taskname,casename,shiftday,calendar_cur,time_units,units,description,nx,ny,nz,r95t_hist):
  import numpy.ma as ma
  from cs_stat import cs_stat
  from netCDF4 import date2num,num2date
  import calendar
  import os.path
  from netCDF4 import Dataset
  from datetime import datetime,timedelta
  from constant import seasonList,monthlyList,dry_lim,qvalue # only several constants
  import time
  from writenc import createsmnc
  dodiag=False
  nctime=rawnc.variables["time"]
  start_ymd=num2date(nctime[0],units=time_units,calendar=calendar_cur)
  end_ymd  =num2date(nctime[-1],units=time_units,calendar=calendar_cur)
  beg_num =date2num( datetime(start_ymd.year,12,1),units=time_units,calendar=calendar_cur)
  end_num =date2num( datetime(  end_ymd.year,12,1),units=time_units,calendar=calendar_cur)

  if beg_num+shiftday>=nctime[0] and end_num<=nctime[-1] and start_ymd.year<end_ymd.year: 
    diagfname="%s_%s_%s.nc"%(casename,taskname,periods)
    diagexist=os.path.isfile(diagfname) 
    if diagexist:
      diagfname="%s_%s_%s.nc"%(casename,taskname,periods)
      diagnc=Dataset(diagfname,'a')
      nctime_diag=diagnc.variables["time"]
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
      diagnc=createsmnc(casename,taskname,periods,fields,nx,ny,nz)
      for field in fields:
        diagnc.variables[field].units=units[field]
        diagnc.variables[field].description=description[field]

  if dodiag:
    Years=range(diag_startyear,diag_endyear+1)
    diagnc.variables['time'][lastindex:]=Years
    for i,year in  enumerate(Years):
      i_cur=i+lastindex
      for j,month in enumerate(monthList):
        if 12 == month[0]:
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
        print("%s analysis beg:%s"%(periods,ymd_datetime))
        ymd_datetime=datetime(int(eyear),int(emonth),int(eday),0,0,0)
        print("%s analysis end:%s"%(periods,ymd_datetime))
        daye=date2num(ymd_datetime,units=time_units,calendar=calendar_cur)
        dayb=dayb-rawnc.variables["time"][0]+shiftday
        daye=daye-rawnc.variables["time"][0]+shiftday
        if taskname=="PR":
          data_daily_ma=ma.masked_values(rawnc.variables["PRAVG"][int(dayb):int(daye),:,:],1.e+20)
          if periods=="seasonal":
            R95T_HIST=r95t_hist.variables["R95T_hist"][j]
          else:
            R95T_HIST=None
          (diagnc.variables["RAINYDAYS"][i_cur,j,:,:],
           diagnc.variables["R10"][i_cur,j,:,:],
           diagnc.variables["R5D"][i_cur,j,:,:],
           diagnc.variables["SDII"][i_cur,j,:,:],
           diagnc.variables["R95T"][i_cur,j,:,:])=cs_stat.precp_extrem(fields=data_daily_ma,r95t_hist=R95T_HIST,dry_lim=dry_lim)
          diagnc.variables["PCT"][i_cur,j,:,:]=cs_stat.quantile_cal(pre_quantile=data_daily_ma,dry_lim=dry_lim,qvalue=qvalue)
          diagnc.variables["CDD"][i_cur,j,:,:]=cs_stat.consective_dry(fields=data_daily_ma,dry_lim=dry_lim)
          diagnc.variables["PRAVG"][i_cur,j,:,:]=np.mean(data_daily_ma,axis=0)
        else:
          for field in fields:
            data_daily_ma=ma.masked_values(rawnc.variables[field][int(dayb):int(daye),:,:],1.e+20)
            diagnc.variables[field][i_cur,j,:,:]=np.mean(data_daily_ma,axis=0)
        print("year %s season %s"% (str(year),str(j)))
    diagnc.history ='Created by Chao Sun sunchao@umd.edu ' + time.ctime(time.time())
    diagnc.source ='CWRF run:%s'%casename
    diagnc.close()

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
