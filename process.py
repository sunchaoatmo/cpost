import numpy as np

def anal_daily(iday,outputdata,wrf_o,wrf_i,taskname,fields,vert_intp,
               outputdim,z_levs,number_of_zlevs,compute_mode,
               wrfncfile_last,wrfncfile_next):
  from constant import RCP,p0,G,Rd,EPS,missing,fill_nocloud,opt_thresh,cwrfnames
  from ARWpost import arwpost
  import time
  pb   =wrf_i.variables['PB'][0,:,:,:]
  phb  =wrf_i.variables['PHB'][0,:,:,:]
  hgt  =wrf_i.variables['HGT'][0,:,:]
  ntime       =wrf_o.dimensions['Time'].size
  ntime_recip =1.0/float(ntime)
  south_north =wrf_o.dimensions['south_north'].size
  west_east   =wrf_o.dimensions['west_east'].size
  bottom_top  =wrf_o.dimensions['bottom_top'].size
  if taskname in ["uv_met","uv_10","Qflux"]:
    cosalpha = wrf_i.variables['COSALPHA'][0]
    sinalpha = wrf_i.variables['SINALPHA'][0]
  if vert_intp=="p":
    tk   =None
    geopt=None
    qv   =None
    psfc =None
    nz,ny,nx=bottom_top,south_north, west_east

    for itime in range(ntime):
      p    =wrf_o.variables['P'][itime,:,:,:]
      pres =p+pb
      if taskname in ["geopt", "height" ,"temp" ,"omega","Qflux"]:
        ph   =wrf_o.variables['PH'][itime,:,:,:]
        psfc =wrf_o.variables['PSFC'][itime,:,:]
        qv   =wrf_o.variables['QVAPOR'][itime,:,:,:]
        geopt_w=ph+phb
        geopt=(geopt_w[1:,:,:]+geopt_w[:-1,:,:])/2.0
        t  =wrf_o.variables['T'][itime,:,:,:]
        tk = (t+300.) * ( pres / p0 )**RCP
        theta = (t+300.)

      if taskname in ["Qflux"]:
        U=wrf_o.variables["U"][itime,:,:,:]
        V=wrf_o.variables["V"][itime,:,:,:]
        U= np.asarray(U)
        V= np.asarray(V)
        V= (V[:,1:,:]+V[:,:-1,:])*0.5
        U= (U[:,:,1:]+U[:,:,:-1])*0.5
        UE = U * cosalpha - V * sinalpha
        VE = V * cosalpha + U * sinalpha
        q  = wrf_o.variables["QVAPOR"][itime,:,:,:]
        q  = np.asarray(q)

      if taskname=="rh3d":
        if itime==0:
          rh3d=np.zeros((bottom_top,south_north,west_east)) #,order='F',dtype=np.float32)
        ph   =wrf_o.variables['PH'][itime,:,:,:]
        psfc =wrf_o.variables['PSFC'][itime,:,:]
        qv   =wrf_o.variables['QVAPOR'][itime,:,:,:]
        geopt_w=ph+phb
        geopt=(geopt_w[1:,:,:]+geopt_w[:-1,:,:])/2.0
        t  =wrf_o.variables['T'][itime,:,:,:]
        tk = (t+300.) * ( pres / p0 )**RCP
        for ilevel in range(bottom_top):
          rh3d[ilevel,:,:]    =arwpost.calc_rh( q2m=qv[ilevel,:,:],t2m=tk[ilevel,:,:],
                                                psfc=pres[ilevel,:,:])

      if taskname=="omega":
        w  =wrf_o.variables['W'][itime,:,:,:]
        www= .5*(w[1:,:,:] + w[:-1,:,:])
        omega = -G*pres/ (Rd*((tk*(EPS + qv))/ (EPS*(1. + qv))))*www

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
        elif field=="omega":
          metfield=omega
        elif field=="rh3d":
          metfield=rh3d
        elif field=="uq":
          metfield=UE*q
        elif field=="vq":
          metfield=VE*q

        outputdata[field][:,:,:]+= arwpost.interp(
             cname=field                     , vertical_type=vert_intp         , 
             data_in=metfield                , 
             z_data=pres                     , 
             number_of_zlevs=number_of_zlevs , z_levs=z_levs                   , 
             psfc=psfc                       , hgt=hgt                         , pres=pres                    , 
             geopt=geopt                     , tk=tk                           , qv=qv                        , 
             nx =nx                          , ny =ny                          , nz=nz                        , 
             bottom_top_dim=bottom_top       , south_north_dim=south_north     , west_east_dim=west_east)
    for field in fields:
      outputdata[field][:,:,:]=outputdata[field][:,:,:]/ntime

    if taskname =="uv_met":
      ur=outputdata["u_met"][:,:,:]
      vr=outputdata["v_met"][:,:,:]
      ue = ur * cosalpha - vr * sinalpha
      ve = vr * cosalpha + ur * sinalpha
      #win= np.sqrt(ur*ur+vr*vr)
      #outputdata["WIN"][:,:,:]=win
      outputdata["u_met"][:,:,:]=ue
      outputdata["v_met"][:,:,:]=ve


  else:
    if taskname=="conv":
      cape=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      cin =np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      for itime in range(ntime):
        p    =wrf_o.variables['P'][itime,:,:,:]
        pres =p+pb
        qv   =wrf_o.variables['QVAPOR'][itime,:,:,:]
        t  =wrf_o.variables['T'][itime,:,:,:]
        tc =-273.15+(t+300.) * ( pres / p0 )**RCP
        tk = (t+300.) * ( pres / p0 )**RCP
        ph   =wrf_o.variables['PH'][itime,:,:,:]
        geopt_w=ph+phb
        psfc =wrf_o.variables['PSFC'][itime,:,:]
        geopt=(geopt_w[1:,:,:]+geopt_w[:-1,:,:])*0.5
        arwpost.calc_cape(cape_out=cape, cin_out=cin,itime=itime+1, 
                      hgt=hgt,  qv_in=qv,  pres_in=pres, tk_in=tk, geopt_in=geopt,psfc=psfc,
             bottom_top_dim=bottom_top       , south_north_dim=south_north     , west_east_dim=west_east,ntime=ntime)
    elif taskname=="RH":
      rh=np.zeros((ntime,south_north,west_east)) #,order='F',dtype=np.float32)
      for itime in range(ntime):
        q2m   =wrf_o.variables['Q2M'][itime,:,:]
        t2m   =wrf_o.variables['AT2M'][itime,:,:]
        #t2m   =wrf_o.variables['T2M'][itime,:,:]
        psfc  =wrf_o.variables['PSFC'][itime,:,:]
        rh[itime,:,:]    =arwpost.calc_rh( q2m=q2m,t2m=t2m,psfc=psfc)
    elif taskname=="cldfrag":
      cldfra_low =np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      cldfra_mid =np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      cldfra_high=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      cldfra_total=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      for itime in range(ntime):
        cldfra=wrf_o.variables['CLDFRA'][itime,:,:]
        p     =wrf_o.variables['P'][itime,:,:,:]
        pres  =p+pb
        presmb=pres/100
        arwpost.calc_cldfra( pres=presmb,cldfra=cldfra,
                             cldfra_total=cldfra_total,
                             cldfra_low=cldfra_low,
                             cldfra_mid=cldfra_mid,
                             cldfra_high=cldfra_high,itime=itime)
    elif taskname=="slp":
      slp=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      for itime in range(ntime):
        p    =wrf_o.variables['P'][itime,:,:,:]
        pres =p+pb
        qv   =wrf_o.variables['QVAPOR'][itime,:,:,:]
        t  =wrf_o.variables['T'][itime,:,:,:]
        tk = (t+300.) * ( pres / p0 )**RCP
        ph   =wrf_o.variables['PH'][itime,:,:,:]
        geopt_w=ph+phb
        geopt=(geopt_w[1:,:,:]+geopt_w[:-1,:,:])*0.5
        arwpost.calc_slp(pres=pres, geopt=geopt,tk=tk, qv=qv,slp=slp,
                         itime=itime,
                         bottom_top_dim =bottom_top , 
                         south_north_dim=south_north,
                         west_east_dim=west_east,ntime=ntime)
    #cloud ralated
    elif taskname in ["TCL","TCR","TCI","TCS","TCG","TPW","TPWmax"]:
      tc_l=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      tc_m=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      tc_h=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      for itime in range(ntime):
        p    =wrf_o.variables['P'][itime,:,:,:]
        pres =p+pb
        q    =wrf_o.variables[cwrfnames[taskname]][itime,:,:,:]
        psfc =wrf_o.variables['PSFC'][itime,:,:]
        arwpost.calc_tpw(pres=pres,qv=q,psfc=psfc,
                         tpw_h=tc_h,
                         tpw_m=tc_m,
                         tpw_l=tc_l,
                         itime=itime,
                         bottom_top_dim =bottom_top , 
                         south_north_dim=south_north,
                         west_east_dim=west_east,ntime=ntime)
    #cloud ralated

    elif taskname=="lwp":
      lwp=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      for itime in range(ntime):
        p    =wrf_o.variables['P'][itime,:,:,:]
        pres =p+pb
        qc   =wrf_o.variables['QCLOUD'][itime,:,:,:]
        qr   =wrf_o.variables['QRAIN'][itime,:,:,:]
        psfc =wrf_o.variables['PSFC'][itime,:,:]
        arwpost.calc_lwp(pres=pres,qc=qc,qr=qr,psfc=psfc,
                         lwp=lwp,
                         itime=itime,
                         bottom_top_dim =bottom_top , 
                         south_north_dim=south_north,
                         west_east_dim=west_east,ntime=ntime)
    elif taskname=="iwp":
      iwp=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      for itime in range(ntime):
        p    =wrf_o.variables['P'][itime,:,:,:]
        pres =p+pb
        qs   =wrf_o.variables['QSNOW'][itime,:,:,:]
        qg   =wrf_o.variables['QGRAUP'][itime,:,:,:]
        qi   =wrf_o.variables['QICE'][itime,:,:,:]
        psfc =wrf_o.variables['PSFC'][itime,:,:]
        arwpost.calc_iwp(pres=pres,qs=qs,qg=qg,qi=qi,psfc=psfc,
                         iwp=iwp,
                         itime=itime,
                         bottom_top_dim =bottom_top , 
                         south_north_dim=south_north,
                         west_east_dim=west_east,ntime=ntime)
    elif taskname=="cto":
      ctt=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      cth=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      for itime in range(ntime):
        p    =wrf_o.variables['P'][itime,:,:,:]
        pres =p+pb
        t  =wrf_o.variables['T'][itime,:,:,:]
        tk = (t+300.) * ( pres / p0 )**RCP
        qc   =wrf_o.variables['QCLOUD'][itime,:,:,:]
        qcw  =qc *1000
        qv   =wrf_o.variables['QVAPOR'][itime,:,:,:]
        qvp  =qv *1000
        pres_hpa =pres*0.01 # convert to hpa for the ctt
        psfc =wrf_o.variables['PSFC'][itime,:,:]
        psfc_hpa =psfc*0.01 # convert to hpa for the ctt
        try:
          qi   =wrf_o.variables['QICE'][itime,:,:,:]
          qice =qi *1000
          haveqci=1
        except:
          haveqci=0
          qice =np.zeros(qv.shape, qv.dtype)
        ph   =wrf_o.variables['PH'][itime,:,:,:]
        geopt_w=ph+phb
        geopt=(geopt_w[1:,:,:]+geopt_w[:-1,:,:])*0.5
        ght  =geopt/9.8

        arwpost.wrfcttcalc(prs=pres,tk=tk,qci=qice,qcw=qcw,qvp=qvp,ght=ght,ter=hgt,
                          psfc=psfc,
                         ctt=ctt,
                         cth=cth,
                         haveqci=haveqci,
                         fill_nocloud=fill_nocloud,
                         missing=missing,
                         opt_thresh=opt_thresh,
                         itime=itime,
                         nz =bottom_top , 
                         ns =south_north,
                         ew =west_east,ntime=ntime)
    elif taskname in ["xsmtg","xsmlg","xsmig"]:
      nsmg=4
      if not hasattr(anal_daily,"xzi"):
        try:
          nl_soil_stag=wrf_o.dimensions['nl_soil_stag'].size
        except:
          nl_soil_stag=11
        nl_soilpsnow=wrf_o.dimensions['nl_soilpsnow_stag'].size
        lb_soil= nl_soilpsnow-nl_soil_stag
        xz=wrf_i.variables['XZ']
        xzi=np.zeros((nl_soil_stag+1,south_north,west_east))
        xzi[1:-1,:,:]=(xz[0,lb_soil+1:,:,:]+xz[0,lb_soil:-1,:,:])*0.5
        xzi[-1,:,:]=5.67643 #xz[0,-1,:,:]
        anal_daily.xzi=xzi*1000 #change unit from m to mm
        anal_daily.xsmtg=np.zeros((ntime,nsmg,south_north,west_east),order='F',dtype=np.float32)
      xsmtg=anal_daily.xsmtg

      for itime in range(ntime):
        xwliq=wrf_o.variables['XWLIQ'][itime,:,:]
        xwice=wrf_o.variables['XWICE'][itime,:,:]
        if taskname=="xsmtg":
          xwater=xwliq+xwice
        elif taskname=="xsmlg":
          xwater=xwliq
        elif taskname=="xsmig":
          xwater=xwice
        arwpost.calc_sm( xwater=xwater,xzsoil=anal_daily.xzi,xsmtg=xsmtg,itime=itime)
    elif taskname=="rh3d":
      rh3d=np.zeros((ntime,bottom_top,south_north,west_east)) #,order='F',dtype=np.float32)

      for itime in range(ntime):
        ph   =wrf_o.variables['PH'][itime,:,:,:]
        psfc =wrf_o.variables['PSFC'][itime,:,:]
        qv   =wrf_o.variables['QVAPOR'][itime,:,:,:]
        geopt_w=ph+phb
        geopt=(geopt_w[1:,:,:]+geopt_w[:-1,:,:])/2.0
        t  =wrf_o.variables['T'][itime,:,:,:]
        p    =wrf_o.variables['P'][itime,:,:,:]
        pres =p+pb
        tk = (t+300.) * ( pres / p0 )**RCP
        for ilevel in range(bottom_top):
          rh3d[itime,ilevel,:,:]    =arwpost.calc_rh( q2m=qv[ilevel,:,:],t2m=tk[ilevel,:,:],
                                                psfc=pres[ilevel,:,:])




    for field in fields:
      if field=="CAPE":
        metfield=cape
      elif field=="CIN":
        metfield=cin
      elif field=="RH":
        metfield=rh
      elif field=="u_10":
        metfield=wrf_o.variables["AU10"][:,:,:]
      elif field=="v_10":
        metfield=wrf_o.variables["AV10"][:,:,:]
      elif field=="cldfra_low":
        metfield=cldfra_low
      elif field=="cldfra_mid":
        metfield=cldfra_mid
      elif field=="cldfra_high":
        metfield=cldfra_high
      elif field=="cldfra_total":
        metfield=cldfra_total
      elif field=="slp":
        metfield=slp
      # cloud ralated 
      elif field in ["TPW_l","TCL_l","TCR_l","TCI_l","TCS_l","TCG_l"]:
        metfield=tc_l*0.1 #convert to cm
      elif field in ["TPW_m","TCL_m","TCR_m","TCI_m","TCS_m","TCG_m"]:
        metfield=tc_m*0.1 # convert to cm
      elif field in ["TPW_h","TCL_h","TCR_h","TCI_h","TCS_h","TCG_h"]:
        metfield=tc_h*0.1 # convert to cm
      # cloud ralated 
      elif field=="lwp":
        metfield=lwp
      elif field=="iwp":
        metfield=iwp
      elif field=="ctt":
        metfield=ctt
      elif field=="cth":
        metfield=cth
      elif field=="rh3d":
        metfield=rh3d
      elif field in ["xsmtg","xsmlg","xsmig"]:
        metfield=xsmtg
      #elif field=="WIN_10":
      #  continue 
      else:
        metfield=wrf_o.variables[field]
      if outputdim==3:
        if compute_mode==1:
          if field in ["ctt","cth"]:
            outputdata[field][:,:]=arwpost.aveexceptmissing(
                         met_3d=metfield,
                         missing=missing,
                         nz =ntime , 
                         ns =south_north,
                         ew =west_east)
          elif wrfncfile_next:
            outputdata[field][:,:]=np.sum(metfield[1:,:,:],axis=0)
            outputdata[field][:,:]+=wrfncfile_next.variables[field][0,:,:]
            outputdata[field][:,:]=outputdata[field][:,:]*ntime_recip
          else:
            outputdata[field][:,:]=np.mean(metfield[:,:,:],axis=0)
        elif compute_mode==2:
          outputdata[field][:,:]=np.mean(np.sum(metfield,axis=1),axis=0)
        elif compute_mode==8:
          outputdata[field][:,:]=np.max(metfield,axis=0)
        elif compute_mode==9:
          outputdata[field][:,:]=np.min(metfield,axis=0)
      elif  outputdim==4:
        if compute_mode==1:
          outputdata[field][:,:,:]=np.mean(metfield,axis=0)
        elif compute_mode==8:
          outputdata[field][:,:,:]=np.max(metfield,axis=0)
        elif compute_mode==9:
          outputdata[field][:,:,:]=np.min(metfield,axis=0)
      else:
        print("can only output to 3 or 4 dim")
    if taskname=="uv_10":
      ur=outputdata["u_10"][:,:]
      vr=outputdata["v_10"][:,:]
      ue = ur * cosalpha - vr * sinalpha
      ve = vr * cosalpha + ur * sinalpha
      #win= np.sqrt(ur*ur+vr*vr)
      #outputdata["WIN_10"][:,:]=win
      outputdata["u_10"][:,:]=ue
      outputdata["v_10"][:,:]=ve
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
  from constant import seasonnames
  import time
  from writenc import createsmnc
  nctime=rawnc.variables["time"]
  start_ymd=num2date(nctime[0],units=time_units,calendar=calendar_cur)
  end_ymd  =num2date(nctime[-1],units=time_units,calendar=calendar_cur)
  beg_num =date2num( datetime(start_ymd.year,12,1),units=time_units,calendar=calendar_cur)
  end_num =date2num( datetime(  end_ymd.year,11,30),units=time_units,calendar=calendar_cur)

  #if beg_num+shiftday>=nctime[0] and end_num<=nctime[-1] and start_ymd.year<end_ymd.year: 
  if periods=="seasonal":
    if start_ymd.year>=end_ymd.year : 
       print("less than one year no seasonal output")
       return None
  if periods=="monthly":
    if start_ymd.month==end_ymd.month and start_ymd.year==end_ymd.year:
       print("less than one month no monthly output")
       return None

  diagfname="%s_%s_%s.nc"%(casename,taskname,periods)
  diagexist=os.path.isfile(diagfname) 
  lastindex=0
  diag_startyear=start_ymd.year#+1
  diag_endyear=end_ymd.year

  try: 
    diagnc=Dataset(diagfname,'a')
    nctime_diag=diagnc.variables["time"]
    lastindex=len(nctime_diag)
    assert  lastindex>0
    if nctime_diag[-1]<end_ymd.year:
      print("Following previous output")
      diag_startyear=nctime_diag[-1]+1
      diag_endyear=end_ymd.year
  except:
    print("create a new")
    diagnc=createsmnc(casename,taskname,periods,fields,nx,ny,nz)
    for field in fields:
      diagnc.variables[field].units=units[field]
      diagnc.variables[field].description=description[field]

  #diag_startyear=1989
  #diag_endyear=2009

  Years=range(diag_startyear,diag_endyear+1)
  diagnc.variables['time'][lastindex:]=Years
  for i,year in  enumerate(Years):
    i_cur=lastindex+i
    for j,month in enumerate(monthList):
      byear=year
      if periods=="seasonal" and 12 == month[0]:
        byear=year-1
      eyear=year
      bmonth=month[0]
      ymd_datetime_b=datetime(int(byear),int(bmonth),1,0,0,0)

      emonth=month[1]
      eday=calendar.monthrange(eyear, emonth)[1]
      if emonth==2:
        if calendar_cur in ["365_day","noleap"]:
           eday=28
        elif calendar_cur in ['all_leap', '366_day']:
           eday=29
      ymd_datetime_e=datetime(int(eyear),int(emonth),int(eday),0,0,0)

      index_b =date2num( ymd_datetime_b,units=time_units,calendar=calendar_cur)
      index_e =date2num( ymd_datetime_e,units=time_units,calendar=calendar_cur)
      index_b = index_b- nctime[0]
      index_e = index_e- nctime[0]
      if "hours" in time_units:
        index_b = index_b/24
        index_e = index_e/24


      if index_e<0 or index_e>len(rawnc.variables["time"]):
        continue
      if index_b<=-35:
        continue
      elif index_b<0:
        print("WRANING your first day is not the begining of the month!!!")
        #ymd_datetime_b=datetime(int(byear),int(bmonth),int(1-dayb),0,0,0) # reset the first day
        diagnc.history = "WARNING the first day is not the begining of the month but "
        index_b=0
      if taskname=="PR":
        data_daily_ma=ma.masked_values(rawnc.variables["PRAVG"][int(index_b):int(index_e),:,:],1.e+20)
        data_daily_ma[data_daily_ma>10000]=0.0
        R95T_HIST=None
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
        diagnc.variables["PMAX"][i_cur,j,:,:]=np.max(data_daily_ma,axis=0)
      else:
        for field in fields:
          temp=np.asarray(rawnc.variables[field][int(index_b):int(index_e),:,:])
          data_daily_ma=ma.masked_values(temp,1.e+20)
          diagnc.variables[field][i_cur,j,:,:]=np.mean(data_daily_ma,axis=0)
      if periods=="monthly":
        print("Year %s %s:"% (str(year),calendar.month_abbr[j+1]))
      else: 
        print("Year %s %s:"% (str(year),seasonnames[j]))
      print("     %s analysis beg:%s"%(periods,ymd_datetime_b))
      print("     %s analysis end:%s"%(periods,ymd_datetime_e))
  return diagnc

def wrftimetodate(wrfstr):
  from datetime import datetime
  YEAR=slice(0,4)
  MONTH=slice(5,7)  
  DAY=slice(8,10)  
  HOUR=slice(11,13)
  MINUTE=slice(14,16)
  SECOND=slice(17,19)
  wrfstr=[x.decode("utf-8")  for x in wrfstr]
  Year=int(''.join(list(wrfstr[YEAR])))
  Month=int(''.join(list(wrfstr[MONTH])))
  Day=int(''.join(list(wrfstr[DAY])))
  Hour=int(''.join(list(wrfstr[HOUR])))
  Minute=int(''.join(list(wrfstr[MINUTE])))
  Second=int(''.join(list(wrfstr[SECOND])))
  wrfdate=datetime(Year, Month, Day, Hour, Minute, Second)
  return wrfdate

def anal_hourly(iday,outputdata,wrf_o,wrf_i,taskname,fields,vert_intp,outputdim,z_levs,number_of_zlevs,compute_mode):
  from constant import RCP,p0
  from ARWpost import arwpost
  import time
  pb   =wrf_i.variables['PB'][0,:,:,:]
  ntime       =wrf_o.dimensions['Time'].size
  south_north =wrf_o.dimensions['south_north'].size
  west_east   =wrf_o.dimensions['west_east'].size
  bottom_top  =wrf_o.dimensions['bottom_top'].size
  if "uv" in taskname:
    cosalpha = wrf_i.variables['COSALPHA'][0]
    sinalpha = wrf_i.variables['SINALPHA'][0]
    cosalpha = np.asarray(cosalpha )
    sinalpha = np.asarray(sinalpha )
  if vert_intp=="p":
    hgt  =wrf_i.variables['HGT'][0,:,:]
    if taskname=="geopt" or taskname=="height" or taskname=="temp" :
      phb  =wrf_i.variables['PHB'][0,:,:,:]
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

        outputdata[field][itime,:,:,:]= arwpost.interp(
             cname=field                     , vertical_type=vert_intp         , 
             data_in=metfield                , 
             z_data=pres                     , 
             number_of_zlevs=number_of_zlevs , z_levs=z_levs                   , 
             psfc=psfc                       , hgt=hgt                         , pres=pres                    , 
             geopt=geopt                     , tk=tk                           , qv=qv                        , 
             nx =nx                          , ny =ny                          , nz=nz                        , 
             bottom_top_dim=bottom_top       , south_north_dim=south_north     , west_east_dim=west_east)
    if taskname=="uv_met":
      ur=outputdata["u_met"][:]
      vr=outputdata["v_met"][:]
      ue = ur * cosalpha - vr * sinalpha
      ve = vr * cosalpha + ur * sinalpha
      #win= np.sqrt(ur*ur+vr*vr)
      #outputdata["WIN"][:]=win
      outputdata["u_met"][:]=ue
      outputdata["v_met"][:]=ve
  else:
    if taskname=="conv":
      cape=np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      cin =np.zeros((ntime,south_north,west_east),order='F',dtype=np.float32)
      hgt  =wrf_i.variables['HGT'][0,:,:]
      phb  =wrf_i.variables['PHB'][0,:,:,:]
      for itime in range(ntime):
        p    =wrf_o.variables['P'][itime,:,:,:]
        pres =p+pb
        qv   =wrf_o.variables['QVAPOR'][itime,:,:,:]
        t  =wrf_o.variables['T'][itime,:,:,:]
        tc =-273.15+(t+300.) * ( pres / p0 )**RCP
        tk = (t+300.) * ( pres / p0 )**RCP
        ph   =wrf_o.variables['PH'][itime,:,:,:]
        geopt_w=ph+phb
        psfc =wrf_o.variables['PSFC'][itime,:,:]
        geopt=(geopt_w[1:,:,:]+geopt_w[:-1,:,:])/2.0
        arwpost.calc_cape(cape_out=cape, cin_out=cin,itime=itime, 
                      hgt=hgt,  qv_in=qv,  pres_in=pres, tk_in=tk, geopt_in=geopt,psfc=psfc,
             bottom_top_dim=bottom_top       , south_north_dim=south_north     , west_east_dim=west_east,ntime=ntime)
    elif taskname=="RH":
      rh=np.zeros((ntime,south_north,west_east)) #,order='F',dtype=np.float32)
      for itime in range(ntime):
        q2m   =wrf_o.variables['Q2M'][itime,:,:]
        t2m   =wrf_o.variables['T2M'][itime,:,:]
        psfc  =wrf_o.variables['PSFC'][itime,:,:]
        rh[itime,:,:]    =arwpost.calc_rh( q2m=q2m,t2m=t2m,psfc=psfc)
    for field in fields:
      if field=="CAPE":
        metfield=cape
      elif field=="CIN":
        metfield=cin
      elif field=="RH":
        metfield=rh
      elif field=="u_10":
        metfield=wrf_o.variables["AU10"][:,:,:]
      elif field=="v_10":
        metfield=wrf_o.variables["AV10"][:,:,:]
      #elif field=="WIN_10":
      #  continue 
      else:
        metfield=wrf_o.variables[field]
      if outputdim==3:
        outputdata[field][:,:,:]=metfield[:,:,:]
      elif  outputdim==4:
        outputdata[field][:,:,:,:]=metfield[:,:,:,:]
      else:
        print("can only output to 3 or 4 dim")
    if taskname=="uv_10":
      ur=outputdata["u_10"][:,:,:]
      vr=outputdata["v_10"][:,:,:]
      ue = ur * cosalpha - vr * sinalpha
      ve = vr * cosalpha + ur * sinalpha
      #win= np.sqrt(ur*ur+vr*vr)
      #outputdata["WIN_10"][:,:,:]=win
      outputdata["u_10"][:,:,:]=ue
      outputdata["v_10"][:,:,:]=ve

  return 


