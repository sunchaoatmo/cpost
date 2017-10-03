def chekitegrty(curtime,nstep,filename,filenames,outputtime,timenum):
  import sys
  if not len(curtime)==nstep and filename != filenames[-1]:
    sys.exit("STOP! one wrfout is incomplete %s "%(filename))
  if not outputtime[0]==timenum:
    sys.exit("STOP! %s is missing in wrfout serial "%(timenum))

def openwrfdata(filename):
  from netCDF4 import Dataset
  try:
    wrfncfile_cur=Dataset(filename,'r')
  except:
    import sys
    sys.exit("can not open:%s "%filename)
  return wrfncfile_cur

def getrawnc(rawfname,casename,taskname,periods,units_cur,calendar_cur,var_parameters,nx,ny,nlev,var_units,var_description ):
  from netCDF4 import Dataset
  import os.path
  from writenc import createnc
  try:
    rawnc=Dataset(rawfname,'a')
    lasttime=rawnc.variables["time"][-1]
    print("Following previous postprocessed results %s"%lasttime)
    raw_exit= True 
  except:
    lasttime=None
  raw_exit= True if lasttime else False
  if not raw_exit:
    if os.path.isfile(rawfname):
      os.remove(rawfname)
    rawnc=createnc(casename,taskname,periods,units_cur,calendar_cur,var_parameters[taskname]["fields"].keys(),nx,ny,nlev  )
    for field in var_parameters[taskname]["fields"]:
      rawnc.variables[field].units=var_units[field]
      rawnc.variables[field].description=var_description[field]
  return rawnc,raw_exit

def updatelastwrfoutfn(rawnc,units_cur,calendar_cur,filenames):
  from netCDF4 import Dataset
  from constant import wrfout_data_fmt   # only several constants
  import sys
  from netCDF4 import date2num,num2date
  lastday=num2date(rawnc.variables["time"][-1],units=units_cur,calendar=calendar_cur)
  lastwrfout="wrfout_d01_%s"%lastday.strftime(wrfout_data_fmt)
  try:
    lastindex=filenames.index(lastwrfout)
    del filenames[:lastindex]
    ncfile_last=Dataset(filenames[0],'r') # replace the ncfile_last from the wrfout
  except:
    sys.exit("STOP! There is a GAP between the record of the last day %s and earlieast wrfout we have in this folder"%(lastday))
  return ncfile_last,lastindex

def processnewwrfout(rawnc,periods,taskname,nstep,filenames,var_parameters ,
               wrfncfile_last ,shiftday,var_units,var_description,
               wrfinputnc,outputdim,z_levs,number_of_zlevs, 
               units_cur,calendar_cur,compute_mode,
               nx,ny,nlev,lastindex):

  if periods=="daily":

    ppdaily(rawnc,nstep,taskname,filenames,var_parameters ,
               wrfncfile_last ,shiftday,var_units,var_description,
               wrfinputnc,outputdim,z_levs,number_of_zlevs,
               units_cur,calendar_cur,compute_mode,
               nx,ny,nlev,lastindex)

  elif periods=="hourly":
    pphourly(rawnc,nstep,taskname,filenames,var_parameters ,
               wrfncfile_last ,shiftday,var_units,var_description,
               wrfinputnc,outputdim,z_levs,number_of_zlevs,
               units_cur,calendar_cur,compute_mode,
               nx,ny,nlev,lastindex)
  else:
    sys.exit("no such option, please choose periods from daily and hourly")
  return rawnc

def ppdaily(rawnc,nstep,taskname,filenames,var_parameters ,
               wrfncfile_last ,shiftday,var_units,var_description,
               wrfinputnc,outputdim,z_levs,number_of_zlevs, 
               units_cur,calendar_cur,compute_mode,
               nx,ny,nlev,lastindex):
    from process import anal_daily, wrftimetodate
    import sys
    from netCDF4 import date2num,num2date
    from netCDF4 import Dataset
    import numpy as np
    from constant import mmstommday
    nt=len(filenames)-shiftday 
    outputtime=np.empty([1])
    outputdata={}
    for field in var_parameters[taskname]["fields"]:
      outputdata[field]=np.empty([nlev,ny,nx])
    simbeg_date=wrftimetodate(Dataset(filenames[shiftday],'r').variables['Times'][0])
    simbeg_num =date2num( simbeg_date,units=units_cur,calendar=calendar_cur)
    for iday,filename in enumerate(filenames[shiftday:]):
      wrfncfile_cur=openwrfdata(filename)
      ntime       =wrfncfile_cur.dimensions['Time'].size
      curtime=wrfncfile_cur.variables['Times']
      date_curstep=wrftimetodate(curtime[0])
      outputtime[0]=date2num( wrftimetodate(curtime[0]),units=units_cur,calendar=calendar_cur)
      # check of the intergrty of data
      timenum=iday+simbeg_num
      chekitegrty(curtime,nstep,filename,filenames,outputtime,timenum)
      # check of the intergrty of data

      if compute_mode==6:
        if taskname=="PR":
          outputdata["PRAVG"][:,:]=np.sum(wrfncfile_cur.variables['PRAVG'][:,:,:],axis=0)*1.0/ntime*mmstommday
          """
          temp_1=wrfncfile_cur.variables['RAINC'][0,:,:] -wrfncfile_last.variables['RAINC'][0,:,:]
          temp_2=wrfncfile_cur.variables['RAINNC'][0,:,:]-wrfncfile_last.variables['RAINNC'][0,:,:]
          if np.all(temp_1>=0) and np.all(temp_1>=0):
            outputdata["PRAVG"][:,:]=temp_1+temp_2
          else:
            outputdata["PRAVG"][:,:]=0.0
          """
        else:
          for field in var_parameters[taskname]["fields"]:
            outputdata[field][:,:]=wrfncfile_cur.variables[taskname][0,:,:]-wrfncfile_last.variables[taskname][0,:,:]
      else:
        anal_daily(iday,outputdata,wrfncfile_cur,wrfinputnc,taskname,
                  var_parameters[taskname]["fields"],var_parameters[taskname]["vert_intp"],outputdim,z_levs,number_of_zlevs,compute_mode)
      wrfncfile_last=wrfncfile_cur
      for field in var_parameters[taskname]["fields"]:
        if outputdim==3:
          rawnc.variables[field][lastindex+iday,:,:]=outputdata[field]
        elif outputdim==4:
          rawnc.variables[field][lastindex+iday,:,:,:]=outputdata[field]
      rawnc.variables["time"][lastindex+iday]=outputtime[0]
      rawnc.sync()
      print(date_curstep)
    return rawnc

def pphourly(rawnc,nstep,taskname,filenames,var_parameters ,
               wrfncfile_last ,shiftday,var_units,var_description,
               wrfinputnc,outputdim,z_levs,number_of_zlevs, 
               units_cur,calendar_cur,compute_mode,
               nx,ny,nlev,lastindex):

    from netCDF4 import Dataset
    from process import anal_hourly, wrftimetodate
    from netCDF4 import date2num,num2date
    import numpy as np
    import sys
    nt=nstep*len(filenames)
    outputtime=np.empty([nstep])
    outputdata={}
    for field in var_parameters[taskname]["fields"]:
      if nlev==1:
        outputdata[field]=np.empty([nstep,ny,nx])
      else:
        outputdata[field]=np.empty([nstep,nlev,ny,nx])
    simbeg_date=wrftimetodate(Dataset(filenames[shiftday],'r').variables['Times'][0])
    simbeg_num =date2num( simbeg_date,units=units_cur,calendar=calendar_cur)
    for iday,filename in enumerate(filenames[shiftday:]):
      wrfncfile_cur=openwrfdata(filename)
      curtime=wrfncfile_cur.variables['Times']
      date_curstep=wrftimetodate(curtime[0])
      for istep in range(nstep):
        outputtime[istep]=date2num( wrftimetodate(curtime[istep]),units=units_cur,calendar=calendar_cur)

      timenum=iday*24+simbeg_num #watchout default the wrfout should be daily!
      chekitegrty(curtime,nstep,filename,filenames,outputtime,timenum)
      anal_hourly(iday,outputdata,wrfncfile_cur,wrfinputnc,taskname,
                  var_parameters[taskname]["fields"],var_parameters[taskname]["vert_intp"],outputdim,z_levs,number_of_zlevs,compute_mode)
      #for field in var_parameters[taskname]["fields"]:
      #  outputdata[field][:]=wrfncfile_cur.variables[taskname][:]
      #wrfncfile_last=wrfncfile_cur
      for field in var_parameters[taskname]["fields"]:
        if outputdim==3:
          rawnc.variables[field][lastindex+iday*nstep:,:,:]=outputdata[field]
        elif outputdim==4:
          rawnc.variables[field][lastindex+iday*nstep:,:,:,:]=outputdata[field]
      rawnc.variables["time"][lastindex+iday*nstep:]=outputtime[:]
      rawnc.sync()
      print(date_curstep)
    return rawnc



 
def fromwrfout(filenames,rawfname,casename,taskname,periods,
               units_cur,calendar_cur,var_parameters,nx,ny,nlev,nstep,compute_mode,
               wrfncfile_last ,shiftday,var_units,var_description,
               wrfinputnc,outputdim,z_levs,number_of_zlevs):
  from netCDF4 import Dataset
  import sys
  import numpy as np
  rawnc,raw_exit=getrawnc(rawfname,casename,taskname,periods,units_cur,calendar_cur,
                        var_parameters,nx,ny,nlev,var_units,var_description )
   
  lastindex=0
  if raw_exit:
    wrfncfile_last,lastindex=updatelastwrfoutfn(rawnc,units_cur,calendar_cur,filenames)

  if len(filenames)>shiftday:
    processnewwrfout(rawnc,periods,taskname,nstep,filenames,var_parameters ,
               wrfncfile_last ,shiftday,var_units,var_description,
               wrfinputnc,outputdim,z_levs,number_of_zlevs,
               units_cur,calendar_cur,compute_mode,
               nx,ny,nlev,lastindex)
  else:
    print("no new wrfout to process")

  return rawnc

