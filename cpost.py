#!/usr/bin/env python
import numpy as np
import sys
from writenc import createnc
from netCDF4 import Dataset
from process import anal_daily,anal_sea_mon, wrftimetodate
from netCDF4 import date2num,num2date
from datetime import datetime,timedelta
from constant import *   # only several constants
from POSTparameter import daily_allvars_2d,daily_allvars_3d,var_parameters
from argument import args,wrfinputnc,r95tnc
import glob
import os.path
import time
from subprocess import check_output,call

# This code assumes each wrfout file contains one full day of output 
# and that all data is consecutive (no gaps)

def setmetadata(filename, rawdata):
  ncfile_last=Dataset(filename,'r')
  var_units,var_description={},{}
  nx,ny,nstep=0,0,0
  outputdim=3
  if rawdata:
    if taskname in ncfile_last.variables:
      var_units[taskname]  =ncfile_last.variables[taskname].units
      var_description[taskname]  =ncfile_last.variables[taskname].description
      varshape =ncfile_last.variables[taskname].shape
      nz       =varshape[1]
      outputdim=len(varshape)
    else:
      outputdim=var_parameters[taskname]['dim']
      for field in var_parameters[taskname]["fields"]:
        var_units[field]=var_parameters[taskname]["fields"][field]['units']
        var_description[field]=var_parameters[taskname]["fields"][field]['description']
      nz=ncfile_last.dimensions['bottom_top'].size
    nstep=ncfile_last.dimensions['Time'].size
    ny=ncfile_last.dimensions['south_north'].size
    nx=ncfile_last.dimensions['west_east'].size
    if outputdim==3:
      nlev=1
    elif outputdim==4:
      nlev =  number_of_zlevs if var_parameters[taskname]['vert_intp'] else nz 
  else:
    for field in var_parameters[taskname]["fields"]:
      var_units[field]=var_parameters[taskname]["fields"][field]['units']
      var_description[field]=var_parameters[taskname]["fields"][field]['description']
    nz=1
    nstep=ncfile_last.dimensions['time'].size
    try:
      ny=ncfile_last.dimensions['lat'].size
      nx=ncfile_last.dimensions['lon'].size
    except:
      ny=ncfile_last.dimensions['south_north'].size
      nx=ncfile_last.dimensions['west_east'].size

    nlev =  number_of_zlevs if var_parameters[taskname]['vert_intp'] else nz 
  return nx,ny,nz,nlev,var_units,var_description,nstep,outputdim,ncfile_last



MATCHINE=check_output("uname -a" , shell=True)
calendar_cur=args.calendar
periods=args.p
casename=args.n
if args.v:
  tasknames=args.v
elif args.dim=='2':
  tasknames=daily_allvars_2d
elif args.dim=='3':
  tasknames=daily_allvars_3d
z_levs=args.z_levs
number_of_zlevs=len(z_levs)

if periods=="daily":
  units_cur = 'days since 0001-01-01 00:00'
elif periods=="hourly":
  units_cur = 'hours since 0001-01-01 00:00'
Oneday=timedelta(days=1)
filenames=sorted(glob.glob(prefix))

if args.sjob:
  nprocs=len(tasknames)
  if 'deepthought2' in MATCHINE: 
    nodes=nprocs/20+1
  else:
    nodes=nprocs/32+1
  job_t="job_temp_%s.pbs"%str(args.p)
  with open('job_default.pbs', 'r') as fin:
    with open(job_t, 'w') as fout:
       for line in fin:
         line=line.replace("NPROCS",str(nprocs))
         line=line.replace("PERIOD",str(args.p))
         line=line.replace("VAR",','.join(tasknames))
         fout.write(line)
  cmd="qsub "+job_t
  call(cmd,shell=True)
else:
  if args.mpi:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank = comm.Get_rank()
    print("my rank %s"%rank)
    if nprocs>=len(tasknames):
      try:
        var_loc =tasknames[rank]
        cmd="cpost.py  -p "+args.p +" -v "+var_loc
        print(cmd,rank)
        call(cmd,shell=True)
        print("call finished on case %s " % (cmd))
      except:
        var_loc =tasknames[rank]
        print("ERROR in rank %s with var: %s"%(rank,var_loc))
    else:
      print("we don't have enough CPU ")
    comm.Barrier()
  else:
    for taskname in tasknames:
      shiftday=var_parameters[taskname]['shiftday'] if periods=="daily"  else 0
      compute_mode=var_parameters[taskname]['compute_mode'] 
      rawfname="%s_%s_%s.nc"%(casename,taskname,periods)
      if filenames:
        filename= filenames[0]
        rawdata = True
      else:
        filename=rawfname
        rawdata=False

      nx,ny,nz,nlev,var_units,var_description,nstep,outputdim,ncfile_last=setmetadata(filename,rawdata)


      lastindex=0


      if rawdata:
        try:
          rawnc=Dataset(rawfname,'a')
          lastday=num2date(rawnc.variables["time"][-1],units=units_cur,calendar=calendar_cur)
          lastwrfout="wrfout_d01_%s"%lastday.strftime(wrfout_data_fmt)
          try:
            lastindex=filenames.index(lastwrfout)
            del filenames[:lastindex]
            ncfile_last=Dataset(filenames[0],'r')
          except:
            sys.exit("STOP! There is a GAP between the record of the last day %s and earlieast wrfout we have in this folder"%(lastday))
        except:
          if os.path.isfile(rawfname):
            os.remove(rawfname)
          rawnc=createnc(casename,taskname,periods,units_cur,calendar_cur,var_parameters[taskname]["fields"].keys(),nx,ny,nlev  )
          for field in var_parameters[taskname]["fields"]:
            rawnc.variables[field].units=var_units[field]
            rawnc.variables[field].description=var_description[field]
        if len(filenames)>shiftday:
          nt=len(filenames)-shiftday if periods=="daily" else nstep*len(filenames)-shiftday
          outputtime=np.empty([1])
          outputdata={}
          for field in var_parameters[taskname]["fields"]:
            outputdata[field]=np.empty([nlev,ny,nx])
          simbeg_date=wrftimetodate(Dataset(filenames[shiftday],'r').variables['Times'][0])
          simbeg_num =date2num( simbeg_date,units=units_cur,calendar=calendar_cur)
          for iday,filename in enumerate(filenames[shiftday:]):
            try:
              ncfile_cur=Dataset(filename,'r')
            except:
              sys.exit("can not open:%s "%filename)
            curtime=ncfile_cur.variables['Times']
            date_curstep=wrftimetodate(curtime[0])
            if periods=="daily":
              outputtime[0]=date2num( wrftimetodate(curtime[0]),units=units_cur,calendar=calendar_cur)
            else:
              for istep in range(nstep):
                outputtime[iday*nstep+istep]=date2num( wrftimetodate(curtime[istep]),units=units_cur,calendar=calendar_cur)

# check of the intergrty of data
            if not len(curtime)==nstep:
              if filename == filenames[-1]:     # it is OK sometimes the last wrfout can be incomplete, the simulation can restart from it later on
                break
              else:
                import sys
                sys.exit("STOP! one wrfout is incomplete %s "%(filename))

            if not outputtime[0]==iday+simbeg_num:
              import sys
              sys.exit("STOP! %s is missing in wrfout serial "%(filename))
# check of the intergrty of data

            if periods=="daily":
              if compute_mode==6:
                if taskname=="PR":
                  temp_1=ncfile_cur.variables['RAINC'][0,:,:] -ncfile_last.variables['RAINC'][0,:,:]
                  temp_2=ncfile_cur.variables['RAINNC'][0,:,:]-ncfile_last.variables['RAINNC'][0,:,:]
                  if np.all(temp_1>=0) and np.all(temp_1>=0):
                    outputdata["PRAVG"][:,:]=temp_1+temp_2
                  else:
                    outputdata["PRAVG"][:,:]=0.0
                else:
                  for field in var_parameters[taskname]["fields"]:
                    outputdata[field][:,:]=ncfile_cur.variables[taskname][0,:,:]-ncfile_last.variables[taskname][0,:,:]
              else:
                anal_daily(iday,outputdata,ncfile_cur,wrfinputnc,taskname,
                          var_parameters[taskname]["fields"],var_parameters[taskname]["vert_intp"],outputdim,z_levs,number_of_zlevs,compute_mode)
            else:
              for field in var_parameters[taskname]["fields"]:
                outputdata[field][iday*nstep:(iday+1)*nstep,:,:]=ncfile_cur.variables[taskname][:,:,:]
            ncfile_last=ncfile_cur
            for field in var_parameters[taskname]["fields"]:
              if outputdim==3:
                rawnc.variables[field][lastindex+iday,:,:]=outputdata[field]
              elif outputdim==4:
                rawnc.variables[field][lastindex+iday,:,:,:]=outputdata[field]
            rawnc.variables["time"][lastindex+iday]=outputtime[0]
            rawnc.sync()
            print(date_curstep)
      else:
        rawnc=Dataset(filename,'r')

########################DIAG PART#################################
      if periods=="daily":
        anal_sea_mon("seasonal",rawnc,seasonList,var_parameters[taskname]["fields"].keys(),
                          taskname,casename,shiftday,calendar_cur,units_cur,var_units,var_description,ny,nx,nlev,r95tnc)
        anal_sea_mon("monthly",rawnc,monthlyList,var_parameters[taskname]["fields"].keys(),
                          taskname,casename,shiftday,calendar_cur,units_cur,var_units,var_description,ny,nx,nlev,r95tnc)
      if rawdata:
        rawnc.history = 'Post-processed Chao Sun sunchao@umd.edu ' + time.ctime(time.time())
        rawnc.source ='CWRF'
      rawnc.close() #flush out rawnc
