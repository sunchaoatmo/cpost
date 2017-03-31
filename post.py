#!/usr/bin/env python
import subprocess
import numpy as np
from writenc import createnc
from netCDF4 import Dataset
from toolkit import diag_season_month, wrftimetodate
from netCDF4 import date2num,num2date
from datetime import datetime,timedelta
from constant import *   # only several constants
from POSTparameter import postList,var_parameters
from argument import args
import glob
import os.path

# This code assumes each wrfout file contains one full day of output 
# and that all data is consecutive (no gaps)


MATCHINE=subprocess.check_output("uname -a" , shell=True)
calendar_cur=args.calendar
periods=args.p[0]
casename=args.n
vnames=args.v if args.v else postList[args.p]

if periods=="daily":
  units_cur = 'days since 0001-01-01 00:00'
elif periods=="hourly":
  units_cur = 'hours since 0001-01-01 00:00'
Oneday=timedelta(days=1)
filenames_all=sorted(glob.glob(prefix))
"""
filesize0=os.path.getsize(filenames_all[0])
filenames=[]
for filename in filenames_all:
  if not filesize0==os.path.getsize(filename):
    break
  filenames.append(filename)
"""
filenames=filenames_all

if args.sjob:
  nprocs=len(vnames)
  if 'deepthought2' in MATCHINE: 
    nodes=nprocs/20+1
  else:
    nodes=nprocs/32+1
  job_t="job_temp_%s.pbs"%str(args.p[0])
  with open('job_default.pbs', 'r') as fin:
    with open(job_t, 'w') as fout:
       for line in fin:
         line=line.replace("NPROCS",str(nprocs))
         line=line.replace("PERIOD",str(args.p[0]))
         line=line.replace("VAR",','.join(vnames))
         fout.write(line)
  cmd="qsub "+job_t
  subprocess.call(cmd,shell=True)
else:
  if args.mpi:
    import imp
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank = comm.Get_rank()
    if nprocs>=len(vnames):
      try:
        var_loc =vnames[rank]
        cmd="post.py -p "+args.p[0] +" -v "+var_loc
        print(cmd,rank)
        subprocess.call(cmd,shell=True)
        print("call finished on case %s " % (cmd))
      except:
        import sys, traceback
        traceback.print_exc(file=sys.stdout)
        sys.exit("ERROR in %s rank %s"%(cmd,rank))
    else:
      print("we don't have enough CPU %s there are %s tasks %s cases %s vars %s years"
                                 % (nprocs,tot_year*tot_case*tot_postvar, tot_case,tot_postvar,tot_year))
    comm.Barrier()
  else:
    for vname in vnames:
      shiftday=var_parameters[vname]['shiftday'] if periods=="daily"  else 0
      compute_mode=var_parameters[vname]['compute_mode'] 
      rawfname="%s_%s_%s.nc"%(casename,vname,periods)
      ncexist=os.path.isfile(rawfname)
      ncfile_last=Dataset(filenames[0],'r')
      var_units  =ncfile_last.variables[vname].units
      var_description  =ncfile_last.variables[vname].description
      var_shape=ncfile_last.variables[vname].shape
      outputdim=len(var_shape)
      if outputdim==3:
        (nstep,nx,ny)=ncfile_last.variables[vname].shape
      elif outputdim==4:
        (nstep,nlev,nx,ny)=ncfile_last.variables[vname].shape
      lastindex=0
      if ncexist:
        rawnc=Dataset(rawfname,'a')
        lastday=num2date(rawnc.variables["time"][-1],units=units_cur,calendar=calendar_cur)
        lastwrfout="wrfout_d01_%s"%lastday.strftime(wrfout_data_fmt)
        try:
          lastindex=filenames.index(lastwrfout)
          del filenames[:lastindex]
        except:
          import sys
          sys.exit("STOP! There is a GAP between the record of the last day %s and earlieast wrfout we have in this folder"%(lastday))
      else:
        rawnc=createnc(casename,vname,periods,units_cur,calendar_cur,nx,ny,dimension=outputdim)
        rawnc.variables[vname].units=var_units
        rawnc.variables[vname].description=var_description
      if len(filenames)>shiftday:
        if periods=="daily":
          outputdata=np.empty([len(filenames)-shiftday,nx,ny])
          outputtime=np.empty([len(filenames)-shiftday])
        else:
          outputdata=np.empty([nstep*len(filenames)-shiftday,nx,ny])
          outputtime=np.empty([nstep*len(filenames)-shiftday])
        simbeg_date=wrftimetodate(Dataset(filenames[shiftday],'r').variables['Times'][0])
        for iday,filename in enumerate(filenames[shiftday:]):
          ncfile_cur=Dataset(filename,'r')
          curtime=ncfile_cur.variables['Times']
          date_curstep=wrftimetodate(curtime[0])
          if date_curstep==iday*Oneday+simbeg_date:
            if len(curtime)==nstep:
              if periods=="daily":
                if compute_mode==6:
                  if vname=="PRAVG":
                    outputdata[iday,:,:]=(ncfile_cur.variables['RAINC'][0,:,:]-ncfile_last.variables['RAINC'][0,:,:]
                                         +ncfile_cur.variables['RAINNC'][0,:,:]-ncfile_last.variables['RAINNC'][0,:,:])
                  else:
                    outputdata[iday,:,:]=ncfile_cur.variables[vname][0,:,:]-ncfile_last.variables[vname][0,:,:]
                elif compute_mode==1:
                  outputdata[iday,:,:]=np.mean(ncfile_cur.variables[vname][:,:,:],axis=0)
                outputtime[iday]=date2num( date_curstep,units=units_cur,calendar=calendar_cur)
              else:
                outputdata[iday*nstep:(iday+1)*nstep,:,:]=ncfile_cur.variables[vname][:,:,:]
                for istep in range(nstep):
                  outputtime[iday*nstep+istep]=date2num( wrftimetodate(curtime[istep]),units=units_cur,calendar=calendar_cur)
              print(date_curstep)
            else:
              if filename == filenames[-1]:
                break
              else:
                import sys
                sys.exit("STOP! one wrfout is incomplete %s ",(filename))
          else:
            import sys
            sys.exit("STOP! one day is missing bettween output %s and %s",(filenames[iday],filenames[iday-1]))
          ncfile_last=ncfile_cur

        if outputdim==3:
          rawnc.variables[vname][lastindex:,:,:]=outputdata

        rawnc.variables["time"][lastindex:]=outputtime

########################DIAG PART#################################
      if periods=="daily":
        if vname=="PRAVG":
          postlist=["PCT","PRAVG","RAINYDAYS","R10","R5D","CDD","R95T","SDII"]
        else:
          postlist=[vname]
        nctime=rawnc.variables["time"]
        start_ymd=num2date(nctime[0],units=units_cur,calendar=calendar_cur)
        end_ymd  =num2date(nctime[-1],units=units_cur,calendar=calendar_cur)
        diag_season_month("seasonal",rawnc,seasonList,start_ymd,end_ymd,postlist,
                          vname,casename,shiftday,calendar_cur,units_cur,var_units,var_description,nx,ny)
        diag_season_month("monthly",rawnc,monthlyList,start_ymd,end_ymd,postlist,
                          vname,casename,shiftday,calendar_cur,units_cur,var_units,var_description,nx,ny)
      rawnc.close() #flush out rawnc
