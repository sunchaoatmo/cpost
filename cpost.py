#!/usr/bin/env python
import numpy as np
from writenc import createnc
from netCDF4 import Dataset
from process import anal_daily,anal_sea_mon, wrftimetodate
from netCDF4 import date2num,num2date
from datetime import datetime,timedelta
from constant import *   # only several constants
from POSTparameter import postList,var_parameters
from argument import args,wrfinputnc,r95tnc
import glob
import os.path
import time
from subprocess import check_output,call

# This code assumes each wrfout file contains one full day of output 
# and that all data is consecutive (no gaps)

MATCHINE=check_output("uname -a" , shell=True)
calendar_cur=args.calendar
periods=args.p
casename=args.n
tasknames=args.v if args.v else postList[periods]
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
  job_t="job_temp_%s.pbs"%str(args.p[0])
  with open('job_default.pbs', 'r') as fin:
    with open(job_t, 'w') as fout:
       for line in fin:
         line=line.replace("NPROCS",str(nprocs))
         line=line.replace("PERIOD",str(args.p[0]))
         line=line.replace("VAR",','.join(tasknames))
         fout.write(line)
  cmd="qsub "+job_t
  call(cmd,shell=True)
else:
  if args.mpi:
    import imp
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank = comm.Get_rank()
    if nprocs>=len(tasknames):
      try:
        var_loc =tasknames[rank]
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
    for taskname in tasknames:
      shiftday=var_parameters[taskname]['shiftday'] if periods=="daily"  else 0
      compute_mode=var_parameters[taskname]['compute_mode'] 
      rawfname="%s_%s_%s.nc"%(casename,taskname,periods)
      ncexist=os.path.isfile(rawfname)
      ncfile_last=Dataset(filenames[0],'r')
      var_units,var_description={},{}
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
      lastindex=0


      if ncexist:
        rawnc=Dataset(rawfname,'a')
        lastday=num2date(rawnc.variables["time"][-1],units=units_cur,calendar=calendar_cur)
        lastwrfout="wrfout_d01_%s"%lastday.strftime(wrfout_data_fmt)
        try:
          lastindex=filenames.index(lastwrfout)
          del filenames[:lastindex]
          ncfile_last=Dataset(filenames[0],'r')
        except:
          import sys
          sys.exit("STOP! There is a GAP between the record of the last day %s and earlieast wrfout we have in this folder"%(lastday))
      else:
        rawnc=createnc(casename,taskname,periods,units_cur,calendar_cur,var_parameters[taskname]["fields"].keys(),nx,ny,nlev  )
        for field in var_parameters[taskname]["fields"]:
          rawnc.variables[field].units=var_units[field]
          rawnc.variables[field].description=var_description[field]
      if len(filenames)>shiftday:
        nt=len(filenames)-shiftday if periods=="daily" else nstep*len(filenames)-shiftday
        outputtime=np.empty([nt])
        outputdata={}
        for field in var_parameters[taskname]["fields"]:
          outputdata[field]=np.empty([nt,nlev,ny,nx])
        simbeg_date=wrftimetodate(Dataset(filenames[shiftday],'r').variables['Times'][0])
        for iday,filename in enumerate(filenames[shiftday:]):
          ncfile_cur=Dataset(filename,'r')
          curtime=ncfile_cur.variables['Times']
          date_curstep=wrftimetodate(curtime[0])
          if date_curstep==iday*Oneday+simbeg_date:
            if len(curtime)==nstep:
              if periods=="daily":
                if compute_mode==6:
                  if taskname=="PR":
                    for field in var_parameters[taskname]["fields"]:
                      if field =="PRAVG":
                        outputdata[field][iday,:,:]=(ncfile_cur.variables['RAINC'][0,:,:]-ncfile_last.variables['RAINC'][0,:,:]
                                           +ncfile_cur.variables['RAINNC'][0,:,:]-ncfile_last.variables['RAINNC'][0,:,:])
                  else:
                    for field in var_parameters[taskname]["fields"]:
                      outputdata[field][iday,:,:]=ncfile_cur.variables[taskname][0,:,:]-ncfile_last.variables[taskname][0,:,:]
                elif compute_mode==1:
                  anal_daily(iday,outputdata,ncfile_cur,wrfinputnc,taskname,
                            var_parameters[taskname]["fields"],var_parameters[taskname]["vert_intp"],outputdim,z_levs,number_of_zlevs)
                outputtime[iday]=date2num( date_curstep,units=units_cur,calendar=calendar_cur)
              else:
                for field in var_parameters[taskname]["fields"]:
                  outputdata[field][iday*nstep:(iday+1)*nstep,:,:]=ncfile_cur.variables[taskname][:,:,:]
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
            print(date_curstep)
            import sys
            sys.exit("STOP! one day is missing bettween output %s and %s",(filenames[iday],filenames[iday-1]))
          ncfile_last=ncfile_cur

        for field in var_parameters[taskname]["fields"]:
          if outputdim==3:
            rawnc.variables[field][lastindex:,:,:]=outputdata[field]
          elif outputdim==4:
            rawnc.variables[field][lastindex:,:,:,:]=outputdata[field]

        rawnc.variables["time"][lastindex:]=outputtime

########################DIAG PART#################################
      if periods=="daily":
        nctime=rawnc.variables["time"]
        start_ymd=num2date(nctime[0],units=units_cur,calendar=calendar_cur)
        end_ymd  =num2date(nctime[-1],units=units_cur,calendar=calendar_cur)
        anal_sea_mon("seasonal",rawnc,seasonList,start_ymd,end_ymd,var_parameters[taskname]["fields"].keys(),
                          taskname,casename,shiftday,calendar_cur,units_cur,var_units,var_description,ny,nx,nlev,r95tnc)
        anal_sea_mon("monthly",rawnc,monthlyList,start_ymd,end_ymd,var_parameters[taskname]["fields"].keys(),
                          taskname,casename,shiftday,calendar_cur,units_cur,var_units,var_description,ny,nx,nlev,r95tnc)
      rawnc.close() #flush out rawnc


