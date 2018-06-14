#!/usr/bin/env python
from netCDF4 import Dataset
from process import anal_sea_mon #, wrftimetodate
from datetime import datetime,timedelta
from POSTparameter import hourly_allvars,daily_allvars_2d,daily_allvars_3d,var_parameters
from argument import args,wrfinputnc,r95tnc
from constant import prefix,seasonList,monthlyList
import glob
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
if args.p=="hourly":
  tasknames=hourly_allvars
if args.v:
  tasknames=args.v
elif args.dim=='2':
  tasknames=daily_allvars_2d
elif args.dim=='3':
  tasknames=daily_allvars_3d
z_levs=args.z_levs
number_of_zlevs=len(z_levs)

if periods=="daily":
  units_cur = 'days since 0001-01-01 00:00:00'
elif periods=="hourly":
  units_cur = 'hours since 0001-01-01 00:00'
Oneday=timedelta(days=1)
if args.inputfname:
  filenames=[args.inputfname]
else:
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
      if rank<len(tasknames):
        var_loc =tasknames[rank]
        cmd="cpost.py -p %s -c %s -v %s"%(args.p,args.calendar,var_loc)
        print(cmd,rank)
        call(cmd,shell=True)
        print("call finished on case %s " % (cmd))
    else:
      print("we don't have enough CPU ")
    comm.Barrier()
  else:
    for taskname in tasknames:
      shiftday=var_parameters[taskname]['shiftday'] if periods=="daily"  else 0
      compute_mode=var_parameters[taskname]['compute_mode'] 
      rawfname="%s_%s_%s.nc"%(casename,taskname,periods)
      if filenames:
        filename= filenames[-1]
        rawdata = True
      else:
        filename=rawfname
        rawdata=False

      nx,ny,nz,nlev,var_units,var_description,nstep,outputdim,ncfile_last=setmetadata(filename,rawdata)




      if rawdata:
        from readin import fromwrfout
        rawnc=fromwrfout(filenames,rawfname,casename,taskname,periods,
               units_cur,calendar_cur,var_parameters,nx,ny,nlev,nstep,compute_mode,
               ncfile_last ,shiftday,var_units,var_description,
               wrfinputnc,outputdim,z_levs,number_of_zlevs)
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
