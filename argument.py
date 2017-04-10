#!/usr/bin/env python
import ConfigParser 
from os import getcwd,path,getenv
from netCDF4 import Dataset
import argparse
import argparse
def loadfile(filepath):
  if path.isfile(filepath):
    filenc=Dataset(filepath, 'r')
  else:
    print("Sorry can not open %s"% r95t_hist)
    sys.exit("please provide correct the path to r95t_hist file in argument")
  return filenc


parser = argparse.ArgumentParser(description='Post pocess the varialbe you choose ')
parser.add_argument("-p",help="run type",default='daily',choices=['daily','hourly'])
parser.add_argument("--calendar",help="the type of calendar which is gonna be used",default='gregorian',
                    choices=['gregorian','noleap'])
parser.add_argument("-i",help="wrfinput location")
parser.add_argument("-v",help="post pocess varialble")
parser.add_argument("-exe",help="location to the post.exe")
parser.add_argument("-n",default=getcwd().split("/")[-1],help="name of the case")
parser.add_argument("-mpi", help="MPI flag, determine whether run the script in MPI mode", action="store_true")
#parser.add_argument("--r95t","-r", help="location of your r95t_hist file", default= "/homes/sunchao/data//OBS_R95T_hist_seasonal.nc" )
#parser.add_argument("--wrfinput","-w", help="location of your wrfinput file", default= "/homes/sunchao/lustre/US/icbc_1999/wrfinput_d01.1998120100" )
parser.add_argument("--sjob",help='create and submit jobs',action="store_true")
z_levs=[1000.0, 975.0, 950.0, 925.0, 900.0, 875.0, 850.0, 825.0, 800.0, 775.0, 750.0, 700.0, 650.0, 600.0, 550.0, 500.0, 450.0, 400.0, 350.0, 300.0, 250.0, 225.0, 200.0, 175.0, 150.0, 125.0, 100.0, 70.0]
z_levs=[x*100 for x in z_levs]
parser.add_argument('--z_levs', '-z', nargs='+',help='level to interpolate',default=z_levs, type=int)

args = parser.parse_args() 
if args.v:
  args.v=args.v.split(",")
if args.p:
  args.p=args.p.split(",")
configfile=getenv("HOME")+'/.post.ini'

if not path.isfile(configfile):
  print "This is your first time to use cpost, we need to do some env setting:"
  target = open(configfile, 'w')
  userpath = raw_input('Enter your r95t file location: ')
  target.write("[PATH]\n")
  target.write("r95t="+userpath+"\n")
  userpath = raw_input('Enter your wrfinput file location: ')
  target.write("wrfinput="+userpath+"\n")
  target.close()
  print "Finished init the configuration file... Thanks! CS"

config = ConfigParser.ConfigParser()
config.read(configfile)
r95t=config.get('PATH','r95t')
wrfinput=config.get('PATH','wrfinput')
r95tnc    =loadfile(r95t)
wrfinputnc=loadfile(wrfinput)
