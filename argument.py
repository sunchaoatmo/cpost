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
    print("Sorry can not open %s"% filepath)
    sys.exit("please provide correct the path to r95t_hist file in argument")
  return filenc

def saveenv(file_hd,field):
  userpath = raw_input('Enter your %s file location: '%field)
  while not path.isfile(userpath):
     print("There is no file in %s" %userpath)
     userpath = raw_input('Please enter CORRECTED  %s file location: '%field)
  file_hd.write(field+"="+userpath+"\n")
  return

parser = argparse.ArgumentParser(description='Post pocess the varialbe you choose ')
parser.add_argument("-p",help="run type",default='daily',choices=['daily','hourly'])
parser.add_argument("--calendar",help="the type of calendar which is gonna be used",default='gregorian',
                    choices=['gregorian','noleap'])
parser.add_argument("-i",help="wrfinput location")
parser.add_argument("-v",help="post pocess varialble")
parser.add_argument("-exe",help="location to the post.exe")
parser.add_argument("-n",default=getcwd().split("/")[-1],help="name of the case")
parser.add_argument("-mpi", help="MPI flag, determine whether run the script in MPI mode", action="store_true")
parser.add_argument("--sjob","-s",help='create and submit jobs',action="store_true")
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
  print "Welcome, this might be your first time to use cpost, we need to do some env setting:"
  target = open(configfile, 'w')
  target.write("[PATH]\n")
  saveenv(target,"r95t")
  saveenv(target,"wrfinput")
  target.close()
  print "Finished init the configuration file... Thanks! CS"

config = ConfigParser.ConfigParser()
config.read(configfile)
r95t=config.get('PATH','r95t')
wrfinput=config.get('PATH','wrfinput')
r95tnc    =loadfile(r95t)
wrfinputnc=loadfile(wrfinput)
