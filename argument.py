#!/usr/bin/env python
from os import getcwd
import argparse
parser = argparse.ArgumentParser(description='Post pocess the varialbe you choose ')
parser.add_argument("-p",help="run type",default='daily',choices=['daily','hourly'])
parser.add_argument("--calendar",help="the type of calendar which is gonna be used",default='gregorian',
                    choices=['gregorian','noleap'])
parser.add_argument("-i",help="wrfinput location")
parser.add_argument("-v",help="post pocess varialble")
parser.add_argument("-b2n",help="whether convert binary to NetCDF, if true, not post pocess" ,action="store_true" )
parser.add_argument("-exe",help="location to the post.exe")
parser.add_argument("-n",default=getcwd().split("/")[-1],help="name of the case")
parser.add_argument("-mpi", help="MPI flag, determine whether run the script in MPI mode", action="store_true")
parser.add_argument("--sjob",help='create and submit jobs',action="store_true")

args = parser.parse_args() 
if args.v:
  args.v=args.v.split(",")
if args.p:
  args.p=args.p.split(",")
