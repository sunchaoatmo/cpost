import numpy as np
seasonList=[ [12,  2], [3, 5], [6,  8], [9,  11]]
monthlyList=[ [12, 1], [1, 2],[2,3],[3,4],[4,5],[5,6],[6, 7],[7, 8],[8,9], [9, 10],[10, 11],[11,12]]
wrfout_data_fmt="%Y-%m-%d_%H:%M:%S"
prefix="wrfout*"
shiftday=1
computemode=2
outputdim=3
dry_lim=1
pct=0.95
G = 9.81
Rd = 287.04
Rv = 461.6
Rm = .608 
Cp = 1004.
Cp = 7.*Rd/2.
Cv = Cp-Rd
CPMD = 0.887
RCP = Rd/Cp
p0 = 100000.
# be awared the z_levs should be sorted from large to small
z_levs=[1000.0, 975.0, 950.0, 925.0, 900.0, 875.0, 850.0, 825.0, 800.0, 775.0, 750.0, 700.0, 650.0, 600.0, 550.0, 500.0, 450.0, 400.0, 350.0, 300.0, 250.0, 225.0, 200.0, 175.0, 150.0, 125.0, 100.0, 70.0]
z_levs=np.array(z_levs)
z_levs=100*z_levs
number_of_zlevs  =len(z_levs)

