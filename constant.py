seasonList=[ [12,  2], [3, 5], [6,  8], [9,  11]]
opt_thresh=1.0
fill_nocloud=1
missing     =-999
monthlyList=[ [1, 1],[2,2],[3,3],[4,4],[5,5],[6, 6],[7, 7],[8,8], [9, 9 ],[10, 10],[11,11],[12,12]]
mmstommday=24.0*60.*60.
wrfout_data_fmt="%Y-%m-%d"
#wrfout_data_fmt="%Y-%m-%d_%H:%M:%S"
prefix="wrfout*"
dry_lim=1
qvalue=0.95
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
EPS = 0.622
cwrfnames={"TPW":"QVAPOR",
           "TCL":"QCLOUD", 
           "TCR":"QRAIN", 
           "TCI":"QICE", 
           "TCS":"QSNOW", 
           "TCG":"QGRAUP"
           }
seasonnames=["DJF","JJA","MAM","SON"]
"""
        fill_nocloud (:obj:`bool`, optional): Set to True to use fill values in 
            regions where clouds are not detected (optical depth less than 1). 
            Otherwise, the output will contain the surface temperature for 
            areas without clouds. Default is False.
            
        missing (:obj:`float`, optional): The fill value to use for areas 
            where no clouds are detected. Only used if *fill_nocloud* is 
            True. Default is 
            :data:`wrf.default_fill(numpy.float64)`. 
            
        opt_thresh (:obj:`float`, optional): The amount of optical 
            depth (integrated from top down) required to trigger a cloud top 
            temperature calculation. The cloud top temperature is calculated at 
            the vertical level where this threshold is met. Vertical columns 
            with less than this threshold will be treated as cloud free areas. 
            In general, the larger the value is for this 
            threshold, the lower the altitude will be for the cloud top 
            temperature calculation, and therefore higher cloud top 
            temperature values. Default is 1.0, which should be sufficient for 
            most users.
""" 
