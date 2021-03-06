var_parameters = {
"QFX"              : { "compute_mode" : 1  , "shiftday" :  0 , "vert_intp" : None ,"fields":{"QFX":None}} , 
"LH"               : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"LH":None}} , 
"GSW"              : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"GSW":None}} , 
"ATH2"             : { "compute_mode" : 1  , "shiftday" :  -1 , "vert_intp" : None ,"fields":{"ATH2":None}} , 
"XLAI"             : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XLAI":None}} , 
"WSAVG"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"WSAVG":None}} , 
"WSMIN"            : { "compute_mode" : 9  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"WSMIN":None}} , 
"WSMAX"            : { "compute_mode" : 8  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"WSMAX":None}} , 
"XFCOV"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XFCOV":None}} , 
"XQRCHRG"          : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XQRCHRG":None}} , 
"XRNOF"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XRNOF":None}} , 
"XRSAT"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XRSAT":None}} , 
"XRDRN"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XRDRN":None}} , 
"XRBAS"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XRBAS":None}} , 
"XRSUR"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XRSUR":None}} , 
"XTAUX"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XTAUX":None}} , 
"XTAUY"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XTAUY":None}} , 
"XUPTKW"           : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XUPTKW":None}} , 
"XRTWS"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XRTWS":None}} , 
"XSWDEPTH"         : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XSWDEPTH":None}} , 
"AGFX"             : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"AGFX":None}} , 
"CANWAT"           : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"CANWAT":None}} , 
"ACSNOW"           : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ACSNOW":None}} , 
"AGFX"             : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"AGFX":None}} , 
"ARBAS"            : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ARBAS":None}} , 
"ATSK"             : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ATSK":None}} , 
"PRMIN"            : { "compute_mode" : 9  , "shiftday" : 1  , "vert_intp" : None ,"fields":{"PRMIN":None}} , 
"PRMAX"            : { "compute_mode" : 8  , "shiftday" : 1  , "vert_intp" : None ,"fields":{"PRMAX":None}} , 
"PRAVG"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"PRAVG":None}} , 
"BIOLAI"           : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"BIOLAI":None}} , 
"STEM"             : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"dim":3,
                           "fields":{"STEM":{"units":"Mg/ha","description":"stem biomas"}, }},
"STOMATAWS"        : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"dim":3,
                           "fields":{"STOMATAWS":{"units":"-","description":"stomatal water stress"}, }},
"LEAF"             : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"dim":3,
                           "fields":{"LEAF":{"units":"Mg/ha","description":"leaf biomas"}, }},
"AT2M"             : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"dim":3,
                           "fields":{"AT2M":{"units":"K","description":"temperatrue at 2 meter "},
                                   #"T2M975":{"units":"K","description":"temperatrue at 2 meter 975 percentile "}
                             }},
"Q2M"             : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"dim":3,
                           "fields":{"Q2M":{"units":"kg kg-1","description":"specific humidity at 2 meter "},
                             }},
"RAINC"            : { "compute_mode" : 6  , "shiftday" : 0   , "vert_intp" : None ,"dim":3,
                           "fields":{"RAINC":{"units":"mm/day","description":"Average convective daily Precip"} }} , 
"PR"            : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"dim":3,
                           "fields":{"PRAVG":{"units":"mm/day","description":"Average daily Precip"},
                                      "PMAX":{"units":"mm/day","description":"Maximum"},
                                       "PCT":{"units":"mm/day","description":"Precip 95 percentile"},
                                 "RAINYDAYS":{"units":"day","description":"No. of days with PR larger than 1mm/day"},
                                       "R10":{"units":"mm/day","description":"No. of days with precipitation >10 mm d-1"},
                                       "R5D":{"units":"mm/day","description":"Maximum 5 d precipitation total"},
                                       "CDD":{"units":"day","description":"Consective dry day"},
                                      "R95T":{"units":"unitless","description":"Fraction of annual total precipitation due to events exceeding the 1961-1990 95th percentile"},
                                      "SDII":{"units":"mm/day","description":"Simple daily intensity index"}}} , 
"AOUTFLOW"         : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"AOUTFLOW":None}} , 
"XSABVG"           : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XSABVG":None}} , 
"ARNOF"            : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ARNOF":None}} , 
"ARSUR"            : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ARSUR":None}} , 
"SNOWNC"           : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"SNOWNC":None}} , 
"SNOWC"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"SNOWC":None}} , 
"EMISS"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"EMISS":None}} , 
"UQg"              : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"UQg":None}} , 
"VQg"              : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"VQg":None}} , 
#"AT2M"             : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"AT2M":None}} , 
"AODVIS"             : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"AODVIS":None}} , 
"AODNIR"             : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"AODNIR":None}} , 
"T2M"              : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"T2M":None}} , 
"SST"              : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"dim":3,
                       "fields":{"SST":{"units":"K","description":"Sea surface temperature "},
                                }
                     },
"TSK"              : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"dim":3,
                       "fields":{"TSK":{"units":"K","description":"Skin temperature "},
                                }
                     },
"XTSS"             : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XTSS":None}} , 
"T2"               : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"T2":None}} , 
"AU2M"             : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"AU2M":None}} , 
"AV2M"             : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"AV2M":None}} , 
"AQ2M"             : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"AQ2M":None}} , 
"T2MAX"            : { "compute_mode" : 8  , "shiftday" : 0  , "vert_intp" : None ,"dim":3,
                       "fields":{"T2MAX":{"units":"K","description":"maxium temperatrue at 2 meter "},
                                }
                     },
"T2MIN"            : { "compute_mode" : 9  , "shiftday" : 0  , "vert_intp" : None ,"dim":3,
                       "fields":{"T2MIN":{"units":"K","description":"minium temperatrue at 2 meter "},
                                }
                     },
"WIN"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"dim":3,
                       "fields":{"WIN_10":{"units":"m/s","description":" Wind Speed at 10m"},
                                }
                     },
"xsmtg"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"dim":4,
                       "nlev":4,
                       "fields":{"xsmtg":{"units":"g","description":" Total grouped SM(liq+ice)"},
                                }
                     },
"xsmlg"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"dim":4,
                       "nlev":4,
                       "fields":{"xsmlg":{"units":"g","description":" liq grouped SM"},
                                }
                     },
"xsmig"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"dim":4,
                       "nlev":4,
                       "fields":{"xsmig":{"units":"g","description":" ice grouped SM"},
                                }
                     },
#"RAINC"            : { "compute_mode" : 6  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"RAINC":None}} , 
"RAINNC"           : { "compute_mode" : 6  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"RAINNC":None}} , 
"CDD"              : { "compute_mode" : 12 , "shiftday" : 0  , "vert_intp" : None ,"fields":{"CDD":None}} , 
"RAINYDAYS"        : { "compute_mode" : 13 , "shiftday" : 0  , "vert_intp" : None ,"fields":{"RAINYDAYS":None}} , 
"R10"              : { "compute_mode" : 13 , "shiftday" : 0  , "vert_intp" : None ,"fields":{"R10":None}} , 
"R5D"              : { "compute_mode" : 14 , "shiftday" : 0  , "vert_intp" : None ,"fields":{"R5D":None}} , 
"SDII"             : { "compute_mode" : 15 , "shiftday" : 0  , "vert_intp" : None ,"fields":{"SDII":None}} , 
"R95T"             : { "compute_mode" : 16 , "shiftday" : 0  , "vert_intp" : None ,"fields":{"R95T":None}} , 
"PSFC"             : { "compute_mode" : 1  , "shiftday" : 0 , "vert_intp" : None ,"fields":{"PSFC":None}} , 
"PMSL"             : { "compute_mode" : 1  , "shiftday" : 0 , "vert_intp" : None ,"fields":{"PMSL":None}} , 
"AU10"             : { "compute_mode" : 1  , "shiftday" : -1 , "vert_intp" : None ,"fields":{"AU10":None}} , 
"AV10"             : { "compute_mode" : 1  , "shiftday" : -1 , "vert_intp" : None ,"fields":{"AV10":None}} , 
"AQ2M"             : { "compute_mode" : 1  , "shiftday" : -1 , "vert_intp" : None ,"fields":{"AQ2M":None}} , 
"ASWUPT"           : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ASWUPT":None}} , 
"ASWUPTC"          : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ASWUPTC":None}} , 
"ASWDNT"           : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ASWDNT":None}} , 
#"ASWDNS"           : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ASWDNS":None}} , 
"ASWDNS"           : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None  ,"dim":3,
                       "fields"       :{"ASWDNS":{"units":"w m-2","description":"short wave downwelling at surface"} }} , 
"ASWDNSC"          : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ASWDNSC":None}} , 
"ASWUPS"           : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ASWUPS":None}} , 
"ASWUPSC"          : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ASWUPSC":None}} , 
"ALWUPT"           : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ALWUPT":None}} , 
"ALWUPTC"          : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ALWUPTC":None}} , 
"ALWDNS"           : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ALWDNS":None}} , 
"ALWDNSC"          : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ALWDNSC":None}} , 
"ALWUPS"           : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ALWUPS":None}} , 
"ALWUPSC"          : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ALWUPSC":None}} , 
"TCWPC"            : { "compute_mode" : 2  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields"       :{"CWPC":{"units":"g m-2","description":"Total cloud liquid water path"} }} , 
"TCWPI"            : { "compute_mode" : 2  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields"       :{"CWPI":{"units":"g m-2","description":"Total rain  water path"} }} , 
"TCWPR"            : { "compute_mode" : 2  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields"       :{"CWPR":{"units":"g m-2","description":"Total ice liquid water path"} }} , 
"CLDFRA"           : { "compute_mode" : 1  , "shiftday" :  0 , "vert_intp" : None ,"fields":{"CLDFRA":None}} , 
"XSWI"             : { "compute_mode" : 1  , "shiftday" :  0 , "vert_intp" : None ,"fields":{"XSWI":None}} , 
"PBLH"             : { "compute_mode" : 1  , "shiftday" :  0 , "vert_intp" : None ,"fields":{"PBLH":None}} , 
"TKES"             : { "compute_mode" : 1  , "shiftday" :  0 , "vert_intp" : None ,"fields":{"TKES":None}} , 
"UST"              : { "compute_mode" : 1  , "shiftday" :  0 , "vert_intp" : None ,"fields":{"UST":None}} , 
"WCSTAR"           : { "compute_mode" : 1  , "shiftday" :  0 , "vert_intp" : None ,"fields":{"WCSTAR":None}} , 
"BR"               : { "compute_mode" : 1  , "shiftday" :  0 , "vert_intp" : None ,"fields":{"BR":None}} , 
"POTEVP"           : { "compute_mode" : 5  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"POTEVP":None}} , 
"SFROFF"           : { "compute_mode" : 5  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"SFROFF":None}} , 
"XRNOF"           : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XRNOF":None}} , 
"UDROFF"           : { "compute_mode" : 5  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"UDROFF":None}} , 
"SNOW"             : { "compute_mode" : 1  , "shiftday" :  0 , "vert_intp" : None ,"fields":{"SNOW":None}} , 
"SNOWFALL"         : { "compute_mode" : 12 , "shiftday" :  0 , "vert_intp" : None ,"fields":{"SNOWFALL":None}} , 
"SNOWH"            : { "compute_mode" : 1  , "shiftday" :  0 , "vert_intp" : None ,"fields":{"SNOWH":None}} , 
"XTSS"             : { "compute_mode" : 1  , "shiftday" :  0 , "vert_intp" : None ,"fields":{"XTSS":None}} , 
"XOUTFLOW"         : { "compute_mode" : 1  , "shiftday" :  1 , "vert_intp" : None ,"fields":{"XOUTFLOW":None}} , 
"XPARSUN"          : { "compute_mode" : 1  , "shiftday" :  1 , "vert_intp" : None ,"fields":{"XPARSUN":None}} , 
"XPARSHA"          : { "compute_mode" : 1  , "shiftday" :  1 , "vert_intp" : None ,"fields":{"XPARSHA":None}} , 
"AHFX"             : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"AHFX":None}} , 
"ALFX"             : { "compute_mode" : 1  , "shiftday" : -1  , "vert_intp" : None ,"fields":{"ALFX":None}} , 
"XETR"             : { "compute_mode" : 1  , "shiftday" :  1 , "vert_intp" : None ,"fields":{"XETR":None}} , 
"XFEVPL"           : { "compute_mode" : 1  , "shiftday" :  1 , "vert_intp" : None ,"fields":{"XFEVPL":None}} , 
"XFEVPG"           : { "compute_mode" : 1  , "shiftday" :  1 , "vert_intp" : None ,"fields":{"XFEVPG":None}} , 
"XZWT"             : { "compute_mode" : 1  , "shiftday" :  1 , "vert_intp" : None ,"fields":{"XZWT":None}} , 
"PS"               : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"PS":None}} , 
"PSL"              : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"PSL":None}} , 
"H"                : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"H":None}} , 
"TPW"              : { "compute_mode" : 1  , "shiftday" :  1 , "vert_intp" : None ,"fields":{"TPW":None}} , 
"AUQV"             : { "compute_mode" : 1  , "shiftday" :  -1 , "vert_intp" : None ,"fields":{"AUQV":None}} , 
"AVQV"             : { "compute_mode" : 1  , "shiftday" :  -1 , "vert_intp" : None ,"fields":{"AVQV":None}} , 
"QVAPOR"           : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"QVAPOR":None}} , 
"QCLOUD"           : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"QCLOUD":None}} , 
"QRAIN"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"QRAIN":None}} , 
"QICE"             : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"QICE":None}} , 
"QSNOW"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"QSNOW":None}} , 
"QGRAUP"           : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"QGRAUP":None}} , 
"T"                : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"fields":{"T":None}} , 
"U"                : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"fields":{"U":None}} , 
"V"                : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"fields":{"V":None}} , 
"W"                : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"fields":{"W":None}} , 
"QA"               : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"fields":{"QA":None}} , 
"UA"               : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"fields":{"UA":None}} , 
"VA"               : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"fields":{"VA":None}} , 
"TA"               : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"fields":{"TA":None}} , 
"ALBEDO"           : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"ALBEDO":None}} , 
"XWLIQ"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XWLIQ":None}} , 
"XWICE"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XWICE":None}} , 
"XALB"             : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XALB":None}} , 
"XALBV"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XALBV":None}} , 
"XALBG"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XALBG":None}} , 
"XSSUN"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XSSUN":None}} , 
"XSSHA"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XSSHA":None}} , 
"XSABG"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XSABG":None}} , 
"XSABVG"           : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XSABVG":None}} , 
"XSABVSUN"         : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XSABVSUN":None}} , 
"XSABVSHA"         : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"XSABVSHA":None}} , 
"TSLB"             : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"TSLB":None}} , 
"SMOIS"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None ,"fields":{"SMOIS":None}} , 
"omega"            : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : 'p'  ,"dim":4,
                       "fields":{
                                "omega":{"units":"Pa s-1","description":"omega"} 
                                }
                     }, 
"RH"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{
                                "RH":{"units":"%","description":"Relative Humidity at 2m"} 
                              }} , 
"Q"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : 'p'  ,"dim":4,
                       "fields":{
                                "Q":{"units":"Kg/Kg","description":"Specific Humidity"} 
                              }} , 
"rh3d"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" :None  ,"dim":4,
                       "fields":{
                                "rh3d":{"units":"%","description":"Relative Humidity"} 
                              }} , 
"temp"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"dim":4,
                       "fields":{
                                "tk":{"units":"K","description":"Temperature"}, 
                             "theta":{"units":"K","description":"Potential Temperature"} }} , 
"conv"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{
                               "CAPE":{"units":"J/kg","description":"Convective Available Potential Energy "}, 
                                "CIN":{"units":"J/kg","description":"Convective Inhibition"} }} , 
"slp"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{ "slp":{"units":"Pa","description":"sea level pressure"} }} , 
"lwp"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{ "lwp":{"units":"Kg","description":"liquid water path"} }} , 
"iwp"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{ "iwp":{"units":"Kg","description":"ice water path"} }} , 
"cto"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{ "ctt":{"units":"K","description":"Cloud top temperature"}, 
                                  "cth":{"units":"meter","description":"Cloud top height(above topography)"} 
                         }} , 
####QU QV
"Qflux"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"dim":4,
                       "fields":{
                                "uq":{"units":"m s-1","description":"Zonal moisture flux "}, 
                                "vq":{"units":"m s-1","description":"Meridional moisture flux"} }} , 
####QU QV
"TPW"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{ "TPW_l":{"units":"cm","description":"totol precipitable water Surface to 680mb"} 
                                 ,"TPW_m":{"units":"cm","description":"totol precipitable water 680 to 440 mb"} 
                                 ,"TPW_h":{"units":"cm","description":"totol precipitable water 440 to model top"} 
                         }} , 
"TPWmax"   : { "compute_mode" : 8  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{ "TPW_l":{"units":"cm","description":"totol precipitable water Surface to 680mb"} 
                                 ,"TPW_m":{"units":"cm","description":"totol precipitable water 680 to 440 mb"} 
                                 ,"TPW_h":{"units":"cm","description":"totol precipitable water 440 to model top"} 
                         }} , 

#for cloud related
"TCL"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{ "TCL_l":{"units":"cm","description":"totol cloud liquid water Surface to 680mb"} 
                                 ,"TCL_m":{"units":"cm","description":"totol cloud liquid water 680 to 440 mb"} 
                                 ,"TCL_h":{"units":"cm","description":"totol cloud liquid water 440 to model top"} }} , 
"TCR"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{ "TCR_l":{"units":"cm","description":"totol cloud rain water Surface to 680mb"} 
                                 ,"TCR_m":{"units":"cm","description":"totol cloud rain water 680 to 440 mb"} 
                                 ,"TCR_h":{"units":"cm","description":"totol cloud rain water 440 to model top"} }} , 
"TCI"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{ "TCI_l":{"units":"cm","description":"totol cloud ice water Surface to 680mb"} 
                                 ,"TCI_m":{"units":"cm","description":"totol cloud ice water 680 to 440 mb"} 
                                 ,"TCI_h":{"units":"cm","description":"totol cloud ice water 440 to model top"} }} , 
"TCS"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{ "TCS_l":{"units":"cm","description":"totol cloud snow water Surface to 680mb"} 
                                 ,"TCS_m":{"units":"cm","description":"totol cloud snow water 680 to 440 mb"} 
                                 ,"TCS_h":{"units":"cm","description":"totol cloud snow water 440 to model top"} }} , 
"TCG"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{ "TCG_l":{"units":"cm","description":"totol cloud graupel water Surface to 680mb"} 
                                 ,"TCG_m":{"units":"cm","description":"totol cloud graupel water 680 to 440 mb"} 
                                 ,"TCG_h":{"units":"cm","description":"totol cloud graupel water 440 to model top"} }} , 
#for cloud related
"cldfrag"      : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{
                               "cldfra_total":{"units":"unitless","description":"total cloud"}, 
                               "cldfra_high" :{"units":"unitless","description":"high cloud"}, 
                               "cldfra_mid"  :{"units":"unitless","description":"mid cloud"}, 
                               "cldfra_low"  :{"units":"unitless","description":"low cloud"}   }} ,
"geopt"     : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"dim":4,
                       "fields":{"geopt":{"units":"m2 s-2","description":"Geopotential"} }} , 
"height"    : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"dim":4,
                       "fields":{"height":{"units":"m","description":"Height"} }} , 
"uv_met"    : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : "p"  ,"dim":4,
                       "fields":{
                              #"WIN":{"units":"m s-1","description":"Wind Speed"} ,
                              "u_met":{"units":"m s-1","description":"Rotated u-wind component"} ,
                              "v_met":{"units":"m s-1","description":"Rotated v-wind component"}}} ,
"uv_10"    : { "compute_mode" : 1  , "shiftday" : 0  , "vert_intp" : None  ,"dim":3,
                       "fields":{
                              #"WIN_10":{"units":"m s-1","description":"Wind Speed at 10m"} ,
                              "u_10":{"units":"m s-1","description":"Rotated u-wind component at 10m"} ,
                              "v_10":{"units":"m s-1","description":"Rotated v-wind component at 10m"}}} ,
} 
#save seasonal data of following vars:

seasonal_allvars  = (         'PRAVG'     , 'PRMAX'   , 'PRMIN'    , 'RAINC'   , 'RAINNC'  , 
                  'AU10'    , 'AV10'      , 'AQ2M'    , 'WSAVG'    , 'WSMAX'   , 
                  'AT2M'    , 'T2MAX'     , 'T2MIN'   , 
                  'ASWUPT'  , 'ASWUPTC'   , 'ASWDNT'  , 'ASWDNS'   , 'ASWDNSC' , 'ASWUPS'  , 'ASWUPSC' , 
                  'ALWUPT'  , 'ALWUPTC'   , 'ALWDNS'  , 'ALWDNSC'  , 'ALWUPS'  , 'ALWUPSC' , 
                  'TCWPC'   , 'TCWPI'     , 'TCWPR'   , 
                  'CLDFRAg' , 'CLDFRApbl' , 'CLDFRA'  , 'CLDCNVg'  , 'CLDSCUg' , 
                  'EMISS'   , 'QFX'       , 'LH'      , 
                  'XSWI'    , 'PBLH'      , 'TKES'    , 'UST'      , 'WCSTAR'  , 
                  'BR'      , 'POTEVP'    , 'SFROFF'  , 'UDROFF'   , 'SNOW'    , 
                  'SNOWH'   , 'XTSS'      , 'XSMTg'   , 'XOUTFLOW' , 'XPARSUN' , 'XPARSHA' , 
                  'AHFX'    , 'ALFX'      , 'XETR'    , 'XFEVPL'   , 'XFEVPG'  , 'XZWT'    , 
                  'ALBEDO'  , 'TSK'       , 
                  'TSLB'    , 'SMOIS'     , 'XQRCHRG' , 
                  'PS'      , 'PSL'       ,
                  #3D
                  'H'       , 'T'         , 'Q'       , 'U'        , 'V'       , 'W',
                  'QCLOUD'  , 'QRAIN'     , 'QICE'     
				 	)


#save monthly data of following vars:
monthly_allvars  = (           'PRAVG'     , 'PRMAX'   , 'PRMIN'    , 'RAINC'   , 'RAINNC'  , 
                  'AU10'    , 'AV10'      , 'AQ2M'    , 'WSAVG'    , 'WSMAX'   , 
                  'AT2M'    , 'T2MAX'     , 'T2MIN'   , 
                  'ASWUPT'  , 'ASWUPTC'   , 'ASWDNT'  , 'ASWDNS'   , 'ASWDNSC' , 'ASWUPS'  , 'ASWUPSC' , 
                  'ALWUPT'  , 'ALWUPTC'   , 'ALWDNS'  , 'ALWDNSC'  , 'ALWUPS'  , 'ALWUPSC' , 
                  'TCWPC'   , 'TCWPI'     , 'TCWPR'   , 
                  'CLDFRAg' , 'CLDFRApbl' , 'CLDFRA'  , 'CLDCNVg'  , 'CLDSCUg' , 
                  'EMISS'   , 'QFX'       , 'LH'      , 
                  'XSWI'    , 'PBLH'      , 'TKES'    , 'UST'      , 'WCSTAR'  , 
                  'BR'      , 'POTEVP'    , 'SFROFF'  , 'UDROFF'   , 'SNOW'    , 
                  'SNOWH'   , 'XTSS'      , 'XSMTg'   , 'XOUTFLOW' , 'XPARSUN' , 'XPARSHA' , 
                  'AHFX'    , 'ALFX'      , 'XETR'    , 'XFEVPL'   , 'XFEVPG'  , 'XZWT'    , 
                  'ALBEDO'  , 'TSK'       , 
                  'TSLB'    , 'SMOIS'     , 'XQRCHRG' , 
                  'PS'      , 'PSL'       ,
                  #3D
                  'H'       , 'T'         , 'Q'       , 'U'        , 'V'       , 'W',
                  'QCLOUD'  , 'QRAIN'     , 'QICE'     
				 	)

# save daily data of following vars 
daily_allvars_2d =[ 'slp', 'PR',     'RAINC',   'RAINNC',
                  'AU10'    , 'AV10'      , 'AQ2M'    , 'WSAVG'    , 'WSMAX'   , 
                    'AT2M',   'T2MAX',   'T2MIN',
                  'ASWUPT',   'ASWUPTC', 'ASWDNT', 'ASWDNS',  'ASWDNSC',  'ASWUPS',  'ASWUPSC',
                  'ALWUPT',   'ALWUPTC', 'ALWDNS', 'ALWDNSC', 'ALWUPS',   'ALWUPSC',
                    'EMISS',   'QFX',    'LH',
                  'XSWI'    , 'PBLH'      , 'TKES'    , 'UST'      , 'WCSTAR'  , 
                  'AHFX'    , 'ALFX'      , 
                  'ALBEDO'  , 'TSK'       , 
                  'XQRCHRG' , 
                  'PSFC'    , 'PMSL'      , 'cldfrag'    ,
                      'SNOW',   'SNOWH', 'uv_10', 'RH',
                   'XZWT',    'XRBAS',   'XRDRN',  'XRSAT', 'TPWmax','lwp','iwp',"TPW","cto"]

#daily_allvars_2d =[ 'RH', 'cldfrag','TPWmax','PBLH','AHFX']

daily_allvars_3d =[ 'uv_met',  'height',"xsmtg",
                   'temp'   , 'XWLIQ',   'XWICE',
                   'QVAPOR',  'QCLOUD'  , 'QRAIN'     , 'QICE'     ]

# diurnal cycle vars:
diurnal_allvars  = (       'PRAVG','RAINC'   ,'RAINNC'  ,
                       'AT2M'    , 'T2MAX'     , 'T2MIN'   , 
                       'AU10'    , 'AV10'      , 'AQ2M'    , 
                       'TCWPC'   , 'TCWPI'     , 
                       'CLDFRAg' , 'CLDFRApbl' , 'CLDCNVg' , 'CLDSCUg' , 
                       'XSWI'    , 'PBLH'      , 'TKES'    , 'UST'     , 'WCSTAR' , 
                       'BR'      , 'SNOW'      , 
                       'SNOWH'   , 'XTSS'      , 'XSMTg'   , 'XPARSUN' , 
                       'XPARSHA' , 'AHFX'      , 'ALFX'    , 'XETR'    , 'XFEVPL' , 
                       'XFEVPG'  , 'XZWT'      , 'ALBEDO'  , 'TSK'
				 	)


# hourly vars:
hourly_allvars  = ('PRAVG' , 'H'      , 'temp'     , 'U'    , 'V' , 'omega' ,'PSFC', 
                        'Q', 'QCLOUD' , 'QRAIN' , 'QICE'
				 	)
#postList={"seasonal":seasonal_allvars,"monthly":monthly_allvars,"daily":daily_allvars,"diurnal":diurnal_allvars} #,
      #  "3hourly":hourly3_allvars}
