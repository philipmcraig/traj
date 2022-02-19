# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 14:38:58 2017

@author: np838619
"""

import pylab as pl
from netCDF4 import Dataset
from kernelfuncs import *

janfiles = []; febfiles = []; marfiles = []; aprfiles = []; mayfiles = []
junfiles = []; julfiles = []; augfiles = []; sepfiles = []; octfiles = []
novfiles = []; decfiles = []

# for 2014 1st Jan 00z 2015 is included and should go into December list
for i in range(len(trajfiles)):
    if trajfiles[i][16:18] == '01':
        janfiles.append(trajfiles[i])
    elif trajfiles[i][16:18] == '02':
        febfiles.append(trajfiles[i])
    elif trajfiles[i][16:18] == '03':
        marfiles.append(trajfiles[i])
    elif trajfiles[i][16:18] == '04':
        aprfiles.append(trajfiles[i])
    elif trajfiles[i][16:18] == '05':
        mayfiles.append(trajfiles[i])
    elif trajfiles[i][16:18] == '06':
        junfiles.append(trajfiles[i])
    elif trajfiles[i][16:18] == '07':
        julfiles.append(trajfiles[i])
    elif trajfiles[i][16:18] == '08':
        augfiles.append(trajfiles[i])
    elif trajfiles[i][16:18] == '09':
        sepfiles.append(trajfiles[i])
    elif trajfiles[i][16:18] == '10':
        octfiles.append(trajfiles[i])
    elif trajfiles[i][16:18] == '11':
        novfiles.append(trajfiles[i])
    elif trajfiles[i][16:18] == '12':
        decfiles.append(trajfiles[i])

if trajfiles[-2][12:16] == '2014' and trajfiles[-1][12:16] == '2015':
    decfiles.append(trajfiles[-1])

# need to use the lengths of each list to split up catchcount lists
janlen = len(janfiles); feblen = len(febfiles); marlen = len(marfiles)
aprlen = len(aprfiles); maylen = len(mayfiles); junlen = len(junfiles)
jullen = len(julfiles); auglen = len(augfiles); seplen = len(sepfiles)
octlen = len(octfiles); novlen = len(novfiles); declen = len(decfiles)

# using the lengths for each month the netcdf files can be written here
# need a function specifying which month to use and size of time dimension
# also needs eralon & eralat as inputs
#MakeNCfile(janfiles,eralon,eralat,shed)
#MakeNCfile(febfiles,eralon,eralat,shed)
#MakeNCfile(marfiles,eralon,eralat,shed)
#MakeNCfile(aprfiles,eralon,eralat,shed)
#MakeNCfile(mayfiles,eralon,eralat,shed)
#MakeNCfile(junfiles,eralon,eralat,shed)
#MakeNCfile(julfiles,eralon,eralat,shed)
#MakeNCfile(augfiles,eralon,eralat,shed)
#MakeNCfile(sepfiles,eralon,eralat,shed)
#MakeNCfile(octfiles,eralon,eralat,shed)
#MakeNCfile(novfiles,eralon,eralat,shed)
#MakeNCfile(decfiles,eralon,eralat,shed)

# split the newlabs & catchcount arrays/lists into months:
janlabs = newlabs[:janlen,:,:]; sf_jan = sf3[:janlen,:,:] # January
jan_atl = atlcount[:janlen]; jan_ind = indcount[:janlen]
jan_pac = paccount[:janlen]; jan_arc = arccount[:janlen]
jan_sou = soucount[:janlen]
oldlen = janlen; newlen = oldlen + feblen
feblabs = newlabs[oldlen:newlen,:,:]; sf_feb = sf3[oldlen:newlen,:,:] # February
feb_atl = atlcount[oldlen:newlen]; feb_ind = indcount[oldlen:newlen]
feb_pac = paccount[oldlen:newlen]; feb_arc = arccount[oldlen:newlen]
feb_sou = soucount[oldlen:newlen]
oldlen = newlen; newlen = oldlen + marlen
marlabs = newlabs[oldlen:newlen,:,:]; sf_mar = sf3[oldlen:newlen,:,:] # March
mar_atl = atlcount[oldlen:newlen]; mar_ind = indcount[oldlen:newlen]
mar_pac = paccount[oldlen:newlen]; mar_arc = arccount[oldlen:newlen]
mar_sou = soucount[oldlen:newlen]
oldlen = newlen; newlen = oldlen + aprlen
aprlabs = newlabs[oldlen:newlen,:,:]; sf_apr = sf3[oldlen:newlen,:,:] # April
apr_atl = atlcount[oldlen:newlen]; apr_ind = indcount[oldlen:newlen]
apr_pac = paccount[oldlen:newlen]; apr_arc = arccount[oldlen:newlen]
apr_sou = soucount[oldlen:newlen]
oldlen = newlen; newlen = oldlen + maylen
maylabs = newlabs[oldlen:newlen,:,:]; sf_may = sf3[oldlen:newlen,:,:] # May
may_atl = atlcount[oldlen:newlen]; may_ind = indcount[oldlen:newlen]
may_pac = paccount[oldlen:newlen]; may_arc = arccount[oldlen:newlen]
may_sou = soucount[oldlen:newlen]
oldlen = newlen; newlen = oldlen + junlen
junlabs = newlabs[oldlen:newlen,:,:]; sf_jun = sf3[oldlen:newlen,:,:] # June
jun_atl = atlcount[oldlen:newlen]; jun_ind = indcount[oldlen:newlen]
jun_pac = paccount[oldlen:newlen]; jun_arc = arccount[oldlen:newlen]
jun_sou = soucount[oldlen:newlen]
oldlen = newlen; newlen = oldlen + jullen
jullabs = newlabs[oldlen:newlen,:,:]; sf_jul = sf3[oldlen:newlen,:,:] # July
jul_atl = atlcount[oldlen:newlen]; jul_ind = indcount[oldlen:newlen]
jul_pac = paccount[oldlen:newlen]; jul_arc = arccount[oldlen:newlen]
jul_sou = soucount[oldlen:newlen]
oldlen = newlen; newlen = oldlen + auglen
auglabs = newlabs[oldlen:newlen,:,:]; sf_aug = sf3[oldlen:newlen,:,:] # August
aug_atl = atlcount[oldlen:newlen]; aug_ind = indcount[oldlen:newlen]
aug_pac = paccount[oldlen:newlen]; aug_arc = arccount[oldlen:newlen]
aug_sou = soucount[oldlen:newlen]
oldlen = newlen; newlen = oldlen + seplen
seplabs = newlabs[oldlen:newlen,:,:]; sf_sep = sf3[oldlen:newlen,:,:] # September
sep_atl = atlcount[oldlen:newlen]; sep_ind = indcount[oldlen:newlen]
sep_pac = paccount[oldlen:newlen]; sep_arc = arccount[oldlen:newlen]
sep_sou = soucount[oldlen:newlen]
oldlen = newlen; newlen = oldlen + octlen
octlabs = newlabs[oldlen:newlen,:,:]; sf_oct = sf3[oldlen:newlen,:,:] # October
oct_atl = atlcount[oldlen:newlen]; oct_ind = indcount[oldlen:newlen]
oct_pac = paccount[oldlen:newlen]; oct_arc = arccount[oldlen:newlen]
oct_sou = soucount[oldlen:newlen]
oldlen = newlen; newlen = oldlen + novlen
novlabs = newlabs[oldlen:newlen,:,:]; sf_nov = sf3[oldlen:newlen,:,:] # November
nov_atl = atlcount[oldlen:newlen]; nov_ind = indcount[oldlen:newlen]
nov_pac = paccount[oldlen:newlen]; nov_arc = arccount[oldlen:newlen]
nov_sou = soucount[oldlen:newlen]
oldlen = newlen; newlen = oldlen + declen
declabs = newlabs[oldlen:newlen,:,:]; sf_dec = sf3[oldlen:newlen,:,:] # December
dec_atl = atlcount[oldlen:newlen]; dec_ind = indcount[oldlen:newlen]
dec_pac = paccount[oldlen:newlen]; dec_arc = arccount[oldlen:newlen]
dec_sou = soucount[oldlen:newlen]

# loop over each month in sequence and stick the data into the netcdf files
# need 2 kernel functions: 1 for full field (newlabs) & partioned fields (catchcount)
# also need option for all trajectories, CATI & CATII


# save origin density & flux density fields to netcdf:
# probably best to do this by month to avoid confusion

yr = trajfiles[0][12:16]
# January:
mnt = janfiles[0][16:18]
# all catchments:
all_DA01, all_FA01 = MainKernelFunc2(janlabs,eralon,eralat,sf_jan,rlslabs,'all',shed,'all',yr,mnt)
all_D101, all_D101 = MainKernelFunc2(janlabs,eralon,eralat,sf_jan,rlslabs,'1',shed,'all',yr,mnt)
all_D201, all_D201 = MainKernelFunc2(janlabs,eralon,eralat,sf_jan,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA01, atl_FA01 = MainKernelFunc(jan_atl,eralon,eralat,sf_jan,rlslabs,'all',shed,'atl',yr,mnt)
atl_D101, atl_F101 = MainKernelFunc(jan_atl,eralon,eralat,sf_jan,rlslabs,'1',shed,'atl',yr,mnt)
atl_D201, atl_F201 = MainKernelFunc(jan_atl,eralon,eralat,sf_jan,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA01, ind_FA01 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'all',shed,'ind',yr,mnt)
ind_D101, ind_F101 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'1',shed,'ind',yr,mnt)
ind_D201, ind_F201 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA01, pac_FA01 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'all',shed,'pac',yr,mnt)
pac_D101, pac_F101 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'1',shed,'pac',yr,mnt)
pac_D201, pac_F201 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA01, arc_FA01 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'all',shed,'arc',yr,mnt)
arc_D101, arc_F101 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'1',shed,'arc',yr,mnt)
arc_D201, arc_F201 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA01, sou_FA01 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'all',shed,'sou',yr,mnt)
sou_D101, sou_F101 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'1',shed,'sou',yr,mnt)
sou_D201, sou_F201 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'2',shed,'sou',yr,mnt)
print '**********January done**********'

# February:
mnt = febfiles[0][16:18]
# all catchments:
all_DA02, all_FA02 = MainKernelFunc2(feblabs,eralon,eralat,sf_feb,rlslabs,'all',shed,'all',yr,mnt)
all_D102, all_D102 = MainKernelFunc2(feblabs,eralon,eralat,sf_feb,rlslabs,'1',shed,'all',yr,mnt)
all_D202, all_D202 = MainKernelFunc2(feblabs,eralon,eralat,sf_feb,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA02, atl_FA02 = MainKernelFunc(feb_atl,eralon,eralat,sf_feb,rlslabs,'all',shed,'atl',yr,mnt)
atl_D102, atl_F102 = MainKernelFunc(feb_atl,eralon,eralat,sf_feb,rlslabs,'1',shed,'atl',yr,mnt)
atl_D202, atl_F202 = MainKernelFunc(feb_atl,eralon,eralat,sf_feb,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA02, ind_FA02 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'all',shed,'ind',yr,mnt)
ind_D102, ind_F102 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'1',shed,'ind',yr,mnt)
ind_D202, ind_F202 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA02, pac_FA02 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'all',shed,'pac',yr,mnt)
pac_D102, pac_F102 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'1',shed,'pac',yr,mnt)
pac_D202, pac_F202 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA02, arc_FA02 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'all',shed,'arc',yr,mnt)
arc_D102, arc_F102 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'1',shed,'arc',yr,mnt)
arc_D202, arc_F202 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA02, sou_FA02 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'all',shed,'sou',yr,mnt)
sou_D102, sou_F102 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'1',shed,'sou',yr,mnt)
sou_D202, sou_F202 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'2',shed,'sou',yr,mnt)
print '**********February done**********'

# March:
mnt = marfiles[0][16:18]
# all catchments:
all_DA03, all_FA03 = MainKernelFunc2(marlabs,eralon,eralat,sf_mar,rlslabs,'all',shed,'all',yr,mnt)
all_D103, all_D103 = MainKernelFunc2(marlabs,eralon,eralat,sf_mar,rlslabs,'1',shed,'all',yr,mnt)
all_D203, all_D203 = MainKernelFunc2(marlabs,eralon,eralat,sf_mar,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA03, atl_FA03 = MainKernelFunc(mar_atl,eralon,eralat,sf_mar,rlslabs,'all',shed,'atl',yr,mnt)
atl_D103, atl_F103 = MainKernelFunc(mar_atl,eralon,eralat,sf_mar,rlslabs,'1',shed,'atl',yr,mnt)
atl_D203, atl_F203 = MainKernelFunc(mar_atl,eralon,eralat,sf_mar,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA03, ind_FA03 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'all',shed,'ind',yr,mnt)
ind_D103, ind_F103 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'1',shed,'ind',yr,mnt)
ind_D203, ind_F203 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA03, pac_FA03 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'all',shed,'pac',yr,mnt)
pac_D103, pac_F103 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'1',shed,'pac',yr,mnt)
pac_D203, pac_F203 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA03, arc_FA03 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'all',shed,'arc',yr,mnt)
arc_D103, arc_F103 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'1',shed,'arc',yr,mnt)
arc_D203, arc_F203 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA03, sou_FA03 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'all',shed,'sou',yr,mnt)
sou_D103, sou_F103 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'1',shed,'sou',yr,mnt)
sou_D203, sou_F203 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'2',shed,'sou',yr,mnt)
print '**********March done**********'

# April:
mnt = aprfiles[0][16:18]
# all catchments:
all_DA04, all_FA04 = MainKernelFunc2(aprlabs,eralon,eralat,sf_apr,rlslabs,'all',shed,'all',yr,mnt)
all_D104, all_D104 = MainKernelFunc2(aprlabs,eralon,eralat,sf_apr,rlslabs,'1',shed,'all',yr,mnt)
all_D204, all_D204 = MainKernelFunc2(aprlabs,eralon,eralat,sf_apr,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA04, atl_FA04 = MainKernelFunc(apr_atl,eralon,eralat,sf_apr,rlslabs,'all',shed,'atl',yr,mnt)
atl_D104, atl_F104 = MainKernelFunc(apr_atl,eralon,eralat,sf_apr,rlslabs,'1',shed,'atl',yr,mnt)
atl_D204, atl_F204 = MainKernelFunc(apr_atl,eralon,eralat,sf_apr,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA04, ind_FA04 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'all',shed,'ind',yr,mnt)
ind_D104, ind_F104 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'1',shed,'ind',yr,mnt)
ind_D204, ind_F204 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA04, pac_FA04 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'all',shed,'pac',yr,mnt)
pac_D104, pac_F104 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'1',shed,'pac',yr,mnt)
pac_D204, pac_F204 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA04, arc_FA04 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'all',shed,'arc',yr,mnt)
arc_D104, arc_F104 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'1',shed,'arc',yr,mnt)
arc_D204, arc_F204 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA04, sou_FA04 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'all',shed,'sou',yr,mnt)
sou_D104, sou_F104 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'1',shed,'sou',yr,mnt)
sou_D204, sou_F204 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'2',shed,'sou',yr,mnt)
print '**********April done**********'

# May:
mnt = mayfiles[0][16:18]
# all catchments:
all_DA05, all_FA05 = MainKernelFunc2(maylabs,eralon,eralat,sf_may,rlslabs,'all',shed,'all',yr,mnt)
all_D105, all_D105 = MainKernelFunc2(maylabs,eralon,eralat,sf_may,rlslabs,'1',shed,'all',yr,mnt)
all_D205, all_D205 = MainKernelFunc2(maylabs,eralon,eralat,sf_may,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA05, atl_FA05 = MainKernelFunc(may_atl,eralon,eralat,sf_may,rlslabs,'all',shed,'atl',yr,mnt)
atl_D105, atl_F105 = MainKernelFunc(may_atl,eralon,eralat,sf_may,rlslabs,'1',shed,'atl',yr,mnt)
atl_D205, atl_F205 = MainKernelFunc(may_atl,eralon,eralat,sf_may,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA05, ind_FA05 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'all',shed,'ind',yr,mnt)
ind_D105, ind_F105 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'1',shed,'ind',yr,mnt)
ind_D205, ind_F205 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA05, pac_FA05 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'all',shed,'pac',yr,mnt)
pac_D105, pac_F105 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'1',shed,'pac',yr,mnt)
pac_D205, pac_F205 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA05, arc_FA05 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'all',shed,'arc',yr,mnt)
arc_D105, arc_F105 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'1',shed,'arc',yr,mnt)
arc_D205, arc_F205 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA05, sou_FA05 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'all',shed,'sou',yr,mnt)
sou_D105, sou_F105 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'1',shed,'sou',yr,mnt)
sou_D205, sou_F205 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'2',shed,'sou',yr,mnt)
print '**********May done**********'

# June:
mnt = junfiles[0][16:18]
# all catchments:
all_DA06, all_FA06 = MainKernelFunc2(junlabs,eralon,eralat,sf_jun,rlslabs,'all',shed,'all',yr,mnt)
all_D106, all_D106 = MainKernelFunc2(junlabs,eralon,eralat,sf_jun,rlslabs,'1',shed,'all',yr,mnt)
all_D206, all_D206 = MainKernelFunc2(junlabs,eralon,eralat,sf_jun,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA06, atl_FA06 = MainKernelFunc(jun_atl,eralon,eralat,sf_jun,rlslabs,'all',shed,'atl',yr,mnt)
atl_D106, atl_F106 = MainKernelFunc(jun_atl,eralon,eralat,sf_jun,rlslabs,'1',shed,'atl',yr,mnt)
atl_D206, atl_F206 = MainKernelFunc(jun_atl,eralon,eralat,sf_jun,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA06, ind_FA06 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'all',shed,'ind',yr,mnt)
ind_D106, ind_F106 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'1',shed,'ind',yr,mnt)
ind_D206, ind_F206 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA06, pac_FA06 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'all',shed,'pac',yr,mnt)
pac_D106, pac_F106 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'1',shed,'pac',yr,mnt)
pac_D206, pac_F206 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA06, arc_FA06 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'all',shed,'arc',yr,mnt)
arc_D106, arc_F106 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'1',shed,'arc',yr,mnt)
arc_D206, arc_F206 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA06, sou_FA06 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'all',shed,'sou',yr,mnt)
sou_D106, sou_F106 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'1',shed,'sou',yr,mnt)
sou_D206, sou_F206 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'2',shed,'sou',yr,mnt)
print '**********June done**********'

# July:
mnt = julfiles[0][16:18]
# all catchments:
all_DA07, all_FA07 = MainKernelFunc2(jullabs,eralon,eralat,sf_jul,rlslabs,'all',shed,'all',yr,mnt)
all_D107, all_D107 = MainKernelFunc2(jullabs,eralon,eralat,sf_jul,rlslabs,'1',shed,'all',yr,mnt)
all_D207, all_D207 = MainKernelFunc2(jullabs,eralon,eralat,sf_jul,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA07, atl_FA07 = MainKernelFunc(jul_atl,eralon,eralat,sf_jul,rlslabs,'all',shed,'atl',yr,mnt)
atl_D107, atl_F107 = MainKernelFunc(jul_atl,eralon,eralat,sf_jul,rlslabs,'1',shed,'atl',yr,mnt)
atl_D207, atl_F207 = MainKernelFunc(jul_atl,eralon,eralat,sf_jul,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA07, ind_FA07 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'all',shed,'ind',yr,mnt)
ind_D107, ind_F107 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'1',shed,'ind',yr,mnt)
ind_D207, ind_F207 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA07, pac_FA07 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'all',shed,'pac',yr,mnt)
pac_D107, pac_F107 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'1',shed,'pac',yr,mnt)
pac_D207, pac_F207 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA07, arc_FA07 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'all',shed,'arc',yr,mnt)
arc_D107, arc_F107 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'1',shed,'arc',yr,mnt)
arc_D207, arc_F207 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA07, sou_FA07 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'all',shed,'sou',yr,mnt)
sou_D107, sou_F107 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'1',shed,'sou',yr,mnt)
sou_D207, sou_F207 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'2',shed,'sou',yr,mnt)
print '**********July done**********'

# August:
mnt = augfiles[0][16:18]
# all catchments:
all_DA08, all_FA08 = MainKernelFunc2(auglabs,eralon,eralat,sf_aug,rlslabs,'all',shed,'all',yr,mnt)
all_D108, all_D108 = MainKernelFunc2(auglabs,eralon,eralat,sf_aug,rlslabs,'1',shed,'all',yr,mnt)
all_D208, all_D208 = MainKernelFunc2(auglabs,eralon,eralat,sf_aug,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA08, atl_FA08 = MainKernelFunc(aug_atl,eralon,eralat,sf_aug,rlslabs,'all',shed,'atl',yr,mnt)
atl_D108, atl_F108 = MainKernelFunc(aug_atl,eralon,eralat,sf_aug,rlslabs,'1',shed,'atl',yr,mnt)
atl_D208, atl_F208 = MainKernelFunc(aug_atl,eralon,eralat,sf_aug,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA08, ind_FA08 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'all',shed,'ind',yr,mnt)
ind_D108, ind_F108 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'1',shed,'ind',yr,mnt)
ind_D208, ind_F208 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA08, pac_FA08 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'all',shed,'pac',yr,mnt)
pac_D108, pac_F108 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'1',shed,'pac',yr,mnt)
pac_D208, pac_F208 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA08, arc_FA08 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'all',shed,'arc',yr,mnt)
arc_D108, arc_F108 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'1',shed,'arc',yr,mnt)
arc_D208, arc_F208 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA08, sou_FA08 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'all',shed,'sou',yr,mnt)
sou_D108, sou_F108 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'1',shed,'sou',yr,mnt)
sou_D208, sou_F208 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'2',shed,'sou',yr,mnt)
print '**********August done**********'

# September:
mnt = sepfiles[0][16:18]
# all catchments:
all_DA09, all_FA09 = MainKernelFunc2(seplabs,eralon,eralat,sf_sep,rlslabs,'all',shed,'all',yr,mnt)
all_D109, all_D109 = MainKernelFunc2(seplabs,eralon,eralat,sf_sep,rlslabs,'1',shed,'all',yr,mnt)
all_D209, all_D209 = MainKernelFunc2(seplabs,eralon,eralat,sf_sep,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA09, atl_FA09 = MainKernelFunc(sep_atl,eralon,eralat,sf_sep,rlslabs,'all',shed,'atl',yr,mnt)
atl_D109, atl_F109 = MainKernelFunc(sep_atl,eralon,eralat,sf_sep,rlslabs,'1',shed,'atl',yr,mnt)
atl_D209, atl_F209 = MainKernelFunc(sep_atl,eralon,eralat,sf_sep,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA09, ind_FA09 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'all',shed,'ind',yr,mnt)
ind_D109, ind_F109 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'1',shed,'ind',yr,mnt)
ind_D209, ind_F209 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA09, pac_FA09 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'all',shed,'pac',yr,mnt)
pac_D109, pac_F109 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'1',shed,'pac',yr,mnt)
pac_D209, pac_F209 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA09, arc_FA09 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'all',shed,'arc',yr,mnt)
arc_D109, arc_F109 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'1',shed,'arc',yr,mnt)
arc_D209, arc_F209 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA09, sou_FA09 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'all',shed,'sou',yr,mnt)
sou_D109, sou_F109 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'1',shed,'sou',yr,mnt)
sou_D209, sou_F209 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'2',shed,'sou',yr,mnt)
print '**********September done**********'

# October:
mnt = octfiles[0][16:18]
# all catchments:
all_DA10, all_FA10 = MainKernelFunc2(octlabs,eralon,eralat,sf_oct,rlslabs,'all',shed,'all',yr,mnt)
all_D110, all_D110 = MainKernelFunc2(octlabs,eralon,eralat,sf_oct,rlslabs,'1',shed,'all',yr,mnt)
all_D210, all_D210 = MainKernelFunc2(octlabs,eralon,eralat,sf_oct,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA10, atl_FA10 = MainKernelFunc(oct_atl,eralon,eralat,sf_oct,rlslabs,'all',shed,'atl',yr,mnt)
atl_D110, atl_F110 = MainKernelFunc(oct_atl,eralon,eralat,sf_oct,rlslabs,'1',shed,'atl',yr,mnt)
atl_D210, atl_F210 = MainKernelFunc(oct_atl,eralon,eralat,sf_oct,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA10, ind_FA10 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'all',shed,'ind',yr,mnt)
ind_D110, ind_F110 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'1',shed,'ind',yr,mnt)
ind_D210, ind_F210 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA10, pac_FA10 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'all',shed,'pac',yr,mnt)
pac_D110, pac_F110 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'1',shed,'pac',yr,mnt)
pac_D210, pac_F210 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA10, arc_FA10 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'all',shed,'arc',yr,mnt)
arc_D110, arc_F110 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'1',shed,'arc',yr,mnt)
arc_D210, arc_F210 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA10, sou_FA10 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'all',shed,'sou',yr,mnt)
sou_D110, sou_F110 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'1',shed,'sou',yr,mnt)
sou_D210, sou_F210 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'2',shed,'sou',yr,mnt)
print '**********October done**********'

# November:
mnt = novfiles[0][16:18]
# all catchments:
all_DA11, all_FA11 = MainKernelFunc2(novlabs,eralon,eralat,sf_nov,rlslabs,'all',shed,'all',yr,mnt)
all_D111, all_D111 = MainKernelFunc2(novlabs,eralon,eralat,sf_nov,rlslabs,'1',shed,'all',yr,mnt)
all_D211, all_D211 = MainKernelFunc2(novlabs,eralon,eralat,sf_nov,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA11, atl_FA11 = MainKernelFunc(nov_atl,eralon,eralat,sf_nov,rlslabs,'all',shed,'atl',yr,mnt)
atl_D111, atl_F111 = MainKernelFunc(nov_atl,eralon,eralat,sf_nov,rlslabs,'1',shed,'atl',yr,mnt)
atl_D211, atl_F211 = MainKernelFunc(nov_atl,eralon,eralat,sf_nov,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA11, ind_FA11 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'all',shed,'ind',yr,mnt)
ind_D111, ind_F111 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'1',shed,'ind',yr,mnt)
ind_D211, ind_F211 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA11, pac_FA11 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'all',shed,'pac',yr,mnt)
pac_D111, pac_F111 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'1',shed,'pac',yr,mnt)
pac_D211, pac_F211 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA11, arc_FA11 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'all',shed,'arc',yr,mnt)
arc_D111, arc_F111 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'1',shed,'arc',yr,mnt)
arc_D211, arc_F211 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA11, sou_FA11 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'all',shed,'sou',yr,mnt)
sou_D111, sou_F111 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'1',shed,'sou',yr,mnt)
sou_D211, sou_F211 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'2',shed,'sou',yr,mnt)
print '**********November done**********'

# December:
mnt = decfiles[0][16:18]
# all catchments:
all_DA12, all_FA12 = MainKernelFunc2(declabs,eralon,eralat,sf_dec,rlslabs,'all',shed,'all',yr,mnt)
all_D112, all_D112 = MainKernelFunc2(declabs,eralon,eralat,sf_dec,rlslabs,'1',shed,'all',yr,mnt)
all_D212, all_D212 = MainKernelFunc2(declabs,eralon,eralat,sf_dec,rlslabs,'2',shed,'all',yr,mnt)
# Atlantic:
atl_DA12, atl_FA12 = MainKernelFunc(dec_atl,eralon,eralat,sf_dec,rlslabs,'all',shed,'atl',yr,mnt)
atl_D112, atl_F112 = MainKernelFunc(dec_atl,eralon,eralat,sf_dec,rlslabs,'1',shed,'atl',yr,mnt)
atl_D212, atl_F212 = MainKernelFunc(dec_atl,eralon,eralat,sf_dec,rlslabs,'2',shed,'atl',yr,mnt)
# Indian:
ind_DA12, ind_FA12 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'all',shed,'ind',yr,mnt)
ind_D112, ind_F112 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'1',shed,'ind',yr,mnt)
ind_D212, ind_F212 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'2',shed,'ind',yr,mnt)
# Pacific:
pac_DA12, pac_FA12 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'all',shed,'pac',yr,mnt)
pac_D112, pac_F112 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'1',shed,'pac',yr,mnt)
pac_D212, pac_F212 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'2',shed,'pac',yr,mnt)
# Arctic:
arc_DA12, arc_FA12 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'all',shed,'arc',yr,mnt)
arc_D112, arc_F112 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'1',shed,'arc',yr,mnt)
arc_D212, arc_F212 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'2',shed,'arc',yr,mnt)
# Southern:
sou_DA12, sou_FA12 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'all',shed,'sou',yr,mnt)
sou_D112, sou_F112 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'1',shed,'sou',yr,mnt)
sou_D212, sou_F212 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'2',shed,'sou',yr,mnt)
print '**********December done**********'

#SaveFields(junfiles,all_DA06,arp_all_FA06,shed,'all','all')
#SaveFields(junfiles,atl_DA06,arp_atl_FA06,shed,'atl','all')

