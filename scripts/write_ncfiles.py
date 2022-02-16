# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 16:20:12 2017

@author: np838619
"""

import pylab as pl
from netCDF4 import Dataset
from kernelfuncs import MakeNCfile

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

year = '2014'
datadir = '/net/glusterfs/scenario/users/np838619/traj/'+year+'/amr/TD/'

trajfiles = PrintFiles(datadir,'utraj-df')
trajfiles = pl.sort(trajfiles)

janfiles = []; febfiles = []; marfiles = []; aprfiles = []; mayfiles = []
junfiles = []; julfiles = []; augfiles = []; sepfiles = []; octfiles = []
novfiles = []; decfiles = []

# for 2014 1st Jan 00z 2015 is included and should go into December list
if year == '2014':
    NRLS = len(trajfiles)-1
else:
    NRLS = len(trajfiles)

for i in range(NRLS):
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

if year == '2014':
    decfiles.append(trajfiles[-1])

nc1 = Dataset('/home/np838619/Watershed/ggis198901010000.nc','r')
eralon = nc1.variables['longitude'][:]; eralat = nc1.variables['latitude'][:]
nc1.close()

sheds = ['amr', 'afr', 'eaa', 'ara', 'ari', 'arp', 'soa', 'soi', 'sop']

for i in range(len(sheds)):
    MakeNCfile(janfiles,eralon,eralat,sheds[i])
    MakeNCfile(febfiles,eralon,eralat,sheds[i])
    MakeNCfile(marfiles,eralon,eralat,sheds[i])
    MakeNCfile(aprfiles,eralon,eralat,sheds[i])
    MakeNCfile(mayfiles,eralon,eralat,sheds[i])
    MakeNCfile(junfiles,eralon,eralat,sheds[i])
    MakeNCfile(julfiles,eralon,eralat,sheds[i])
    MakeNCfile(augfiles,eralon,eralat,sheds[i])
    MakeNCfile(sepfiles,eralon,eralat,sheds[i])
    MakeNCfile(octfiles,eralon,eralat,sheds[i])
    MakeNCfile(novfiles,eralon,eralat,sheds[i])
    MakeNCfile(decfiles,eralon,eralat,sheds[i])