# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 20:16:05 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from trajfuncs import BilinInterp, Add360, NearestIndex

sheddir = '/home/np838619/Watershed/shed_defs/'
lsmdir = '/home/np838619/Watershed/'

ncfile = Dataset(lsmdir+'ggis198901010000.nc','r')
lsm = ncfile.variables['LSM'][:]
eralon = ncfile.variables['longitude'][:]
eralat = ncfile.variables['latitude'][:]
ncfile.close()
lsm = pl.squeeze(lsm)

sheds = ['Am','AfMe','EAA','ArA','ArI','ArP','SOA','SOI','SOP']


lsm_sheds = []

for i in range(len(sheds)):
    rlspts = pl.genfromtxt(sheddir+sheds[i]+'_traj_release_new.txt',skip_header=5)
    rlspts = Add360(rlspts)
    msk_intrp = pl.zeros([rlspts.shape[0]])
    for j in range(rlspts.shape[0]):
        msk_intrp[j] = BilinInterp(rlspts[j],eralon,eralat,lsm)
    lsm_sheds.append(msk_intrp)


fig, ax = pl.subplots(3,3)

for i in range(9):
    axx = pl.subplot(3,3,i+1)
    axx.plot(lsm_sheds[i])
    axx.axhline(y=0.5,ls='--',color='k')
    pl.ylim(0,1)

pl.tight_layout()