# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 18:39:09 2019

@author: qx911590
"""

from __future__ import division
import pylab as pl
import xarray as xr
#from netCDF4 import Dataset
import glob

clusdir = '/storage/shared/glusterfs/scenario/users/qx911590/'

years = pl.linspace(2010,2014,5).astype(int)
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
sheds = ['amr','afr','eaa','ara','ari','arp','soa','soi','sop']
basins = ['all','atl','ind','pac','arc','sou']
ortypes = ['all','CAT I','CAT II']

# loop over years:
    # loop over sheds:
        # loop over files (can do this by month, or just in order):
            # open a netcdf file with xarray
            # open each array iteratively:
                # do an 'origin' and a 'flux' density array at same time
                # loop over basins:
                    # loop over ortypes:
                        #data1 - ncfile.variables['origin density '+sheds[i]+basins[j]+ ' +ortypes[k]]
                        #data2 = ncfile.variables['flux density '+sheds[i]+basins[j]+ ' +ortypes[k]]
                        #data1 = xr.DataArray(data1)
                        #data2 = xr.DataArray(data2)

                        # might be too much looping
                        # could just open all the arrays in ncfile together
                        # may consider using multiprocessing

filenames = glob.glob(clusdir+'traj/monthly_means/'+str(years[0])+'/'+sheds[0]+'/*')

ncfile = xr.open_dataset(filenames[0])

data1 = ncfile.variables['origin density amrall all'][:]
data1 = xr.DataArray(data1)

data2 = ncfile.variables['flux density amrall all'][:]
data2 = xr.DataArray(data2)