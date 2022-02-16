# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 10:53:14 2017

@author: np838619
"""

clusdir = '/net/glusterfs/scenario/users/np838619/traj/'
newnc = Dataset(clusdir+'test.nc',mode='a',format='NETCDF4')
#lat_dim = newnc.createDimension('lat',256)
#lon_dim = newnc.createDimension('lon',512)
#lat_in = newnc.createVariable('lat',pl.float64,('lat',))
#lat_in.units = 'degrees_north'
#lat_in.long_name = 'latitude'
#lon_in = newnc.createVariable('lon',pl.float64,('lon',))
#lon_in.units = 'degrees_east'
#lon_in.long_name = 'longitude'
#lat_in = eralat
#lon_in = eralon

#time_dim = newnc.createDimension('time',60)
#time = newnc.createVariable('time',pl.float64,('time',))
#time.units = 'hours'
#time.long_name = 'release date/time'

ordens = newnc.createVariable('origin density sou CAT2',pl.float64,('time','lat','lon'))
ordens.units = 'trajectories per steradian'
ordens.standard_name = 'density of trajectory origins'
ordens[:,:,:] = Kz[:60]

fxdens = newnc.createVariable('flux density sou CAT2',pl.float64,('time','lat','lon'))
fxdens.units = 'Sv per steradian'
fxdens.standard_name = 'flux weighted density'
fxdens[:,:,:] = T
newnc.close()