# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 12:19:35 2017

@author: np838619
"""

import pylab as pl
from netCDF4 import Dataset
#from kernelfuncs import SaveFields

def SavetoNC(densfield,fluxfield,shed,catch,cat):
    # which category of trajectories? all, CAT 1 or CAT 2?
    if cat == 'all':
        field_ex = ' all'
    elif cat == '1':
        field_ex = ' CAT I'
    elif cat == '2':
        field_ex = ' CAT II'
    
    clusdir = '/net/glusterfs/scenario/users/np838619/traj/meandens/'
    ncfile = Dataset(clusdir+shed+'/ordens_'+shed+'12.nc','a')
    
    D = ncfile.createVariable('origin density ' + shed + catch + field_ex,pl.float64,
                                                          ('time','lat','lon'))
    D.units = 'trajectories per steradian'
    D.standard_name = 'trajectory origin density ' + shed + catch + field_ex
    D[:,:,:] = densfield
    
    F = ncfile.createVariable('flux density ' + shed + catch + field_ex,pl.float64,
                                                              ('time','lat','lon'))
    F.units = 'kg/s per steradian'
    F.standard_name = 'moisture flux density ' + shed + catch + field_ex
    F[:,:,:] = fluxfield
    
    ncfile.close()

clusdir = '/net/glusterfs_essc/scenario/users/np838619/traj/'

febrel = pl.zeros([5,63,256,512])
fd_all_all = pl.zeros([5,63,256,512]); ff_all_all = pl.zeros_like(fd_all_all)
fd_all_c1 = pl.zeros_like(fd_all_all); ff_all_c1 = pl.zeros_like(fd_all_all)
fd_all_c2 = pl.zeros_like(fd_all_all); ff_all_c2 = pl.zeros_like(fd_all_all)
fd_atl_all = pl.zeros([5,63,256,512]); ff_atl_all = pl.zeros_like(fd_all_all)
fd_atl_c1 = pl.zeros_like(fd_all_all); ff_atl_c1 = pl.zeros_like(fd_all_all)
fd_atl_c2 = pl.zeros_like(fd_all_all); ff_atl_c2 = pl.zeros_like(fd_all_all)
fd_ind_all = pl.zeros([5,63,256,512]); ff_ind_all = pl.zeros_like(fd_all_all)
fd_ind_c1 = pl.zeros_like(fd_all_all); ff_ind_c1 = pl.zeros_like(fd_all_all)
fd_ind_c2 = pl.zeros_like(fd_all_all); ff_ind_c2 = pl.zeros_like(fd_all_all)
fd_pac_all = pl.zeros([5,63,256,512]); ff_pac_all = pl.zeros_like(fd_all_all)
fd_pac_c1 = pl.zeros_like(fd_all_all); ff_pac_c1 = pl.zeros_like(fd_all_all)
fd_pac_c2 = pl.zeros_like(fd_all_all); ff_pac_c2 = pl.zeros_like(fd_all_all)
fd_arc_all = pl.zeros([5,63,256,512]); ff_arc_all = pl.zeros_like(fd_all_all)
fd_arc_c1 = pl.zeros_like(fd_all_all); ff_arc_c1 = pl.zeros_like(fd_all_all)
fd_arc_c2 = pl.zeros_like(fd_all_all); ff_arc_c2 = pl.zeros_like(fd_all_all)
fd_sou_all = pl.zeros([5,63,256,512]); ff_sou_all = pl.zeros_like(fd_all_all)
fd_sou_c1 = pl.zeros_like(fd_all_all); ff_sou_c1 = pl.zeros_like(fd_all_all)
fd_sou_c2 = pl.zeros_like(fd_all_all); ff_sou_c2 = pl.zeros_like(fd_all_all)

decrel = pl.zeros([62,256,512])

years = ['2010','2011','2012','2013','2014']
shed = 'soi'

ntraj10 = pl.genfromtxt(clusdir+'2010/amr/NTRAJ/amr_all201012_all.csv')
ntraj11 = pl.genfromtxt(clusdir+'2011/amr/NTRAJ/amr_all201112_all.csv')
ntraj12 = pl.genfromtxt(clusdir+'2012/amr/NTRAJ/amr_all201212_all.csv')
ntraj13 = pl.genfromtxt(clusdir+'2013/amr/NTRAJ/amr_all201312_all.csv')
ntraj14 = pl.genfromtxt(clusdir+'2014/amr/NTRAJ/amr_all201412_all.csv')
N = pl.sum(ntraj10) + pl.sum(ntraj11) + pl.sum(ntraj12) + pl.sum(ntraj13) + pl.sum(ntraj14)

for i in range(len(years)):
    febfile = Dataset(clusdir+years[i]+'/'+shed+'/density/ordens_'+shed+years[i]+'12.nc','r')
    den_all_all = febfile.variables['origin density '+shed+'all all'][:]
    flx_all_all = febfile.variables['flux density '+shed+'all all'][:]
    den_all_c1 = febfile.variables['origin density '+shed+'all CAT I'][:]
    flx_all_c1 = febfile.variables['flux density '+shed+'all CAT I'][:]
    den_all_c2 = febfile.variables['origin density '+shed+'all CAT II'][:]
    flx_all_c2 = febfile.variables['flux density '+shed+'all CAT II'][:]
    den_atl_all = febfile.variables['origin density '+shed+'atl all'][:]
    flx_atl_all = febfile.variables['flux density '+shed+'atl all'][:]
    den_atl_c1 = febfile.variables['origin density '+shed+'atl CAT I'][:]
    flx_atl_c1 = febfile.variables['flux density '+shed+'atl CAT I'][:]
    den_atl_c2 = febfile.variables['origin density '+shed+'atl CAT II'][:]
    flx_atl_c2 = febfile.variables['flux density '+shed+'atl CAT II'][:]
    den_ind_all = febfile.variables['origin density '+shed+'ind all'][:]
    flx_ind_all = febfile.variables['flux density '+shed+'ind all'][:]
    den_ind_c1 = febfile.variables['origin density '+shed+'ind CAT I'][:]
    flx_ind_c1 = febfile.variables['flux density '+shed+'ind CAT I'][:]
    den_ind_c2 = febfile.variables['origin density '+shed+'ind CAT II'][:]
    flx_ind_c2 = febfile.variables['flux density '+shed+'ind CAT II'][:]
    den_pac_all = febfile.variables['origin density '+shed+'pac all'][:]
    flx_pac_all = febfile.variables['flux density '+shed+'pac all'][:]
    den_pac_c1 = febfile.variables['origin density '+shed+'pac CAT I'][:]
    flx_pac_c1 = febfile.variables['flux density '+shed+'pac CAT I'][:]
    den_pac_c2 = febfile.variables['origin density '+shed+'pac CAT II'][:]
    flx_pac_c2 = febfile.variables['flux density '+shed+'pac CAT II'][:]
    den_arc_all = febfile.variables['origin density '+shed+'arc all'][:]
    flx_arc_all = febfile.variables['flux density '+shed+'arc all'][:]
    den_arc_c1 = febfile.variables['origin density '+shed+'arc CAT I'][:]
    flx_arc_c1 = febfile.variables['flux density '+shed+'arc CAT I'][:]
    den_arc_c2 = febfile.variables['origin density '+shed+'arc CAT II'][:]
    flx_arc_c2 = febfile.variables['flux density '+shed+'arc CAT II'][:]
    den_sou_all = febfile.variables['origin density '+shed+'sou all'][:]
    flx_sou_all = febfile.variables['flux density '+shed+'sou all'][:]
    den_sou_c1 = febfile.variables['origin density '+shed+'sou CAT I'][:]
    flx_sou_c1 = febfile.variables['flux density '+shed+'sou CAT I'][:]
    den_sou_c2 = febfile.variables['origin density '+shed+'sou CAT II'][:]
    flx_sou_c2 = febfile.variables['flux density '+shed+'sou CAT II'][:]
    febfile.close()
    if years[i] == '2014':
        fd_all_all[i,:,:,:] = den_all_all[:,:,:]
        ff_all_all[i,:,:,:] = flx_all_all[:,:,:]
        fd_all_c1[i,:,:,:] = den_all_c1[:,:,:]
        ff_all_c1[i,:,:,:] = flx_all_c1[:,:,:]
        fd_all_c2[i,:,:,:] = den_all_c2[:,:,:]
        ff_all_c2[i,:,:,:] = flx_all_c2[:,:,:]
        fd_atl_all[i,:,:,:] = den_atl_all[:,:,:]
        ff_atl_all[i,:,:,:] = flx_atl_all[:,:,:]
        fd_atl_c1[i,:,:,:] = den_atl_c1[:,:,:]
        ff_atl_c1[i,:,:,:] = flx_atl_c1[:,:,:]
        fd_atl_c2[i,:,:,:] = den_atl_c2[:,:,:]
        ff_atl_c2[i,:,:,:] = flx_atl_c2[:,:,:]
        fd_pac_all[i,:,:,:] = den_pac_all[:,:,:]
        ff_pac_all[i,:,:,:] = flx_pac_all[:,:,:]
        fd_pac_c1[i,:,:,:] = den_pac_c1[:,:,:]
        ff_pac_c1[i,:,:,:] = flx_pac_c1[:,:,:]
        fd_pac_c2[i,:,:,:] = den_pac_c2[:,:,:]
        ff_pac_c2[i,:,:,:] = flx_pac_c2[:,:,:]
        fd_arc_all[i,:,:,:] = den_arc_all[:,:,:]
        ff_arc_all[i,:,:,:] = flx_arc_all[:,:,:]
        fd_arc_c1[i,:,:,:] = den_arc_c1[:,:,:]
        ff_arc_c1[i,:,:,:] = flx_arc_c1[:,:,:]
        fd_arc_c2[i,:,:,:] = den_arc_c2[:,:,:]
        ff_arc_c2[i,:,:,:] = flx_arc_c2[:,:,:]
        fd_sou_all[i,:,:,:] = den_sou_all[:,:,:]
        ff_sou_all[i,:,:,:] = flx_sou_all[:,:,:]
        fd_sou_c1[i,:,:,:] = den_sou_c1[:,:,:]
        ff_sou_c1[i,:,:,:] = flx_sou_c1[:,:,:]
        fd_sou_c2[i,:,:,:] = den_sou_c2[:,:,:]
        ff_sou_c2[i,:,:,:] = flx_sou_c2[:,:,:]
    else:
        fd_all_all[i,:-1,:,:] = den_all_all[:,:,:]
        ff_all_all[i,:-1,:,:] = flx_all_all[:,:,:]
        fd_all_c1[i,:-1,:,:] = den_all_c1[:,:,:]
        ff_all_c1[i,:-1,:,:] = flx_all_c1[:,:,:]
        fd_all_c2[i,:-1,:,:] = den_all_c2[:,:,:]
        ff_all_c2[i,:-1,:,:] = flx_all_c2[:,:,:]
        fd_atl_all[i,:-1,:,:] = den_atl_all[:,:,:]
        ff_atl_all[i,:-1,:,:] = flx_atl_all[:,:,:]
        fd_atl_c1[i,:-1,:,:] = den_atl_c1[:,:,:]
        ff_atl_c1[i,:-1,:,:] = flx_atl_c1[:,:,:]
        fd_atl_c2[i,:-1,:,:] = den_atl_c2[:,:,:]
        ff_atl_c2[i,:-1,:,:] = flx_atl_c2[:,:,:]
        fd_pac_all[i,:-1,:,:] = den_pac_all[:,:,:]
        ff_pac_all[i,:-1,:,:] = flx_pac_all[:,:,:]
        fd_pac_c1[i,:-1,:,:] = den_pac_c1[:,:,:]
        ff_pac_c1[i,:-1,:,:] = flx_pac_c1[:,:,:]
        fd_pac_c2[i,:-1,:,:] = den_pac_c2[:,:,:]
        ff_pac_c2[i,:-1,:,:] = flx_pac_c2[:,:,:]
        fd_arc_all[i,:-1,:,:] = den_arc_all[:,:,:]
        ff_arc_all[i,:-1,:,:] = flx_arc_all[:,:,:]
        fd_arc_c1[i,:-1,:,:] = den_arc_c1[:,:,:]
        ff_arc_c1[i,:-1,:,:] = flx_arc_c1[:,:,:]
        fd_arc_c2[i,:-1,:,:] = den_arc_c2[:,:,:]
        ff_arc_c2[i,:-1,:,:] = flx_arc_c2[:,:,:]
        fd_sou_all[i,:-1,:,:] = den_sou_all[:,:,:]
        ff_sou_all[i,:-1,:,:] = flx_sou_all[:,:,:]
        fd_sou_c1[i,:-1,:,:] = den_sou_c1[:,:,:]
        ff_sou_c1[i,:-1,:,:] = flx_sou_c1[:,:,:]
        fd_sou_c2[i,:-1,:,:] = den_sou_c2[:,:,:]
        ff_sou_c2[i,:-1,:,:] = flx_sou_c2[:,:,:]

fd_all_all = pl.nansum(fd_all_all,axis=0); ff_all_all = pl.nansum(ff_all_all,axis=0)
fd_all_c1 = pl.nansum(fd_all_c1,axis=0); ff_all_c1 = pl.nansum(ff_all_c1,axis=0)
fd_all_c2 = pl.nansum(fd_all_c2,axis=0); ff_all_c2 = pl.nansum(ff_all_c2,axis=0)
fd_atl_all = pl.nansum(fd_atl_all,axis=0); ff_atl_all = pl.nansum(ff_atl_all,axis=0)
fd_atl_c1 = pl.nansum(fd_atl_c1,axis=0); ff_atl_c1 = pl.nansum(ff_atl_c1,axis=0)
fd_atl_c2 = pl.nansum(fd_atl_c2,axis=0); ff_atl_c2 = pl.nansum(ff_atl_c2,axis=0)
fd_ind_all = pl.nansum(fd_ind_all,axis=0); ff_ind_all = pl.nansum(ff_ind_all,axis=0)
fd_ind_c1 = pl.nansum(fd_ind_c1,axis=0); ff_ind_c1 = pl.nansum(ff_ind_c1,axis=0)
fd_ind_c2 = pl.nansum(fd_ind_c2,axis=0); ff_ind_c2 = pl.nansum(ff_ind_c2,axis=0)
fd_pac_all = pl.nansum(fd_pac_all,axis=0); ff_pac_all = pl.nansum(ff_pac_all,axis=0)
fd_pac_c1 = pl.nansum(fd_pac_c1,axis=0); ff_pac_c1 = pl.nansum(ff_pac_c1,axis=0)
fd_pac_c2 = pl.nansum(fd_pac_c2,axis=0); ff_pac_c2 = pl.nansum(ff_pac_c2,axis=0)
fd_arc_all = pl.nansum(fd_arc_all,axis=0); ff_arc_all = pl.nansum(ff_arc_all,axis=0)
fd_arc_c1 = pl.nansum(fd_arc_c1,axis=0); ff_arc_c1 = pl.nansum(ff_arc_c1,axis=0)
fd_arc_c2 = pl.nansum(fd_arc_c2,axis=0); ff_arc_c2 = pl.nansum(ff_arc_c2,axis=0)
fd_sou_all = pl.nansum(fd_sou_all,axis=0); ff_sou_all = pl.nansum(ff_sou_all,axis=0)
fd_sou_c1 = pl.nansum(fd_sou_c1,axis=0); ff_sou_c1 = pl.nansum(ff_sou_c1,axis=0)
fd_sou_c2 = pl.nansum(fd_sou_c2,axis=0); ff_sou_c2 = pl.nansum(ff_sou_c2,axis=0)


ncfile = Dataset('/home/np838619/Watershed/quv_annmean.nc','r')
eralat = ncfile.variables['lat'][:]
eralon = ncfile.variables['lon'][:]
ncfile.close()

newnc = Dataset(clusdir+'meandens/'+shed+'/ordens_'+shed+'12.nc','w')

lat_dim = newnc.createDimension('lat',eralat.size)
lat_in = newnc.createVariable('lat',pl.float64,('lat',))
lat_in.units = 'degrees_north'
lat_in.long_name = 'latitude'
lat_in[:] = eralat

lon_dim = newnc.createDimension('lon',eralon.size)
lon_in = newnc.createVariable('lon',pl.float64,('lon',))
lon_in.units = 'degrees_east'
lon_in.long_name = 'longitude'
lon_in[:] = eralon

time_dim = newnc.createDimension('time',fd_all_all.shape[0])
time = newnc.createVariable('time',pl.float64,('time',))
time.units = 'hours'
time.long_name = 'release date/time'

newnc.close()

SavetoNC(fd_all_all,ff_all_all,shed,'all','all')
SavetoNC(fd_all_c1,ff_all_c1,shed,'all','1')
SavetoNC(fd_all_c2,ff_all_c2,shed,'all','2')
SavetoNC(fd_atl_all,ff_atl_all,shed,'atl','all')
SavetoNC(fd_atl_c1,ff_atl_c1,shed,'atl','1')
SavetoNC(fd_atl_c2,ff_atl_c2,shed,'atl','2')
SavetoNC(fd_ind_all,ff_ind_all,shed,'ind','all')
SavetoNC(fd_ind_c1,ff_ind_c1,shed,'ind','1')
SavetoNC(fd_ind_c2,ff_ind_c2,shed,'ind','2')
SavetoNC(fd_pac_all,ff_pac_all,shed,'pac','all')
SavetoNC(fd_pac_c1,ff_pac_c1,shed,'pac','1')
SavetoNC(fd_pac_c2,ff_pac_c2,shed,'pac','2')
SavetoNC(fd_arc_all,ff_arc_all,shed,'arc','all')
SavetoNC(fd_arc_c1,ff_arc_c1,shed,'arc','1')
SavetoNC(fd_arc_c2,ff_arc_c2,shed,'arc','2')
SavetoNC(fd_sou_all,ff_sou_all,shed,'sou','all')
SavetoNC(fd_sou_c1,ff_sou_c1,shed,'sou','1')
SavetoNC(fd_sou_c2,ff_sou_c2,shed,'sou','2')